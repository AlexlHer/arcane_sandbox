// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// computeVelocity() and computePhi() and end of startInit() are based on
// AMReX example 06_Advection_Amr available here :
// https://github.com/atmyers/ecp-tutorials/blob/main/06_Advection_Amr
// _cutDim part is based on Berger-Rigoutsos algo.

#include "SayHelloModule.h"

#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/core/IMesh.h>
#include <arcane/utils/StringBuilder.h>
#include <arcane/core/IParallelMng.h>
#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/Vector2.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/Array3View.h>

#include "arcane/core/ICartesianMeshGenerationInfo.h"
#include "arcane/cartesianmesh/ICartesianMeshAMRPatchMng.h"
#include "arcane/cartesianmesh/CartesianMeshUtils.h"

#include <arcane/cartesianmesh/CartesianPatch.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/FaceDirectionMng.h>

using namespace Arcane;


constexpr Integer MIN_SIZE = 4;
constexpr Integer TARGET_SIZE = 8;
constexpr Real TARGET_SIZE_WEIGHT_IN_EFFICACITY = 0.5;
constexpr Integer MAX_NB_CUT = 5;
constexpr Real TARGET_EFFICACITY = 0.90;

struct Patch
{
  Patch()
  : m_level(-1)
  , m_min(-1)
  , m_max(-1)
  {}

  std::pair<Patch, Patch> cut(Integer dim, Integer cut_point) const
  {
    Patch p0;
    p0.m_level = m_level;
    p0.m_min = m_min;
    p0.m_max = m_max;

    Patch p1;
    p1.m_level = m_level;
    p1.m_min = m_min;
    p1.m_max = m_max;

    if (dim == MD_DirX) {
      p0.m_max.x = cut_point;
      p1.m_min.x = cut_point;

      if (p0.m_max.x <= p0.m_min.x) {
        ARCANE_FATAL("Bad cut_point ori");
      }
      if (p1.m_max.x <= p1.m_min.x) {
        ARCANE_FATAL("Bad cut_point p1");
      }
    }

    else if (dim == MD_DirY) {
      p0.m_max.y = cut_point;
      p1.m_min.y = cut_point;
      if (p0.m_max.y <= p0.m_min.y) {
        ARCANE_FATAL("Bad cut_point ori");
      }
      if (p1.m_max.y <= p1.m_min.y) {
        ARCANE_FATAL("Bad cut_point p1");
      }
    }
    else {
      p0.m_max.z = cut_point;
      p1.m_min.z = cut_point;
      if (p0.m_max.z <= p0.m_min.z) {
        ARCANE_FATAL("Bad cut_point ori");
      }
      if (p1.m_max.z <= p1.m_min.z) {
        ARCANE_FATAL("Bad cut_point p1");
      }
    }
    return {p0, p1};
  }

  Int32x3 length() const
  {
    return m_max - m_min;
  }

  Int32x3 minWithMargin(Integer level) const
  {
    if (level == m_level) {
      return m_min - 1;
    }
    if (level == m_level-1) {
      return (m_min - 1) / 2;
    }
    ARCANE_FATAL("Pas utile");
  }

  Int32x3 maxWithMargin(Integer level) const
  {
    if (level == m_level) {
      return m_max + 1;
    }
    if (level == m_level-1) {
      Int32x3 max = m_max + 1;
      return { static_cast<Int32>(std::ceil(max.x / 2.)), static_cast<Int32>(std::ceil(max.y / 2.)), static_cast<Int32>(std::ceil(max.z / 2.)) };
    }
    ARCANE_FATAL("Pas utile");
  }

  Int32x3 min(Integer level) const
  {
    if (level == m_level) {
      return m_min;
    }
    if (level == m_level+1) {
      return m_min * 2;
    }
    if (level == m_level-1) {
      return m_min / 2;
    }
    if (level < m_level) {
      Int32 dif = static_cast<Int32>(math::pow(2., static_cast<Real>(m_level - level)));
      return m_min / dif;
    }

    Int32 dif = static_cast<Int32>(math::pow(2., static_cast<Real>(level - m_level)));
    return m_min * dif;
  }

  Int32x3 max(Integer level) const
  {
    if (level == m_level) {
      return m_max;
    }
    if (level == m_level+1) {
      return m_max * 2;
    }
    if (level == m_level-1) {
      return { static_cast<Int32>(std::ceil(m_max.x / 2.)), static_cast<Int32>(std::ceil(m_max.y / 2.)), static_cast<Int32>(std::ceil(m_max.z / 2.)) };
    }
    if (level < m_level) {
      Int32 dif = static_cast<Int32>(math::pow(2., static_cast<Real>(level - m_level)));
      return { static_cast<Int32>(std::ceil(m_max.x / dif)), static_cast<Int32>(std::ceil(m_max.y / dif)), static_cast<Int32>(std::ceil(m_max.z / dif)) };
    }
    Int32 dif = static_cast<Int32>(math::pow(2., static_cast<Real>(level - m_level)));
    return m_max * dif;
  }

  bool isIn(Integer x, Integer y) const
  {
    return x >= m_min.x && x < m_max.x && y >= m_min.y && y < m_max.y;
  }

  bool isInWithMargin(Integer level, Integer x, Integer y) const
  {
    Int32x3 level_min = minWithMargin(level);
    Int32x3 level_max = maxWithMargin(level);
    return x >= level_min.x && x < level_max.x && y >= level_min.y && y < level_max.y;
  }

  Integer m_level;
  Int32x3 m_min;
  Int32x3 m_max;
};




struct LevelPatches
{
  LevelPatches(Integer max_level, IMesh* mesh)
  : m_max_level(max_level)
  , m_patches(max_level+1)
  , m_mesh(mesh)
  {}

  Integer m_max_level;
  UniqueArray<UniqueArray<Patch>> m_patches;
  IMesh* m_mesh;
};

/*
struct PatchCells
{
  PatchCells(IMesh* mesh, Integer level, Int32x3 nb_cells)
  : m_level(level)
  , m_nb_cells(nb_cells)
  , m_mesh(mesh)
  {}

  void computeBoolCells()
  {
    m_cells.resize(m_nb_cells.x * m_nb_cells.y * m_nb_cells.z, false);
    Array3View<bool> av(m_cells.data(), m_nb_cells.x, m_nb_cells.y, m_nb_cells.z);

    ENUMERATE_(Cell, icell, m_mesh->ownCells()) {
      if (icell->level() == m_level && icell->hasFlags(ItemFlags::II_Refine)) {
        Integer pos_x = icell->uniqueId().asInt32() % m_nb_cells.x; // Pas bon level != 0
        Integer pos_y = icell->uniqueId().asInt32() / m_nb_cells.x;
        av(pos_x, pos_y, 0) = true;
      }
    }
  }

  void computePatchCells(LevelPatches* all_patches)
  {
    UniqueArray<Patch>& patches = all_patches->m_patches[m_level];
    Array3View<bool> av(m_cells.data(), m_nb_cells.x, m_nb_cells.y, m_nb_cells.z);

    for (auto patch : patches) {
      Int32x3 min = patch.m_min;
      Int32x3 max = patch.m_max;

      for (Integer i = min.x; i < max.x; ++i) {
        for (Integer j = min.y; j < max.y; ++j) {
          for (Integer k = min.z; k < max.z; ++k) {
            av(i, j, k) = true;
          }
        }
      }
    }
  }

  void addCellsFromLevel(LevelPatches* all_patches)
  {
    if (m_level == all_patches->m_max_level) return;

    UniqueArray<Patch>& patches = all_patches->m_patches[m_level+1];
    Array3View<bool> av(m_cells.data(), m_nb_cells.x, m_nb_cells.y, m_nb_cells.z);

    for (auto patch : patches) {
      Int32x3 min = patch.min(m_level);
      Int32x3 max = patch.max(m_level);

      if (min.x > 0) min.x--;
      if (min.y > 0) min.y--;
      if (min.z > 0) min.z--;

      if (max.x < m_nb_cells.x-1) max.x++;
      if (max.y < m_nb_cells.y-1) max.y++;
      if (max.z < m_nb_cells.z-1) max.z++;

      for (Integer i = min.x; i < max.x; ++i) {
        for (Integer j = min.y; j < max.y; ++j) {
          for (Integer k = min.z; k < max.z; ++k) {
            av(i, j, k) = true;
          }
        }
      }
    }
  }

  UniqueArray<bool> m_cells;
  Integer m_level;
  Int32x3 m_nb_cells;
  IMesh* m_mesh;
};
*/

struct Signature
{
  Signature()
  : m_is_null(true)
  , m_mesh(nullptr)
  , m_nb_cut(0)
  , m_stop_cut(false)
  , m_numbering(nullptr)
  , m_have_cells(false)
  , m_is_computed(false)
  , m_all_patches(nullptr)
  {}

  Signature(Patch patch, ICartesianMeshNumberingMng* numbering, LevelPatches* all_patches)
  : m_is_null(false)
  , m_patch(patch)
  , m_mesh(all_patches->m_mesh)
  , m_nb_cut(0)
  , m_stop_cut(false)
  , m_numbering(numbering)
  , m_have_cells(false)
  , m_is_computed(false)
  , m_sig_x(patch.m_max.x - patch.m_min.x, 0)
  , m_sig_y(patch.m_max.y - patch.m_min.y, 0)
  , m_all_patches(all_patches)
  {}

  Signature(Patch patch, ICartesianMeshNumberingMng* numbering, LevelPatches* all_patches, Integer nb_cut)
  : m_is_null(false)
  , m_patch(patch)
  , m_mesh(all_patches->m_mesh)
  , m_nb_cut(nb_cut)
  , m_stop_cut(false)
  , m_numbering(numbering)
  , m_have_cells(false)
  , m_is_computed(false)
  , m_sig_x(patch.m_max.x - patch.m_min.x, 0)
  , m_sig_y(patch.m_max.y - patch.m_min.y, 0)
  , m_all_patches(all_patches)
  {}

  void compress()
  {
    if (!m_have_cells) {
      return;
    }

    Integer reduce_x_min = 0;
    if (m_sig_x[0] == 0) {
      for (; reduce_x_min < m_sig_x.size(); ++reduce_x_min) {
        if (m_sig_x[reduce_x_min] != 0) {
          break;
        }
      }
      if (reduce_x_min == m_sig_x.size()) {
        ARCANE_FATAL("aaa");
      }
    }
    Integer reduce_y_min = 0;
    if (m_sig_y[0] == 0) {
      for (; reduce_y_min < m_sig_y.size(); ++reduce_y_min) {
        if (m_sig_y[reduce_y_min] != 0) {
          break;
        }
      }
      if (reduce_y_min == m_sig_y.size()) {
        ARCANE_FATAL("bbb");
      }
    }

    Integer reduce_x_max = m_sig_x.size()-1;
    if (m_sig_x[reduce_x_max] == 0) {
      for (; reduce_x_max >= 0; --reduce_x_max) {
        if (m_sig_x[reduce_x_max] != 0) {
          break;
        }
      }
      if (reduce_x_max < reduce_x_min) {
        ARCANE_FATAL("ccc");
      }
    }
    Integer reduce_y_max = m_sig_y.size()-1;
    if (m_sig_y[reduce_y_max] == 0) {
      for (; reduce_y_max >= 0; --reduce_y_max) {
        if (m_sig_y[reduce_y_max] != 0) {
          break;
        }
      }
      if (reduce_y_max < reduce_y_min) {
        ARCANE_FATAL("ddd");
      }
    }

    if (reduce_x_min != 0 || reduce_x_max != m_sig_x.size()-1) {
      reduce_x_max++;
      UniqueArray tmp = m_sig_x.subView(reduce_x_min, reduce_x_max - reduce_x_min);
      m_sig_x = tmp;
      m_patch.m_min.x += reduce_x_min;
      m_patch.m_max.x = m_patch.m_min.x + (reduce_x_max - reduce_x_min);
    }
    if (reduce_y_min != 0 || reduce_y_max != m_sig_y.size()-1) {
      reduce_y_max++;
      UniqueArray tmp = m_sig_y.subView(reduce_y_min, reduce_y_max - reduce_y_min);
      m_sig_y = tmp;
      m_patch.m_min.y += reduce_y_min;
      m_patch.m_max.y = m_patch.m_min.y + (reduce_y_max - reduce_y_min);
    }
  }

  void fillSig()
  {
    m_sig_x.fill(0);
    m_sig_y.fill(0);
    ENUMERATE_ (Cell, icell, m_mesh->ownCells()) {
      if (!icell->hasFlags(ItemFlags::II_Refine) || icell->level() != m_patch.m_level) {
        continue;
      }

      Integer pos_x = m_numbering->cellUniqueIdToCoordX(*icell);
      Integer pos_y = m_numbering->cellUniqueIdToCoordY(*icell);

      if (pos_x < m_patch.m_min.x || pos_x >= m_patch.m_max.x || pos_y < m_patch.m_min.y || pos_y >= m_patch.m_max.y ) {
        continue;
      }
      m_have_cells = true;
      m_sig_x[pos_x - m_patch.m_min.x]++;
      m_sig_y[pos_y - m_patch.m_min.y]++;
    }

    if (m_all_patches->m_max_level > m_patch.m_level) {
      Int32x3 min = m_patch.m_min;
      Integer nb_proc = m_mesh->parallelMng()->commSize();
      Integer my_proc = m_mesh->parallelMng()->commRank();

      Integer base = m_sig_x.size() / nb_proc;
      Integer reste = m_sig_x.size() % nb_proc;
      Integer size = base + (my_proc < reste ? 1 : 0);
      Integer begin = my_proc * base + std::min(my_proc, reste);
      Integer end = begin + size;

      for (Integer j = 0; j < m_sig_y.size(); ++j) {
        Integer pos_y = min.y + j;
        for (Integer i = begin; i < end; ++i) {
          Integer pos_x = min.x + i;
          for (auto elem : m_all_patches->m_patches[m_patch.m_level+1]) {
            if (elem.isInWithMargin(m_patch.m_level, pos_x, pos_y)) {
              m_have_cells = true;
              m_sig_x[i]++;
              m_sig_y[j]++;
              break;
            }
          }
        }
      }
    }


    m_mesh->parallelMng()->reduce(MessagePassing::ReduceSum, m_sig_x);
    m_mesh->parallelMng()->reduce(MessagePassing::ReduceSum, m_sig_y);
    m_have_cells = m_mesh->parallelMng()->reduce(MessagePassing::ReduceMax, m_have_cells);
  }

  bool isValid() const
  {
    if (m_is_null) {
      return false;
    }
    if (m_sig_x.size() < MIN_SIZE || m_sig_y.size() < MIN_SIZE) {
      return false;
    }
    return true;
  }

  bool canBeCut() const
  {
    m_mesh->traceMng()->info() << "canBeCut() -- m_sig_x.size : " << m_sig_x.size()
      << " -- m_sig_y.size : " << m_sig_y.size()
      << " -- min = " << m_patch.m_min
      << " -- max = " << m_patch.m_max
      << " -- length = " << m_patch.length()
      << " -- isValid : " << isValid()
      << " -- efficacity : " << efficacity() << " / " << TARGET_EFFICACITY
      << " -- m_nb_cut : " << m_nb_cut << " / " << MAX_NB_CUT
      << " -- m_stop_cut : " << m_stop_cut
    ;

    if (!isValid()) {
      return false;
    }

    if (m_stop_cut) {
      return false;
    }

    if (efficacity() > TARGET_EFFICACITY) {
      return false;
    }
    if (MAX_NB_CUT != -1 && m_nb_cut >= MAX_NB_CUT) {
      return false;
    }
    return true;
  }

  void compute()
  {
    m_mesh->traceMng()->info() << "Compute() -- Patch before compute : min = " << m_patch.m_min << " -- max = " << m_patch.m_max << " -- length = " << m_patch.length();
    fillSig();
    //m_mesh->traceMng()->info() << "Compute() -- Signature : x = " << m_sig_x << " -- y = " << m_sig_y ;
    compress();
    m_mesh->traceMng()->info() << "Compute() -- Compress : min = " << m_patch.m_min << " -- max = " << m_patch.m_max << " -- x = " << m_sig_x << " -- y = " << m_sig_y ;
    m_mesh->traceMng()->info() << "Compute() -- Patch computed :       min = " << m_patch.m_min << " -- max = " << m_patch.m_max << " -- length = " << m_patch.length();

    m_is_computed = true;
  }

  Real efficacity() const
  {
    if (!m_is_computed) {
      // Sans compression, pas terrible.
      m_mesh->traceMng()->warning() << "Need to be computed";
    }
    Integer sum = 0;
    for (const Integer elem : m_sig_x) {
      sum += elem;
    }

    Real eff = static_cast<Real>(sum) / (m_sig_x.size() * m_sig_y.size());

    if constexpr (TARGET_SIZE == -1 || TARGET_SIZE_WEIGHT_IN_EFFICACITY == 0) {
      return eff;
    }

    Real eff_xy = 0;
    if (m_sig_x.size() <= TARGET_SIZE) {
      eff_xy = static_cast<Real>(m_sig_x.size()) / TARGET_SIZE;
    }
    else if (m_sig_x.size() < TARGET_SIZE*2) {
      Real size_x = math::abs(m_sig_x.size() - TARGET_SIZE*2);
      eff_xy = size_x / TARGET_SIZE;
    }

    if (m_sig_y.size() <= TARGET_SIZE) {
      eff_xy = (eff_xy + (static_cast<Real>(m_sig_y.size()) / TARGET_SIZE)) / 2;
    }
    else if (m_sig_y.size() < TARGET_SIZE*2) {
      Real size_y = math::abs(m_sig_y.size() - TARGET_SIZE*2);
      eff_xy = (eff_xy + (size_y / TARGET_SIZE)) / 2;
    }
    else {
      eff_xy /= 2;
    }
    eff_xy *= TARGET_SIZE_WEIGHT_IN_EFFICACITY;

    return (eff+eff_xy)/(1+TARGET_SIZE_WEIGHT_IN_EFFICACITY);
  }

  std::pair<Signature, Signature> cut(Integer dim, Integer cut_point) const
  {
    auto [fst, snd] = m_patch.cut(dim, cut_point);
    return {Signature(fst, m_numbering, m_all_patches, m_nb_cut+1), Signature(snd, m_numbering, m_all_patches, m_nb_cut+1)};
  }

  bool isIn(Integer x, Integer y) const
  {
    return m_patch.isIn(x, y);
  }

  bool m_is_null;
  Patch m_patch;
  IMesh* m_mesh;
  Integer m_nb_cut;
  bool m_stop_cut;

  ICartesianMeshNumberingMng* m_numbering;


  bool m_have_cells;
  bool m_is_computed;

  UniqueArray<Integer> m_sig_x;
  UniqueArray<Integer> m_sig_y;

  LevelPatches* m_all_patches;
};



struct CutPatch
{
  static Integer _cutDim(ConstArrayView<Integer> sig)
  {
    if (sig.size() < MIN_SIZE*2) {
      return -1;
    }

    Integer cut_point = -1;
    Integer mid = sig.size() / 2;

    {
      for (Integer i = 0; i < sig.size(); ++i) {
        if (sig[i] == 0) {
          cut_point = i;
          break;
        }
      }

      if (cut_point == 0) {
        ARCANE_FATAL("Call compress before");
      }
      if (cut_point != -1 && cut_point >= MIN_SIZE && sig.size() - cut_point >= MIN_SIZE) {
        return cut_point;
      }
    }

    {
      UniqueArray<Integer> dsec_sig(sig.size(), 0);

      Integer max = -1;
      for (Integer i = 1; i < dsec_sig.size()-1; ++i) {
        dsec_sig[i] = sig[i+1] - 2 * sig[i] + sig[i-1];
        Integer dif = math::abs(dsec_sig[i-1] - dsec_sig[i]);
        if (dif > max) {
          cut_point = i;
          max = dif;
        }
        else if (dif == max && math::abs(cut_point - mid) > math::abs(i - mid)) {
          cut_point = i;
        }
      }

      if (cut_point != -1 && cut_point >= MIN_SIZE && sig.size() - cut_point >= MIN_SIZE) {
        return cut_point;
      }
    }

    {
      cut_point = mid;

      if (cut_point != -1 && cut_point >= MIN_SIZE && sig.size() - cut_point >= MIN_SIZE) {
        return cut_point;
      }
    }

    return -1;
  }

  static std::pair<Signature, Signature> cut(const Signature& sig)
  {
    Integer cut_point_x = _cutDim(sig.m_sig_x);
    Integer cut_point_y = _cutDim(sig.m_sig_y);


    if (cut_point_x == -1 && cut_point_y == -1) {
      return {};
    }
    if (cut_point_x != -1) {
      cut_point_x += sig.m_patch.m_min.x;
    }
    if (cut_point_y != -1) {
      cut_point_y += sig.m_patch.m_min.y;
    }

    if (cut_point_x != -1 && cut_point_y != -1) {
      Real x_efficacity = 0;
      auto [fst_x, snd_x] = sig.cut(MD_DirX, cut_point_x);
      {
        sig.m_mesh->traceMng()->info() << "Cut() -- Compute X -- Cut Point : " << cut_point_x;

        fst_x.compute();
        snd_x.compute();
        if (fst_x.isValid() && snd_x.isValid()) {

          sig.m_mesh->traceMng()->info() << "Cut() -- X.fst_x"
            << " -- min = " << fst_x.m_patch.m_min
            << " -- max = " << fst_x.m_patch.m_max
            << " -- efficacity : " << fst_x.efficacity()
          ;
          sig.m_mesh->traceMng()->info() << "Cut() -- X.snd_x"
            << " -- min = " << snd_x.m_patch.m_min
            << " -- max = " << snd_x.m_patch.m_max
            << " -- efficacity : " << snd_x.efficacity()
          ;

          x_efficacity = (fst_x.efficacity() + snd_x.efficacity()) / 2;
          sig.m_mesh->traceMng()->info() << "Cut() -- efficacity X : " << x_efficacity;
        }
        else {
          sig.m_mesh->traceMng()->info() << "Cut() -- Compute X invalid (too small) -- fst_x.length() : " << fst_x.m_patch.length() << " -- snd_x.length() : " << snd_x.m_patch.length();
        }
      }

      Real y_efficacity = 0;
      auto [fst_y, snd_y] = sig.cut(MD_DirY, cut_point_y);
      {
        sig.m_mesh->traceMng()->info() << "Cut() -- Compute Y -- Cut Point : " << cut_point_y;

        fst_y.compute();
        snd_y.compute();
        if (fst_y.isValid() && snd_y.isValid()) {

          sig.m_mesh->traceMng()->info() << "Cut() -- Y.fst_y"
            << " -- min = " << fst_y.m_patch.m_min
            << " -- max = " << fst_y.m_patch.m_max
            << " -- efficacity : " << fst_y.efficacity()
          ;
          sig.m_mesh->traceMng()->info() << "Cut() -- Y.snd_y"
            << " -- min = " << snd_y.m_patch.m_min
            << " -- max = " << snd_y.m_patch.m_max
            << " -- efficacity : " << snd_y.efficacity()
          ;

          y_efficacity = (fst_y.efficacity() + snd_y.efficacity()) / 2;
          sig.m_mesh->traceMng()->info() << "Cut() -- efficacity Y : " << y_efficacity;
        }
        else {
          sig.m_mesh->traceMng()->info() << "Cut() -- Compute Y invalid (too small) -- fst_y.length() : " << fst_y.m_patch.length() << " -- snd_y.length() : " << snd_y.m_patch.length();
        }
      }

      if (sig.efficacity() > x_efficacity && sig.efficacity() > y_efficacity) {
        return {};
      }

      if (x_efficacity >= y_efficacity && x_efficacity != 0) {
        return {fst_x, snd_x};
      }
      if (y_efficacity == 0) {
        ARCANE_FATAL("Invalid cut");
      }
      return {fst_y, snd_y};
    }

    if (cut_point_x != -1) {
      Real x_efficacity = 0;
      auto [fst_x, snd_x] = sig.cut(MD_DirX, cut_point_x);

      sig.m_mesh->traceMng()->info() << "Cut() -- Compute X -- Cut Point : " << cut_point_x;

      fst_x.compute();
      snd_x.compute();
      if (fst_x.isValid() && snd_x.isValid()) {

        sig.m_mesh->traceMng()->info() << "Cut() -- efficacity X.fst_x : " << fst_x.efficacity();
        sig.m_mesh->traceMng()->info() << "Cut() -- efficacity X.snd_x : " << snd_x.efficacity();
        x_efficacity = (fst_x.efficacity() + snd_x.efficacity()) / 2;
        sig.m_mesh->traceMng()->info() << "Cut() -- efficacity X : " << x_efficacity;
      }
      else {
        sig.m_mesh->traceMng()->info() << "Cut() -- Compute X invalid (too small) -- fst_x.length() : " << fst_x.m_patch.length() << " -- snd_x.length() : " << snd_x.m_patch.length();
      }
      if (sig.efficacity() > x_efficacity) {
        return {};
      }
      return {fst_x, snd_x};
    }

    Real y_efficacity = 0;
    auto [fst_y, snd_y] = sig.cut(MD_DirY, cut_point_y);

    sig.m_mesh->traceMng()->info() << "Cut() -- Compute Y -- Cut Point : " << cut_point_y;

    fst_y.compute();
    snd_y.compute();
    if (fst_y.isValid() && snd_y.isValid()) {

      sig.m_mesh->traceMng()->info() << "Cut() -- efficacity Y.fst_y : " << fst_y.efficacity();
      sig.m_mesh->traceMng()->info() << "Cut() -- efficacity Y.snd_y : " << snd_y.efficacity();
      y_efficacity = (fst_y.efficacity() + snd_y.efficacity()) / 2;
      sig.m_mesh->traceMng()->info() << "Cut() -- efficacity Y : " << y_efficacity;
    }
    else {
      sig.m_mesh->traceMng()->info() << "Cut() -- Compute Y invalid (too small) -- fst_y.length() : " << fst_y.m_patch.length() << " -- snd_y.length() : " << snd_y.m_patch.length();
    }
    if (sig.efficacity() > y_efficacity) {
      return {};
    }
    return {fst_y, snd_y};

  }

  static void cut(UniqueArray<Signature>& sig_array_a)
  {
    UniqueArray<Signature> sig_array_b;
    bool a_a_b_b = false;
    bool need_cut = true;

    while (need_cut) {
      need_cut = false;
      a_a_b_b = !a_a_b_b;

      UniqueArray<Signature>& array_in = a_a_b_b ? sig_array_a : sig_array_b;
      UniqueArray<Signature>& array_out = a_a_b_b ? sig_array_b : sig_array_a;

      for (Integer i = 0; i < array_in.size(); ++i) {
        Signature sig = array_in[i];
        sig.m_mesh->traceMng()->info() << "Cut() -- i : " << i;

        if (!sig.m_stop_cut) {
          if (!sig.m_is_computed) {
            sig.compute();
          }
          if (sig.canBeCut()) {
            auto [fst, snd] = cut(sig);

            if (fst.isValid()) {
              need_cut = true;
              array_out.add(fst);
              array_out.add(snd);
              sig.m_mesh->traceMng()->info() << "First Signature :";
              sig.m_mesh->traceMng()->info() << "\tmin = " << fst.m_patch.m_min << " -- max = " << fst.m_patch.m_max << " -- length = " << fst.m_patch.length();
              sig.m_mesh->traceMng()->info() << "Second Signature :";
              sig.m_mesh->traceMng()->info() << "\tmin = " << snd.m_patch.m_min << " -- max = " << snd.m_patch.m_max << " -- length = " << snd.m_patch.length();
              continue;
            }
            sig.m_mesh->traceMng()->info() << "Invalid Signature";
            sig.m_stop_cut = true;
          }
          else {
            sig.m_stop_cut = true;
          }
        }
        sig.m_mesh->traceMng()->info() << "No Update";
        sig.m_mesh->traceMng()->info() << "\tmin = " << sig.m_patch.m_min << " -- max = " << sig.m_patch.m_max;
        array_out.add(sig);
      }
      array_in.clear();
    }
    if (a_a_b_b) {
      sig_array_a.clear();
      sig_array_a = sig_array_b;
    }
  }
};



/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";

  m_cartesian_mesh = ICartesianMesh::getReference(mesh());

  Ref<ICartesianMeshAMRPatchMng> coarser = CartesianMeshUtils::cartesianMeshAMRPatchMng(m_cartesian_mesh);
  coarser->createSubLevel();
  m_cartesian_mesh->computeDirections();

  m_numbering = CartesianMeshUtils::cartesianMeshNumberingMng(m_cartesian_mesh);

  const auto* m_generation_info = ICartesianMeshGenerationInfo::getReference(mesh(), true);
  m_global_deltat = 0.01;

  m_global_length = m_generation_info->globalLength();
  m_origin = m_generation_info->globalOrigin();

  {
    Int64ConstArrayView nb_cells = m_generation_info->globalNbCells();
    m_nb_cells.x = nb_cells[MD_DirX];
    m_nb_cells.y = nb_cells[MD_DirY];
    m_nb_cells.z = (mesh()->dimension() == 2 ? 1 : nb_cells[MD_DirZ]);
  }

  m_cell_size.x = m_global_length.x / static_cast<Real>(m_nb_cells.x);
  m_cell_size.y = m_global_length.y / static_cast<Real>(m_nb_cells.y);
  m_cell_size.z = m_global_length.z / static_cast<Real>(m_nb_cells.z);

  info() << "Global nb_cells : " << m_nb_cells;
  info() << "Global length : " << m_global_length;
  info() << "Global cell_size : " << m_cell_size;
  info() << "Global origin : " << m_origin;


  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  ENUMERATE_ (Cell, icell, ownCells()) {
    Real3 n0_coord = node_coords[icell->node(0)];
    Real3 n1_coord = node_coords[(mesh()->dimension() == 2 ? icell->node(2) : icell->node(6))];

    Real3 center = n0_coord + ((n1_coord - n0_coord) / 2);

    Real r2 = (std::pow(center.x - 0.5, 2) + std::pow((center.y - 0.75), 2)) / 0.01;
    m_phi[icell] = 1.0 + std::exp(-r2);
  }
  m_phi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  m_cartesian_mesh->computeDirections();

  computeVelocity(m_global_time());

  computePhi();

  markCellsToRefine();

  refine();

  if (globalIteration() > options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computeVelocity(Real time)
{
  constexpr Real PI = 3.1415926535897932384626;

  FaceDirectionMng fdm_x{ m_cartesian_mesh->faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ m_cartesian_mesh->faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ m_cartesian_mesh->cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ m_cartesian_mesh->cellDirection(MD_DirY) };

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  ENUMERATE_ (Cell, icell, ownCells()) {

    Real3 n0_coord = node_coords[icell->node(0)];
    Real3 n1_coord = node_coords[(mesh()->dimension() == 2 ? icell->node(2) : icell->node(6))];

    Real3 center = n0_coord + ((n1_coord - n0_coord) / 2);

    m_psi[icell] = std::pow(std::sin(PI * center.x), 2) * std::pow(std::sin(PI * center.y), 2) * std::cos(PI * time / 2.0) * 1.0 / PI;
  }

  m_psi.synchronize();

  {
    /*
     * |------|------|
     * |  01  |  11  |
     * |------|------|
     * | prev # next |  (iface = #)
     * |------|------|
     * |  00  |  10  |
     * |------|------|
     */
    ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
      DirFace cells_of_face(fdm_x[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      DirCell dir_y_prev_cell(cdm_y[prev]);
      Cell c00 = dir_y_prev_cell.previous();
      Cell c01 = dir_y_prev_cell.next();
      if (c00.null() || c01.null())
        continue;

      DirCell dir_y_next_cell(cdm_y[next]);
      Cell c10 = dir_y_next_cell.previous();
      Cell c11 = dir_y_next_cell.next();
      if (c10.null() || c11.null())
        continue;

      m_velocity[iface] = -((m_psi[c11] + m_psi[c01]) - (m_psi[c10] + m_psi[c00])) * (0.25 / m_cell_size.y);
    }

    /*
     * |------|------|------|
     * |  01  | next |  11  |
     * |------|--##--|------|  (iface = --##--)
     * |  00  | prev |  10  |
     * |------|------|------|
     */
    ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
      DirFace cells_of_face(fdm_y[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      DirCell dir_x_prev_cell(cdm_x[prev]);
      Cell c00 = dir_x_prev_cell.previous();
      Cell c10 = dir_x_prev_cell.next();
      if (c00.null() || c10.null())
        continue;

      DirCell dir_x_next_cell(cdm_x[next]);
      Cell c01 = dir_x_next_cell.previous();
      Cell c11 = dir_x_next_cell.next();
      if (c01.null() || c11.null())
        continue;

      m_velocity[iface] = ((m_psi[c11] + m_psi[c10]) - (m_psi[c01] + m_psi[c00])) * (0.25 / m_cell_size.x);
    }
  }
  m_velocity.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computePhi()
{
  Real dtdx = m_global_deltat() / m_cell_size.x;
  Real dtdy = m_global_deltat() / m_cell_size.y;

  FaceDirectionMng fdm_x{ m_cartesian_mesh->faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ m_cartesian_mesh->faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ m_cartesian_mesh->cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ m_cartesian_mesh->cellDirection(MD_DirY) };

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  VariableFaceReal phixy = VariableBuildInfo(mesh(), "Phixy");

  VariableCellReal phi = VariableBuildInfo(mesh(), "Phi");
  phi.copy(m_phi);

  {
    VariableCellReal slope2 = VariableBuildInfo(mesh(), "Slope2");
    VariableCellReal slope4 = VariableBuildInfo(mesh(), "Slope4");
    {
      ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
        DirCell dir_cell(cdm_x[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        slope2[icell] = dsgn * std::min(dlim, abs(dcen));
      }
      slope2.synchronize();

      ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
        DirCell dir_cell(cdm_x[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        Real dq1 = 4.0 / 3.0 * dcen - (1.0 / 6.0) * (slope2[next] + slope2[prev]);
        slope4[icell] = dsgn * std::min(dlim, abs(dq1));
      }
      slope4.synchronize();

      ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
        DirFace cells_of_face(fdm_x[iface]);
        Cell prev = cells_of_face.previousCell();
        Cell next = cells_of_face.nextCell();

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi[next] - slope4[next] * (0.5 + 0.5 * dtdx * m_velocity[iface]) :
          phi[prev] + slope4[prev] * (0.5 - 0.5 * dtdx * m_velocity[iface]));
      }
    }

    {
      ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
        DirCell dir_cell(cdm_y[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        slope2[icell] = dsgn * std::min(dlim, abs(dcen));
      }
      slope2.synchronize();

      ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
        DirCell dir_cell(cdm_y[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        Real dq1 = 4.0 / 3.0 * dcen - (1.0 / 6.0) * (slope2[next] + slope2[prev]);
        slope4[icell] = dsgn * std::min(dlim, abs(dq1));
      }
      slope4.synchronize();

      ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
        DirFace cells_of_face(fdm_y[iface]);
        Cell prev = cells_of_face.previousCell();
        Cell next = cells_of_face.nextCell();

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi[next] - slope4[next] * (0.5 + 0.5 * dtdy * m_velocity[iface]) :
          phi[prev] + slope4[prev] * (0.5 - 0.5 * dtdy * m_velocity[iface]));
      }
    }
  }
  phixy.synchronize();

  VariableFaceReal tflux = VariableBuildInfo(mesh(), "Tflux");
  {
    ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
      DirFace cells_of_face(fdm_x[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      Face next_f2 = next.face(2);
      Face next_f0 = next.face(0);

      Face prev_f2 = prev.face(2);
      Face prev_f0 = prev.face(0);

      tflux[iface] = ((m_velocity[iface] < 0) ?
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[next_f2] + m_velocity[next_f0]) * (phixy[next_f2] - phixy[next_f0]))) * m_velocity[iface] :
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[prev_f2] + m_velocity[prev_f0]) * (phixy[prev_f2] - phixy[prev_f0]))) * m_velocity[iface]);
    }

    ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
      DirFace cells_of_face(fdm_y[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      Face next_f1 = next.face(1);
      Face next_f3 = next.face(3);

      Face prev_f1 = prev.face(1);
      Face prev_f3 = prev.face(3);

      tflux[iface] = ((m_velocity[iface] < 0) ?
        (phixy[iface] - 0.5 * dtdx * (0.5 * (m_velocity[next_f1] + m_velocity[next_f3]) * (phixy[next_f1] - phixy[next_f3]))) * m_velocity[iface] :
        (phixy[iface] - 0.5 * dtdx * (0.5 * (m_velocity[prev_f1] + m_velocity[prev_f3]) * (phixy[prev_f1] - phixy[prev_f3]))) * m_velocity[iface]);
    }
  }

  {
    ENUMERATE_ (Cell, icell, ownCells()) {
      Face f0 = icell->face(0);
      Face f1 = icell->face(1);
      Face f2 = icell->face(2);
      Face f3 = icell->face(3);

      m_phi[icell] = phi[icell] + ((tflux[f3] - tflux[f1]) * dtdx + (tflux[f0] - tflux[f2]) * dtdy);
    }
  }

  m_phi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
markCellsToRefine()
{

  ENUMERATE_ (Cell, icell, ownCells()) {
    if (icell->level() != 1) continue;
    Integer pos_x = m_numbering->cellUniqueIdToCoordX(*icell);
    Integer pos_y = m_numbering->cellUniqueIdToCoordY(*icell);

    // if (pos_x >= 2 && pos_x < 6 && pos_y >= 2 && pos_y < 5) {
    //   icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    // }
    //
    // if (pos_x >= 7 && pos_x < 11 && pos_y >= 6 && pos_y < 9) {
    //   icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    // }

    if (pos_x >= 3 && pos_x < 11 && pos_y >= 25 && pos_y < 37) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }

    if (pos_x >= 19 && pos_x < 27 && pos_y >= 2 && pos_y < 19) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 19 && pos_x < 27 && pos_y >= 43 && pos_y < 60) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }

    if (pos_x >= 5 && pos_x < 12 && pos_y >= 19 && pos_y < 29) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 7 && pos_x < 13 && pos_y >= 17 && pos_y < 26) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 9 && pos_x < 15 && pos_y >= 15 && pos_y < 23) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 11 && pos_x < 16 && pos_y >= 13 && pos_y < 22) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 15 && pos_x < 18 && pos_y >= 11 && pos_y < 21) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 18 && pos_x < 21 && pos_y >= 11 && pos_y < 20) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }

    if (pos_x >= 5 && pos_x < 12 && pos_y >= 33 && pos_y < 43) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 7 && pos_x < 13 && pos_y >= 36 && pos_y < 45) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 9 && pos_x < 15 && pos_y >= 39 && pos_y < 47) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 11 && pos_x < 16 && pos_y >= 40 && pos_y < 49) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 15 && pos_x < 18 && pos_y >= 41 && pos_y < 51) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
    if (pos_x >= 18 && pos_x < 21 && pos_y >= 42 && pos_y < 51) {
      icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
    }
  }

  Integer nb_cell_x = m_numbering->globalNbCellsX(1);

  StringBuilder str = "";
  ENUMERATE_(Cell, icell, ownCells()) {
    if (icell->level() != 1) continue;
    if (icell->uniqueId().asInt32() % nb_cell_x == 0) {
      str += "\n";
    }
    if (icell->hasFlags(ItemFlags::II_Refine)) {
      str += "[XX]";
    }
    else {
      str += "[..]";
    }
  }
  info() << str;

  nb_cell_x = m_numbering->globalNbCellsX(0);
  str = "";
  ENUMERATE_(Cell, icell, ownCells()) {
    if (icell->level() != 0) continue;
    if (icell->uniqueId().asInt32() % nb_cell_x == 0) {
      str += "\n";
    }
    if (icell->hasFlags(ItemFlags::II_Refine)) {
      str += "[XX]";
    }
    else {
      str += "[..]";
    }
  }
  info() << str;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
refine()
{
  Integer min_level = 0;
  Integer max_level = -1;

  ENUMERATE_(Cell, icell, ownCells()) {
    if (icell->hasFlags(ItemFlags::II_Refine)) {
      Integer level = icell->level();
      if (level > max_level) max_level = level;

      // if (level > 0) {
      //   Cell parent = icell->hParent();
      //   if (!parent.hasFlags(ItemFlags::II_Refine)) {
      //     parent.mutableItemBase().addFlags(ItemFlags::II_Refine); // TODO : Correction ou crash ?
      //   }
      //   Integer nb_children = parent.nbHChildren();
      //   for (Integer i = 0; i < nb_children; ++i) {
      //     Cell child = parent.hChild(i);
      //     if (!child.hasFlags(ItemFlags::II_Refine)) {
      //       child.mutableItemBase().addFlags(ItemFlags::II_Refine); // TODO : Correction ou crash ?
      //     }
      //   }
      // }
    }
  }
  max_level = parallelMng()->reduce(MessagePassing::ReduceMax, max_level);

  info() << "Min level : " << min_level << " -- Max level : " << max_level;

  LevelPatches all_patches(max_level, mesh());

  for (Integer level = max_level; level >= min_level; --level) {
    if (level != max_level) {
      propage(level, &all_patches);
      info() << "All patch level+1 with margin (can be overlap) : ";
      for (auto& elem : all_patches.m_patches[level+1]) {
        info() << "\tPatch -- min = " << elem.minWithMargin(level) << " -- max = " << elem.maxWithMargin(level);
      }
    }

    Patch all_level;
    all_level.m_level = level;
    all_level.m_min = {0,0,0};
    all_level.m_max = { static_cast<Integer>(m_numbering->globalNbCellsX(level)), static_cast<Integer>(m_numbering->globalNbCellsY(level)), 0};

    Signature sig(all_level, m_numbering.get(), &all_patches);
    UniqueArray<Signature> sig_array;
    sig_array.add(sig);

    CutPatch::cut(sig_array);

    UniqueArray<Patch>& level_patches = all_patches.m_patches[level];
    for (const auto& elem : sig_array) {
      level_patches.add(elem.m_patch);
    }

    Real global_efficacity = 0;
    info() << "All patch : ";
    for (auto& elem : sig_array) {
      info() << "\tPatch -- min = " << elem.m_patch.m_min << " -- max = " << elem.m_patch.m_max << " -- Efficacité : " << elem.efficacity();
      global_efficacity += elem.efficacity();
    }
    global_efficacity /= sig_array.size();
    info() << "Global efficacity : " << global_efficacity;

    StringBuilder str = "";
    ENUMERATE_(Cell, icell, ownCells()) {
      if (icell->level() != level) continue;
      Integer pos_x = m_numbering->cellUniqueIdToCoordX(*icell);
      Integer pos_y = m_numbering->cellUniqueIdToCoordY(*icell);
      Integer patch = -1;
      for (Integer i = 0; i < sig_array.size(); ++i) {
        const Signature& elem = sig_array[i];
        if (elem.isIn(pos_x, pos_y)) {
          if (patch != -1) {
            ARCANE_FATAL("ABCDEFG -- old : {0} -- new : {1}", patch, i);
          }
          patch = i;
        }
      }
      if (patch == -1 && icell->hasFlags(ItemFlags::II_Refine)) {
        ARCANE_FATAL("Bad Patch");
      }

      if (pos_x == 0) {
        str += "\n";
      }
      if (patch != -1) {
        str += "[";
        if (patch < 10) {
          str += " ";
        }
        str += patch;
        str += "]";
      }
      else {
        str += "[..]";
      }
    }
    info() << str;

    // Construction des patchs / Ajout des II_Refine
    // Si niveau != 0
      // Sur level-1, ajout des II_Refine autour des patchs de level (on ne doit pas avoir plus de deux niveaux de différences entre deux mailles).
  }

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
propage(Integer level, LevelPatches* all_patches)
{
  IMesh* mesh = all_patches->m_mesh;

  // D'abord, on vire les mailles flaggées qui sont dans des patchs (pour éviter les doublons de signature).
  ENUMERATE_(Cell, icell, mesh->ownCells()) {
    if (icell->level() == level && icell->hasFlags(ItemFlags::II_Refine)) {
      Integer pos_x = m_numbering->cellUniqueIdToCoordX(*icell);
      Integer pos_y = m_numbering->cellUniqueIdToCoordY(*icell);
      for (auto patch : all_patches->m_patches[level]) {
        if (patch.isInWithMargin(level, pos_x, pos_y)) {
          icell->mutableItemBase().removeFlags(ItemFlags::II_Refine);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
endModule()
{
  info() << "Module SayHello END";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SAYHELLO(SayHelloModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/