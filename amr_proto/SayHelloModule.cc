// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// computeVelocity() and computePhi() and end of startInit() are based on
// AMReX example 06_Advection_Amr available here :
// https://github.com/atmyers/ecp-tutorials/blob/main/06_Advection_Amr
// _cutDim part is based on Berger-Rigoutsos algo.

#include "SayHelloModule.h"

#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/core/IMesh.h>
#include <arcane/core/ArcaneTypes.h>
#include <arcane/utils/StringBuilder.h>
#include <arcane/core/IParallelMng.h>
#include <arcane/utils/ITraceMng.h>
#include <arcane/utils/Vector2.h>
#include <arcane/utils/NumArray.h>
#include <arcane/utils/Array3View.h>

#include "arcane/core/ICartesianMeshGenerationInfo.h"
#include "arcane/core/Directory.h"
#include "arcane/core/IVariableMng.h"
#include "arcane/utils/List.h"
#include "arcane/cartesianmesh/ICartesianMeshAMRPatchMng.h"
#include "arcane/cartesianmesh/CartesianMeshPatchListView.h"
#include "arcane/cartesianmesh/CartesianMeshUtils.h"

#include <arcane/cartesianmesh/CartesianPatch.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/FaceDirectionMng.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";

  m_cartesian_mesh = ICartesianMesh::getReference(mesh());

  // Ref<ICartesianMeshAMRPatchMng> coarser = CartesianMeshUtils::cartesianMeshAMRPatchMng(m_cartesian_mesh);
  // coarser->createSubLevel();
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

  m_inout.fill(0);
  m_indexpatch.fill(0);
  m_velocity.fill(0);
  m_psi.fill(0);

  // Level 0
  {
    VariableCellReal phi_tmp = VariableBuildInfo(mesh(), "Phi_tmp");
    phi_tmp.copy(m_phi);
    for (Integer p = 0; p < m_cartesian_mesh->nbPatch(); ++p) {
      auto patch = m_cartesian_mesh->amrPatch(p);
      if (patch.level() == 0) {
        computePsi(m_global_time(), patch);
        computeVelocity(patch);
        computePhi(patch, phi_tmp);
      }
    }
  }



  // Level l
  for (Integer l = 0; l < 2; ++l) {
    if (markCellsToRefine(l)) {
      info() << "NbPatches before refine : " << m_cartesian_mesh->nbPatch();
      m_cartesian_mesh->refine();
      info() << "NbPatches after refine : " << m_cartesian_mesh->nbPatch();
      syncUp(l, m_phi);

      VariableCellReal phi_tmp = VariableBuildInfo(mesh(), "Phi_tmp");
      phi_tmp.copy(m_phi);

      for (Integer p = 0; p < m_cartesian_mesh->nbPatch(); ++p) {
        auto patch = m_cartesian_mesh->amrPatch(p);
        if (patch.level() == l+1) {
          computePsi(m_global_time(), patch);
        }
      }
      for (Integer p = 0; p < m_cartesian_mesh->nbPatch(); ++p) {
        auto patch = m_cartesian_mesh->amrPatch(p);
        if (patch.level() == l+1) {
          info() << "computeVelocity() with patch : " << p << " -- level : " << l+1 << " -- index : " << patch.index();
          computeVelocity(patch);
          computePhi(patch, phi_tmp);
        }
        else {
          info() << "Found patch with level : " << patch.level();
        }
      }
    }
  }

  // info() << "FaceUID;Cell0;Cell1;Velocity;";
  //
  // ENUMERATE_(Face, iface, allFaces())
  // {
  //   info() << iface->uniqueId() << ";" << iface->cell(0).uniqueId() << ";" << (iface->nbCell() == 2 ? iface->cell(1).uniqueId() : -1l) << ";" << m_velocity[iface] << ";";
  // }
  // info() << "CellUID;Psi;Phi;";
  // ENUMERATE_(Cell, icell, allCells())
  // {
  //   info() << icell->uniqueId() << ";" << m_psi[icell] << ";" << m_phi[icell] << ";";
  // }

  //
  // m_cartesian_mesh->refine();
  //
  info() << "Post-process AMR";
  IPostProcessorWriter* post_processor = options()->postProcessor();
  Directory output_directory = Directory(subDomain()->exportDirectory(),"amrtestpost1");
  output_directory.createDirectory();
  info() << "Creating output dir '" << output_directory.path() << "' for export";
  times.add(m_global_time());
  post_processor->setTimes(times);
  post_processor->setBaseDirectoryName(output_directory.path());

  ItemGroupList groups;
  // groups.add(allCells());
  for (Integer p = 0; p < m_cartesian_mesh->nbPatch(); ++p) {
    auto patch = m_cartesian_mesh->amrPatch(p);
    groups.add(patch.cells());
  }

  post_processor->setGroups(groups);

  IVariableMng* vm = subDomain()->variableMng();

  vm->writePostProcessing(post_processor);

  if (globalIteration() > options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
syncUp(Integer level_down, VariableCellReal& var)
{
  ENUMERATE_(Cell, icell, mesh()->allLevelCells(level_down)){
    if (icell->hasHChildren()) {
      Real value = var[icell];
      for (Integer i = 0; i < icell->nbHChildren(); ++i) {
        Cell child = icell->hChild(i);
        var[child] = value;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
syncDown(Integer level_down, VariableCellReal& var)
{
  ENUMERATE_(Cell, icell, mesh()->allLevelCells(level_down)){
    if (icell->hasHChildren()) {
      Real value = 0;
      for (Integer i = 0; i < icell->nbHChildren(); ++i) {
        Cell child = icell->hChild(i);
        value += var[child];
      }
      value /= icell->nbHChildren();
      var[icell] = value;
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computePsi(Real time, CartesianPatch& patch)
{
  constexpr Real PI = 3.1415926535897932384626;

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  ENUMERATE_ (Cell, icell, patch.cells()) {

    Real3 n0_coord = node_coords[icell->node(0)];
    Real3 n1_coord = node_coords[(mesh()->dimension() == 2 ? icell->node(2) : icell->node(6))];

    Real3 center = n0_coord + ((n1_coord - n0_coord) / 2);

    m_psi[icell] = std::pow(std::sin(PI * center.x), 2) * std::pow(std::sin(PI * center.y), 2) * std::cos(PI * time / 2.0) * 1.0 / PI;

    if (icell->uniqueId() == 3487) {
      info() << "m_psi[3487] : " << m_psi[icell];
    }
  }

  m_psi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computeVelocity(CartesianPatch& patch)
{
  FaceDirectionMng fdm_x{ patch.faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ patch.faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ patch.cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ patch.cellDirection(MD_DirY) };


  // // Utiliser ARCANE_DATA_INIT_POLICY=DEFAULT
  // ENUMERATE_ (Face, iface, patch.cells().faceGroup()){
  //   m_velocity[iface] = 0;
  // }

  ENUMERATE_ (Cell, icell, patch.cells()){
    m_inout[icell] = 0;
    m_uidcells[icell] = icell->uniqueId();
    m_indexpatch[icell] = patch.index();
  }

  ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
    m_inout[icell] += 1;
  }
  ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
    m_inout[icell] += 1;
  }

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

      // if (prev.uniqueId() == 3417 && next.uniqueId() == 3418) {
      //   info() << "3417 3418 : " << iface->uniqueId();
      // }

      // TODO Normalement, en innerFaces, pas besoin.
      if (prev.null() || next.null())
        continue;

      Cell c00;
      Cell c01;
      Cell c10;
      Cell c11;

      auto flemme = [&](Integer i) {
        if (iface->uniqueId() == 7371) {
          info() << i << " c00 : " << c00.null() << " -- c01 : " << c01.null() << " -- cells_of_face.isPreviousCellOwn() : " << cells_of_face.isPreviousCellOwn();
          info() << i << " c10 : " << c10.null() << " -- c11 : " << c11.null() << " -- cells_of_face.isNextCellOwn() : " << cells_of_face.isNextCellOwn();
        }
      };

      // Si next n'est pas dans notre patch, il faut passer par prev+c00 pour
      // avoir c10 et prev+c01 pour avoir c11.
      // TODO : Trouver une solution viable à mettre dans Arcane !
      if (!cells_of_face.isNextCellOwn()) {
        if (!cells_of_face.isPreviousCellOwn()) {
          ARCANE_FATAL("Impossible");
        }
        DirCell dir_y_prev_cell(cdm_y[prev]);
        c00 = dir_y_prev_cell.previous();
        c01 = dir_y_prev_cell.next();
        if (c00.null() || c01.null())
          continue;
        flemme(0);

        DirCell dir_y_c00_cell(cdm_y[c00]);
        c10 = dir_y_c00_cell.next();
        DirCell dir_y_c01_cell(cdm_y[c01]);
        c11 = dir_y_c01_cell.next();
        if (c10.null() || c11.null())
          continue;
      }
      // Si prev n'est pas dans notre patch, il faut passer par next+c10 pour
      // avoir c00 et next+c11 pour avoir c01.
      else if (!cells_of_face.isPreviousCellOwn()) {
        if (!cells_of_face.isNextCellOwn()) {
          ARCANE_FATAL("Impossible");
        }
        DirCell dir_y_next_cell(cdm_y[next]);
        c10 = dir_y_next_cell.previous();
        c11 = dir_y_next_cell.next();
        if (c10.null() || c11.null())
          continue;
        flemme(1);

        DirCell dir_y_c10_cell(cdm_y[c10]);
        c00 = dir_y_c10_cell.previous();
        DirCell dir_y_c11_cell(cdm_y[c11]);
        c01 = dir_y_c11_cell.previous();
        flemme(5);
        if (c00.null() || c01.null())
          continue;
      }
      else {
        DirCell dir_y_prev_cell(cdm_y[prev]);
        c00 = dir_y_prev_cell.previous();
        c01 = dir_y_prev_cell.next();
        if (c00.null() || c01.null())
          continue;
        flemme(2);
        DirCell dir_y_next_cell(cdm_y[next]);
        c10 = dir_y_next_cell.previous();
        c11 = dir_y_next_cell.next();
        flemme(4);
        if (c10.null() || c11.null())
          continue;
      }

      flemme(3);

      m_velocity[iface] = -((m_psi[c11] + m_psi[c01]) - (m_psi[c10] + m_psi[c00])) * (0.25 / m_cell_size.y);

      if (iface->uniqueId() == 7371) {
        info() << "X 7371 " << m_velocity[iface];
      }
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
      // info() << "---";
      // info() << "Face UID : " << iface->uniqueId();
      // info() << "Cell prev UID : " << prev.uniqueId();
      // info() << "Cell next UID : " << next.uniqueId();
      // info() << "---";
      // TODO Normalement, en innerFaces, pas besoin.
      if (prev.null() || next.null())
        continue;

      Cell c00;
      Cell c01;
      Cell c10;
      Cell c11;

      // Si next n'est pas dans notre patch, il faut passer par prev+c00 pour
      // avoir c01 et prev+c10 pour avoir c11.
      // TODO : Trouver une solution viable à mettre dans Arcane !
      if (!cells_of_face.isNextCellOwn()) {
        if (!cells_of_face.isPreviousCellOwn()) {
          ARCANE_FATAL("Impossible");
        }
        DirCell dir_x_prev_cell(cdm_x[prev]);
        c00 = dir_x_prev_cell.previous();
        c10 = dir_x_prev_cell.next();
        if (c00.null() || c10.null())
          continue;

        DirCell dir_x_c00_cell(cdm_x[c00]);
        c01 = dir_x_c00_cell.next();
        DirCell dir_x_c10_cell(cdm_x[c10]);
        c11 = dir_x_c10_cell.next();
        if (c01.null() || c11.null())
          continue;
      }
      // Si prev n'est pas dans notre patch, il faut passer par next+c01 pour
      // avoir c00 et next+c11 pour avoir c10.
      else if (!cells_of_face.isPreviousCellOwn()) {
        if (!cells_of_face.isNextCellOwn()) {
          ARCANE_FATAL("Impossible");
        }
        DirCell dir_x_next_cell(cdm_x[next]);
        c01 = dir_x_next_cell.previous();
        c11 = dir_x_next_cell.next();
        if (c01.null() || c11.null())
          continue;

        DirCell dir_x_c01_cell(cdm_x[c01]);
        c00 = dir_x_c01_cell.previous();
        DirCell dir_x_c11_cell(cdm_x[c11]);
        c10 = dir_x_c11_cell.previous();
        if (c00.null() || c10.null())
          continue;
      }
      else {
        DirCell dir_x_prev_cell(cdm_x[prev]);
        c00 = dir_x_prev_cell.previous();
        c10 = dir_x_prev_cell.next();
        if (c00.null() || c10.null())
          continue;
        DirCell dir_x_next_cell(cdm_x[next]);
        c01 = dir_x_next_cell.previous();
        c11 = dir_x_next_cell.next();
        if (c01.null() || c11.null())
          continue;
      }

      m_velocity[iface] = ((m_psi[c11] + m_psi[c10]) - (m_psi[c01] + m_psi[c00])) * (0.25 / m_cell_size.x);
    }
  }
  m_velocity.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computePhi(CartesianPatch& patch, VariableCellReal& phi_tmp)
{
  Real dtdx = m_global_deltat() / m_cell_size.x;
  Real dtdy = m_global_deltat() / m_cell_size.y;

  FaceDirectionMng fdm_x{ patch.faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ patch.faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ patch.cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ patch.cellDirection(MD_DirY) };

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  VariableFaceReal phixy = VariableBuildInfo(mesh(), "Phixy");

  {
    VariableCellReal slope2 = VariableBuildInfo(mesh(), "Slope2");
    VariableCellReal slope4 = VariableBuildInfo(mesh(), "Slope4");
    {
      ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
        DirCell dir_cell(cdm_x[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        // if (prev.null() || next.null())
        //   continue;

        Real dlft = phi_tmp[icell] - phi_tmp[prev];
        Real drgt = phi_tmp[next] - phi_tmp[icell];
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

        // if (prev.null() || next.null())
        //   continue;

        Real dlft = phi_tmp[icell] - phi_tmp[prev];
        Real drgt = phi_tmp[next] - phi_tmp[icell];
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
        // TODO Normalement, en innerFaces, pas besoin.
        if (prev.null() || next.null())
          continue;

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi_tmp[next] - slope4[next] * (0.5 + 0.5 * dtdx * m_velocity[iface]) :
          phi_tmp[prev] + slope4[prev] * (0.5 - 0.5 * dtdx * m_velocity[iface]));

        if (iface->uniqueId() == 8657) {
          info() << "phixy[8657] = " << phixy[iface];
        }
      }
    }

    {
      ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
        DirCell dir_cell(cdm_y[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        // if (prev.null() || next.null())
        //   continue;

        Real dlft = phi_tmp[icell] - phi_tmp[prev];
        Real drgt = phi_tmp[next] - phi_tmp[icell];
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

        // if (prev.null() || next.null())
        //   continue;

        Real dlft = phi_tmp[icell] - phi_tmp[prev];
        Real drgt = phi_tmp[next] - phi_tmp[icell];
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

        // TODO Normalement, en innerFaces, pas besoin.
        if (prev.null() || next.null())
          continue;

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi_tmp[next] - slope4[next] * (0.5 + 0.5 * dtdy * m_velocity[iface]) :
          phi_tmp[prev] + slope4[prev] * (0.5 - 0.5 * dtdy * m_velocity[iface]));
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
      // TODO Normalement, en innerFaces, pas besoin.
      if (prev.null() || next.null())
        continue;

      Face next_f2 = next.face(2);
      Face next_f0 = next.face(0);

      Face prev_f2 = prev.face(2);
      Face prev_f0 = prev.face(0);

      tflux[iface] = ((m_velocity[iface] < 0) ?
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[next_f2] + m_velocity[next_f0]) * (phixy[next_f2] - phixy[next_f0]))) * m_velocity[iface] :
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[prev_f2] + m_velocity[prev_f0]) * (phixy[prev_f2] - phixy[prev_f0]))) * m_velocity[iface]);

      if (iface->uniqueId() == 8657) {
        info() << "tflux[8657] = " << tflux[iface];
      }
    }

    ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
      DirFace cells_of_face(fdm_y[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();
      // TODO Normalement, en innerFaces, pas besoin.
      if (prev.null() || next.null())
        continue;

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
    ENUMERATE_ (Cell, icell, patch.cells()) {
      Face f0 = icell->face(0);
      Face f1 = icell->face(1);
      Face f2 = icell->face(2);
      Face f3 = icell->face(3);

      m_phi[icell] = phi_tmp[icell] + ((tflux[f3] - tflux[f1]) * dtdx + (tflux[f0] - tflux[f2]) * dtdy);
    }
  }

  m_phi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
testMarkCellsToRefine()
{
  ENUMERATE_ (Cell, icell, mesh()->allLevelCells(0)) {
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

bool SayHelloModule::
markCellsToRefine(Integer max_level)
{
  info() << "markCellsToRefine(" << max_level << ")";
  Real ref = 0.05;
  bool is_edit = false;

  for (Integer i = 0; i <= max_level; ++i) {
    Real ref_level = 1 + ref * pow(10, i);
    info() << "ref_level : " << ref_level;
    ENUMERATE_ (Cell, icell, mesh()->allLevelCells(i)) {
      if (m_phi[icell] > ref_level) {
        icell->mutableItemBase().addFlags(ItemFlags::II_Refine);
        if (max_level == i) {
          is_edit = true;
        }
      }
    }
  }
  if (!is_edit) {
    return false;
  }


  Integer nb_cell_x = m_numbering->globalNbCellsX(max_level);

  StringBuilder str = "";
  ENUMERATE_(Cell, icell, mesh()->allLevelCells(max_level)) {
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

  return true;
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