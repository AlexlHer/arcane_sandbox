// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H

#include <arcane/core/ITimeLoopMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/utils/Vector3.h>
#include "arcane/core/IPostProcessorWriter.h"

#include "SayHello_axl.h"

#include <arcane/cartesianmesh/ICartesianMeshNumberingMng.h>

struct LevelPatches;
using namespace Arcane;

class SayHelloModule
: public ArcaneSayHelloObject
{
 public:

  explicit SayHelloModule(const ModuleBuildInfo& mbi)
  : ArcaneSayHelloObject(mbi)
  {
  }

 public:

  void startInit() override;
  void compute() override;
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 public:

  void syncUp(Integer level_down, VariableCellReal& var);
  void syncDown(Integer level_down, VariableCellReal& var);
  void computePsi(Real time, CartesianPatch& patch);
  void computeVelocity(CartesianPatch& patch);
  void computePhi(CartesianPatch& patch, VariableCellReal& phi_tmp);
  void testMarkCellsToRefine();
  bool markCellsToRefine(Integer max_level);

 private:

  Real3 m_global_length;
  Real3 m_origin;
  Int32x3 m_nb_cells;
  ICartesianMesh* m_cartesian_mesh;
  Ref<ICartesianMeshNumberingMng> m_numbering;
  Real3 m_cell_size;
  UniqueArray<Real> times;
};

#endif