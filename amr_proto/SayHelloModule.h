// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H

#include <arcane/core/ITimeLoopMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/utils/Vector3.h>

#include "SayHello_axl.h"


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
  void computeVelocity(Real time);
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

private:

  Real3 m_global_length;
  Real3 m_origin;
  Int64x3 m_nb_cells;
  ICartesianMesh* m_cartesian_mesh;
  Real3 m_cell_size;
};

#endif