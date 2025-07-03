// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H

#include <arcane/core/ITimeLoopMng.h>
#include "SayHello_axl.h"

#include <arcane/cartesianmesh/ICartesianMesh.h>

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

private:

  Real3 m_center;
  Real m_radius;
  ICartesianMesh* m_cartesian_mesh;
};

#endif