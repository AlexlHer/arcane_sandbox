// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H
 
#include <arcane/ITimeLoopMng.h>
#include <arcane/cartesianmesh/ICartesianMesh.h>
#include "arcane/core/IPostProcessorWriter.h"

#include "SayHello_axl.h"

using namespace Arcane;
 
class SayHelloModule
: public ArcaneSayHelloObject
{
 public:
  explicit SayHelloModule(const ModuleBuildInfo& mbi) 
  : ArcaneSayHelloObject(mbi), m_cartesian_mesh(nullptr) { }
 
 public:
  void init() override;
  void compute() override;
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

  ICartesianMesh* m_cartesian_mesh;
  IItemFamily* m_particle_family;
};
 
#endif
