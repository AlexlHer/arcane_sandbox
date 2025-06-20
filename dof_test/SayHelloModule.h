// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H

#include <arcane/core/ITimeLoopMng.h>

#include "SayHello_axl.h"

using namespace Arcane;

class SayHelloModule
: public ArcaneSayHelloObject
{
 public:

  explicit SayHelloModule(const ModuleBuildInfo& mbi)
  : ArcaneSayHelloObject(mbi)
  , m_nb_dof_per_cell(5)
  , m_var_dim2(3)
  {}

 public:

  void buildModule() override;
  void startInit() override;
  void compute() override;
  void computeOdd();
  void computeEven();
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:

  Ref<VariableDoFArrayReal3x3> m_var_test;
  Ref<IIndexedIncrementalItemConnectivity> m_cell_to_dof;
  Integer m_nb_dof_per_cell;
  Integer m_var_dim2;
};

#endif
