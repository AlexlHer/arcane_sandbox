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
  : ArcaneSayHelloObject(mbi) { }

public:
  void startInit() override;
  void compute() override;
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }
};

#endif
