// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H
 
#include <arcane/ITimeLoopMng.h>

#include <arcane/accelerator/core/IAcceleratorMng.h>
#include <arcane/accelerator/core/RunQueue.h>
#include <arcane/accelerator/RunCommandEnumerate.h>
#include <arcane/accelerator/RunCommandLoop.h>

#include "arcane/accelerator/AsyncRunQueuePool.h"
 
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
  void compute1();
  void compute2();
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 public:
  Ref<Arcane::Accelerator::AsyncRunQueuePool> m_asyncQueue;
};
 
#endif
