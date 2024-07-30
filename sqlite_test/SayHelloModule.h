// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef SAYHELLOMODULE_H
#define SAYHELLOMODULE_H
 
#include <arcane/ITimeLoopMng.h>
#include <sqlite3.h>
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
  void createTable() override;
  void compute() override;
  void printResults() override;
  void endModule() override;
  VersionInfo versionInfo() const override { return VersionInfo(1, 0, 0); }

 private:
  sqlite3* m_db;
};
 
#endif
