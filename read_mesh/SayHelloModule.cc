// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/core/IMesh.h"
#include "arcane/core/IMeshUtilities.h"
#include "arcane/core/Connectivity.h"

#include <arcane/core/IParallelMng.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
buildModule()
{
}

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";
}

void SayHelloModule::
endModule()
{
  IParallelMng* pm = mesh()->parallelMng();
  mesh()->utilities()->writeToFile(String("test") + pm->commRank() + String(".vtk"), "VtkLegacyMeshWriter");
  info() << "Module SayHello END";
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_SAYHELLO(SayHelloModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
