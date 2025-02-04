// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/cartesianmesh/CartesianMeshUtils.h"
#include "arcane/cartesianmesh/CartesianMeshCoarsening2.h"
#include "arcane/cartesianmesh/CartesianMeshPatchListView.h"
#include "arcane/cartesianmesh/CartesianMeshRenumberingInfo.h"
#include "arcane/core/IMesh.h"
#include "arcane/core/IItemFamily.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
init()
{
  info() << "Module SayHello INIT";

}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";
  subDomain()->timeLoopMng()->stopComputeLoop(true);

  info() << "testOption : " << options()->testOption();
  info() << "boundaryCondition : " << options()->boundaryCondition();
  info() << "pdesRandomNumberGenerator meshName : " << options()->pdesRandomNumberGenerator.meshName();
  info() << "pdesRandomNumberGenerator final meshName : " << options()->pdesRandomNumberGenerator.toICaseOptions()->meshHandle().meshName();
  info() << "pdesRandomNumberGenerator name : " << options()->pdesRandomNumberGenerator.name();
}

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
