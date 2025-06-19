// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* main.cc                                                     (C) 2000-2022 */
/*                                                                           */
/* Main direct_cartesian sample.                                             */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/launcher/ArcaneLauncher.h>

#include <arcane/core/IMesh.h>
#include <arcane/core/ISubDomain.h>
#include <arcane/core/VariableTypes.h>
#include <arcane/core/IMeshModifier.h>
#include <arcane/utils/List.h>

#include <arcane/utils/ITraceMng.h>

using namespace Arcane;

int main(int argc,char* argv[])
{
  ArcaneLauncher::init(CommandLineArguments());
  StandaloneSubDomain launcher{ ArcaneLauncher::createStandaloneSubDomain(String{}) };
  ISubDomain* sd = launcher.subDomain();
  ITraceMng* tm = launcher.traceMng();
  IParallelMng* pm = sd->parallelMng();
  MeshHandle mh = sd->defaultMeshHandle();

  tm->info() << "created : " << mh.hasMesh();

  if (!mh.hasMesh()) return 1;

  IMesh* mesh = mh.mesh();
  mesh->modifier()->setDynamic(true);
  tm->info() << "nbcell : " <<  mesh->nbCell();

  return 0;
}
