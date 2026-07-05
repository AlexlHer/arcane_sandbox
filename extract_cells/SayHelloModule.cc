// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/core/IMesh.h"

#include <arcane/core/Directory.h>
#include <arcane/core/IVariableMng.h>
#include <arccore/common/List.h>

#include <arcane/core/MeshBuildInfo.h>
#include <arcane/core/IMeshMng.h>
#include <arcane/core/IMeshFactoryMng.h>
#include <arcane/core/IPrimaryMesh.h>
#include <arcane/core/IMeshModifier.h>
#include <arcane/core/ServiceBuilder.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
  m_loop_sum = 0;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  // Paramètres du plan : point appartenant au plan et vecteur normal.
  constexpr Int32 size_array = 2;
  Real3 p0[size_array]{ { -0.2, 0, 0 }, { 0.2, 0, 0 } };
  Real3 normal[size_array]{ { 1, 0, 0 }, { -1, 0, 0 } };

  // Normalise le vecteur normal pour garantir un calcul correct des distances.
  for (Real3& n : normal) {
    n = math::normalizeReal3(n);
  }

  std::unordered_map<Int64, Real3> coord_map;
  UniqueArray<Int64> cells_infos;
  cells_infos.reserve(10000);
  Int32 nb_cell = 0;

  {
    VariableNodeReal3& node_coord = mesh()->nodesCoordinates();

    ENUMERATE_ (Cell, icell, allCells()) {
      Real3 b{ 0 };
      for (Node node : icell->nodes()) {
        b += node_coord[node];
      }
      b /= icell->nbNode();

      bool in_cut = true;
      for (Integer i = 0; i < size_array; ++i) {
        const Real dist = math::dot({ b - p0[i] }, normal[i]);
        if (dist < 0) {
          in_cut = false;
          break;
        }
      }

      if (!in_cut)
        continue;

      Int16 cell_type = icell->itemTypeId();
      cells_infos.add(cell_type);

      Int64 cell_uid = icell->uniqueId().asInt64();
      cells_infos.add(cell_uid);

      for (Node node : icell->nodes()) {
        Int64 node_uid = node.uniqueId().asInt64();
        cells_infos.add(node_uid);
        coord_map[node.uniqueId()] = node_coord[node];
      }
      ++nb_cell;
    }
  }

  // On crée un nouveau maillage 3D pour stocker la section transversale extraite.
  IPrimaryMesh* primary_cloned_mesh;
  IMeshMng* mm = subDomain()->meshMng();
  IParallelMng* pm = parallelMng();
  MeshBuildInfo mbi("Mesh1");
  mbi.addParallelMng(makeRef(pm));
  primary_cloned_mesh = mm->meshFactoryMng()->createMesh(mbi);
  primary_cloned_mesh->modifier()->setDynamic(true);
  primary_cloned_mesh->setDimension(3);
  primary_cloned_mesh->endAllocate();
  primary_cloned_mesh->modifier()->addCells(nb_cell, cells_infos);
  primary_cloned_mesh->modifier()->endUpdate();

  {
    VariableNodeReal3& node_coords(primary_cloned_mesh->nodesCoordinates());
    ENUMERATE_ (Node, inode, primary_cloned_mesh->allNodes()) {
      node_coords[inode] = coord_map[inode->uniqueId()];
    }
  }

  info() << "New mesh -- NbNode : " << primary_cloned_mesh->nbNode() << " -- NbCells : " << primary_cloned_mesh->nbCell();

  /////////////////////

  times.add(m_global_time());
  {
    ServiceBuilder<IPostProcessorWriter> spp(primary_cloned_mesh->handle());
    Ref<IPostProcessorWriter> pp = spp.createReference("VtkHdfV2PostProcessor");
    Directory output_directory = Directory(subDomain()->exportDirectory(), "amrtestpost1");
    pp->setBaseDirectoryName(output_directory.path());
    IPostProcessorWriter* post_processor = pp.get();
    post_processor->setTimes(times);

    VariableList variables;
    variables.add(primary_cloned_mesh->nodesCoordinates().variable());
    post_processor->setVariables(variables);

    ItemGroupList groups;
    groups.add(primary_cloned_mesh->allNodes());
    post_processor->setGroups(groups);

    IVariableMng* vm = primary_cloned_mesh->variableMng();
    vm->writePostProcessing(post_processor);
  }

  {
    ServiceBuilder<IPostProcessorWriter> spp(mesh()->handle());
    Ref<IPostProcessorWriter> pp = spp.createReference("VtkHdfV2PostProcessor");
    Directory output_directory = Directory(subDomain()->exportDirectory(), "amrtestpost1");
    pp->setBaseDirectoryName(output_directory.path());
    IPostProcessorWriter* post_processor = pp.get();
    post_processor->setTimes(times);

    VariableList variables;
    variables.add(mesh()->nodesCoordinates().variable());
    post_processor->setVariables(variables);

    ItemGroupList groups;
    groups.add(mesh()->allNodes());
    post_processor->setGroups(groups);

    IVariableMng* vm = mesh()->variableMng();
    vm->writePostProcessing(post_processor);
  }

  subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

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
