// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/core/IMesh.h"

#include <map>
#include <arcane/core/Directory.h>
#include <arcane/core/IVariableMng.h>
#include <arccore/common/List.h>

#include <arcane/core/MeshBuildInfo.h>
#include <arcane/core/IMeshMng.h>
#include <arcane/core/IMeshFactoryMng.h>
#include <arcane/core/IPrimaryMesh.h>
#include <arcane/core/IMeshModifier.h>
#include <arcane/core/ServiceBuilder.h>
using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
  m_loop_sum = 0;
}

struct Pair
{
  Pair(Int64 aa, Int64 bb)
  {
    if (aa < bb) {
      a = aa;
      b = bb;
    }
    else {
      a = bb;
      b = aa;
    }
  }
  Pair()
  : a(0)
  , b(0)
  {}
  bool operator<(const Pair& other) const
  {
    if (a < other.a && b < other.b) {
      return true;
    }
    return false;
  }
  Int64 a;
  Int64 b;
};

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  constexpr Int64 MAX_NB_NODE = 4;
  Int64 node_uid = 0;

  Real3 p0{ 1.5, 0, 0 };
  Real3 n{ 1, 0, 0 };

  VariableNodeReal3& node_coord = mesh()->nodesCoordinates();
  VariableNodeReal node_dist(VariableBuildInfo(mesh(), "NodeDist"));

  ENUMERATE_ (Node, inode, allNodes()) {
    node_dist[inode] = math::dot({ node_coord[inode] - p0 }, n);
  }

  std::map<Pair, Int64> mapp;
  std::map<Int64, Real3> mapp2;

  UniqueArray<Int64> cells_infos;
  cells_infos.reserve(10000);

  Int32 nb_cell = 0;
  UniqueArray<Real3> tmp;
  UniqueArray<Int64> tmp2;
  ENUMERATE_ (Cell, icell, allCells()) {

    Node node0 = icell->node(icell->nbNode() - 1);
    for (Integer i = 0; i < icell->nbNode(); ++i) {
      Node node1 = icell->node(i);

      if (node_dist[node0] * node_dist[node1] < 0) {
        if ((std::abs(node_dist[node0]) + std::abs(node_dist[node1])) == 0) {
          ARCANE_FATAL("A faire");
        }
        Real t = std::abs(node_dist[node0]) / (std::abs(node_dist[node0]) + std::abs(node_dist[node1]));
        Real3 p;
        p.x = node_coord[node0].x + t * (node_coord[node1].x - node_coord[node0].x);
        p.y = node_coord[node0].y + t * (node_coord[node1].y - node_coord[node0].y);
        p.z = node_coord[node0].z + t * (node_coord[node1].z - node_coord[node0].z);

        tmp.add(p);
        if (!mapp.contains({ node0.uniqueId(), node1.uniqueId() })) {
          mapp[{ node0.uniqueId(), node1.uniqueId() }] = node_uid;
          mapp2[node_uid] = p;
          node_uid++;
        }
        tmp2.add(mapp[{ node0.uniqueId(), node1.uniqueId() }]);
      }

      node0 = node1;
    }
    if (!tmp.empty()) {
      info() << "CellUID : " << icell->uniqueId() << " -- Tmp : " << tmp;
      info() << "CellUID : " << icell->uniqueId() << " -- Tmp2 : " << tmp2;
      cells_infos.add(ITI_Quad4);
      cells_infos.add(icell->uniqueId());
      cells_infos.add(tmp2[0]);
      cells_infos.add(tmp2[1]);
      cells_infos.add(tmp2[3]);
      cells_infos.add(tmp2[2]);

      // for (Integer i = 0; i < tmp.size(); ++i) {
      //   cells_infos.add(tmp2[i]);
      // }
      nb_cell++;
    }
    tmp.clear();
    tmp2.clear();
  }

  IPrimaryMesh* primary_cloned_mesh;
  IMeshMng* mm = subDomain()->meshMng();
  IParallelMng* pm = parallelMng();
  MeshBuildInfo mbi("Mesh1");
  mbi.addParallelMng(makeRef(pm));
  primary_cloned_mesh = mm->meshFactoryMng()->createMesh(mbi);
  primary_cloned_mesh->modifier()->setDynamic(true);
  primary_cloned_mesh->setDimension(2);
  primary_cloned_mesh->endAllocate();
  primary_cloned_mesh->modifier()->addCells(nb_cell, cells_infos);
  primary_cloned_mesh->modifier()->endUpdate();

  VariableNodeReal3& node_coords(primary_cloned_mesh->nodesCoordinates());
  ENUMERATE_ (Node, inode, primary_cloned_mesh->allNodes()) {
    node_coords[inode] = mapp2[inode->uniqueId().asInt64()];
  }

  info() << "New mesh -- NbNode : " << primary_cloned_mesh->nbNode() << " -- NbCells : " << primary_cloned_mesh->nbCell();

  /////////////////////

  ServiceBuilder<IPostProcessorWriter> spp(primary_cloned_mesh->handle());
  Ref<IPostProcessorWriter> pp = spp.createReference("VtkHdfV2PostProcessor");
  Directory output_directory = Directory(subDomain()->exportDirectory(), "amrtestpost1");
  pp->setBaseDirectoryName( output_directory.path());
  IPostProcessorWriter* post_processor = pp.get();
  times.add(m_global_time());
  post_processor->setTimes(times);

  VariableList variables;
  variables.add(primary_cloned_mesh->nodesCoordinates().variable());
  post_processor->setVariables(variables);

  ItemGroupList groups;
  groups.add(primary_cloned_mesh->allNodes());
  post_processor->setGroups(groups);

  IVariableMng* vm = primary_cloned_mesh->variableMng();
  vm->writePostProcessing(post_processor);

  subDomain()->timeLoopMng()->stopComputeLoop(true);
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