// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"

#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/core/IMesh.h>

#include "arcane/core/ICartesianMeshGenerationInfo.h"

#include <arcane/cartesianmesh/CartesianPatch.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
  m_loop_sum = 0;

  const auto* m_generation_info = ICartesianMeshGenerationInfo::getReference(mesh(), true);
  m_global_deltat = 1;

  m_center = m_generation_info->globalLength() / 2;
  m_radius = 2;

  info() << "Global length : " << m_generation_info->globalLength();
  info() << "Global center : " << m_center;
  info() << "Global radius : " << m_radius;

  m_cartesian_mesh = ICartesianMesh::getReference(mesh());
}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  Real3 center_large_circle = m_center + Real3{ std::cos(globalIteration()), std::sin(globalIteration()), 0 } * m_radius;
  Real radius_large_circle = 5;

  info() << "Large center : " << center_large_circle;
  info() << "Large radius : " << radius_large_circle;

  VariableNodeReal3& node_coords = mesh()->nodesCoordinates();

  info() << "Nb cells : " << m_cartesian_mesh->amrPatch(0).cells().size();

  ENUMERATE_ (Cell, icell, allCells()) {
    bool sup = false;
    bool inf = false;
    for (Node node : icell->nodes()) {
      Real node_dist = math::normL2(node_coords[node] - center_large_circle);
      if (node_dist < radius_large_circle)
        inf = true;
      else
        sup = true;
    }
    Int16 old_amr = m_amr[icell];
    if (sup && inf) {
      m_amr[icell] = 1;
    }
    else {
      m_amr[icell] = 0;
    }

    Real2 min_pos = Real2{node_coords[icell->node(0)]};

    // if (old_amr != m_amr[icell]) {
    //   AMRZonePosition aaaa(min_pos, {1, 1});
    //   if (old_amr == 0) {
    //     m_cartesian_mesh->refinePatch(aaaa);
    //   }
    //   else {
    //     m_cartesian_mesh->coarseZone(aaaa);
    //   }
    // }
  }

  m_loop_sum = m_loop_sum() + m_global_iteration();
  if (globalIteration() > options()->getNSteps())
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