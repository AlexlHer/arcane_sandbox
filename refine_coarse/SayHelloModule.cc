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
  m_cartesian_mesh = ICartesianMesh::getReference(mesh());
  m_cartesian_mesh->computeDirections();

  if (subDomain()->isContinue()) {
    m_cartesian_mesh->recreateFromDump();
  }
  else {
    Ref<CartesianMeshCoarsening2> coarser = CartesianMeshUtils::createCartesianMeshCoarsening2(m_cartesian_mesh);
    coarser->createCoarseCells();

    {
      CartesianMeshPatchListView patches = m_cartesian_mesh->patches();
      Int32 nb_patch = patches.size();
      Int32 index = 0;
      info() << "NB_PATCH=" << nb_patch;
      for( Integer i=0; i<nb_patch; ++i ){
        ICartesianMeshPatch* p = m_cartesian_mesh->patch(i);
        info() << "Patch i=" << index << " nb_cell=" << p->cells().size();
        ++index;
      }
    }

    // m_cartesian_mesh->coarseZone2D({2.0, 2.0}, {4.0, 4.0});
    // m_cartesian_mesh->computeDirections();

    {
      CartesianMeshPatchListView patches = m_cartesian_mesh->patches();
      Int32 nb_patch = patches.size();
      Int32 index = 0;
      info() << "NB_PATCH=" << nb_patch;
      for( Integer i=0; i<nb_patch; ++i ){
        ICartesianMeshPatch* p = m_cartesian_mesh->patch(i);
        info() << "Patch i=" << index << " nb_cell=" << p->cells().size();
        ++index;
      }
    }

    m_cartesian_mesh->computeDirections();
    CartesianMeshRenumberingInfo renumbering_info;
    renumbering_info.setRenumberPatchMethod(1);
    renumbering_info.setSortAfterRenumbering(true);
    renumbering_info.setParentPatch(m_cartesian_mesh->amrPatch(1));
    m_cartesian_mesh->renumberItemsUniqueId(renumbering_info);
  }

}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  {
    CartesianMeshPatchListView patches = m_cartesian_mesh->patches();
    Int32 nb_patch = patches.size();
    Int32 index = 0;
    info() << "NB_PATCH=" << nb_patch;
    for( Integer i=0; i<nb_patch; ++i ){
      ICartesianMeshPatch* p = m_cartesian_mesh->patch(i);
      info() << "Patch i=" << index << " nb_cell=" << p->cells().size();
      ++index;
    }
  }

  ENUMERATE_(Cell, icell, allCells()){
    info() << "Cell : " << icell->uniqueId() << " -- level : " << icell->level();
  }
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
