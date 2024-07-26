// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/core/IMesh.h"
#include "arcane/core/UnstructuredMeshConnectivity.h"
#include "arcane/core/Connectivity.h"
#include "arcane/core/IMeshModifier.h"

#include "arcane/cea/ICartesianMesh.h"
#include "arcane/cea/ICartesianMeshPatch.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"
#include "arcane/cea/NodeDirectionMng.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
buildModule()
{
  Connectivity c(mesh()->connectivity());
  c.enableConnectivity(Connectivity::CT_HasEdge + Connectivity::CT_FaceToEdge + Connectivity::CT_CellToEdge);
}

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
}

void SayHelloModule::
compute()
{
  UnstructuredMeshConnectivityView m_connectivity_view;
  m_connectivity_view.setMesh(mesh());
  Cell cell;
  ENUMERATE_(Cell, icell, allCells()){
    cell = *icell;
    break;
  }
  CellLocalId cid(cell);

  Node node;
  ENUMERATE_(Node, inode, allNodes()){
    node = *inode;
    break;
  }
  NodeLocalId nid(node);

  Face face;
  ENUMERATE_(Face, iface, allFaces()){
    face = *iface;
    break;
  }
  FaceLocalId fid(face);

  {
    auto cfe = m_connectivity_view.cellEdge();
    const Integer nbedgereal = cfe.nbEdge(cid);
    info() << "Truc : " << nbedgereal;
  }
  {
    auto cfe = m_connectivity_view.cellNode();
    const Integer nbedgereal = cfe.nbNode(cid);
    info() << "Truc : " << nbedgereal;
  }

  {
    auto cfe = m_connectivity_view.faceEdge();
    const Integer nbedgereal = cfe.nbEdge(fid);
    info() << "Truc : " << nbedgereal;
  }

  {
    auto cfe = m_connectivity_view.nodeEdge();
    const Integer nbedgereal = cfe.nbEdge(nid);
    info() << "Truc : " << nbedgereal;
  }



  // ICartesianMesh* m_cartesian_mesh = ICartesianMesh::getReference(this->mesh());
  // m_cartesian_mesh->computeDirections();

  // {
  //   CellDirectionMng cdm_x {m_cartesian_mesh->cellDirection(MD_DirX)};
  //   CellDirectionMng cdm_y {m_cartesian_mesh->cellDirection(MD_DirY)};

  //   info() << "cdm_x.allCells().size() : " << cdm_x.allCells().size();
  //   info() << "cdm_y.allCells().size() : " << cdm_y.allCells().size();

  //   ENUMERATE_(Cell, icell, cdm_x.allCells()){
  //     DirCell dir_cell(cdm_x[icell]);
  //     Cell next = dir_cell.next();
  //     Cell prev = dir_cell.previous();
  //     info() << "cdm_x -- prev : " << prev.uniqueId() << " -- actual : " << icell->uniqueId() << " -- next : " << next.uniqueId();
  //   }

  //   ENUMERATE_(Cell, icell, cdm_y.allCells()){
  //     DirCell dir_cell(cdm_y[icell]);
  //     Cell next = dir_cell.next();
  //     Cell prev = dir_cell.previous();
  //     info() << "cdm_y -- prev : " << prev.uniqueId() << " -- actual : " << icell->uniqueId() << " -- next : " << next.uniqueId();
  //   }
  // }

  // {
  //   NodeDirectionMng cdm_x {m_cartesian_mesh->nodeDirection(MD_DirX)};
  //   NodeDirectionMng cdm_y {m_cartesian_mesh->nodeDirection(MD_DirY)};

  //   info() << "cdm_x.allNodes().size() : " << cdm_x.allNodes().size();
  //   info() << "cdm_y.allNodes().size() : " << cdm_y.allNodes().size();

  //   ENUMERATE_(Node, inode, cdm_x.allNodes()){
  //     DirNode dir_cell(cdm_x[inode]);
  //     Node next = dir_cell.next();
  //     Node prev = dir_cell.previous();
  //     info() << "cdm_x -- prev : " << prev.uniqueId() << " -- actual : " << inode->uniqueId() << " -- next : " << next.uniqueId();
  //   }

  //   ENUMERATE_(Node, inode, cdm_y.allNodes()){
  //     DirNode dir_cell(cdm_y[inode]);
  //     Node next = dir_cell.next();
  //     Node prev = dir_cell.previous();
  //     info() << "cdm_y -- prev : " << prev.uniqueId() << " -- actual : " << inode->uniqueId() << " -- next : " << next.uniqueId();
  //   }
  // }

  // {
  //   FaceDirectionMng cdm_x {m_cartesian_mesh->faceDirection(MD_DirX)};
  //   FaceDirectionMng cdm_y {m_cartesian_mesh->faceDirection(MD_DirY)};

  //   info() << "cdm_x.allFaces().size() : " << cdm_x.allFaces().size();
  //   info() << "cdm_y.allFaces().size() : " << cdm_y.allFaces().size();

  //   ENUMERATE_(Face, iface, cdm_x.allFaces()){
  //     DirFace dir_cell(cdm_x[iface]);
  //     Cell next = dir_cell.nextCell();
  //     Cell prev = dir_cell.previousCell();
  //     info() << "cdm_x -- prev : " << prev.uniqueId() << " -- actual : " << iface->uniqueId() << " -- next : " << next.uniqueId();
  //   }

  //   ENUMERATE_(Face, iface, cdm_y.allFaces()){
  //     DirFace dir_cell(cdm_y[iface]);
  //     Cell next = dir_cell.nextCell();
  //     Cell prev = dir_cell.previousCell();
  //     info() << "cdm_y -- prev : " << prev.uniqueId() << " -- actual : " << iface->uniqueId() << " -- next : " << next.uniqueId();
  //   }
  // }
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
