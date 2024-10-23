// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include "arcane/cartesianmesh/CartesianMeshUtils.h"
#include "arcane/cartesianmesh/CartesianMeshCoarsening2.h"
#include "arcane/cartesianmesh/CartesianMeshPatchListView.h"
#include "arcane/cartesianmesh/CartesianMeshRenumberingInfo.h"
#include "arcane/core/IMesh.h"
#include "arcane/core/IItemFamily.h"
#include <arcane/core/IParticleFamily.h>
#include "arcane/utils/MD5HashAlgorithm.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
init()
{
  info() << "Module SayHello INIT";
  m_cartesian_mesh = ICartesianMesh::getReference(mesh());

  if (subDomain()->isContinue()) {
    m_cartesian_mesh->recreateFromDump();
  }
  else {
    m_cartesian_mesh->computeDirections(); // A ne pas appeler avant recreateFromDump().

    {
      Ref<CartesianMeshCoarsening2> coarser = CartesianMeshUtils::createCartesianMeshCoarsening2(m_cartesian_mesh);
      coarser->createCoarseCells();


      m_cartesian_mesh->computeDirections();
      CartesianMeshRenumberingInfo renumbering_info;
      renumbering_info.setRenumberPatchMethod(1);
      renumbering_info.setSortAfterRenumbering(true);
      renumbering_info.setParentPatch(m_cartesian_mesh->amrPatch(1));
      m_cartesian_mesh->renumberItemsUniqueId(renumbering_info);

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

      for(Integer i = 2; i <= 26; i += 8) {
        for(Integer j = 2; j <= 26; j += 8) {
          m_cartesian_mesh->coarseZone2D({i, j}, {2.0, 2.0});
        }
      }
      for(Integer i = 6; i <= 22; i += 8) {
        for(Integer j = 6; j <= 22; j += 8) {
          m_cartesian_mesh->refinePatch2D({i, j}, {2.0, 2.0});
        }
      }
      m_cartesian_mesh->computeDirections();

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

    }
  }




  m_particle_family = mesh()->findItemFamily("ArcaneParticles", true);


  if (subDomain()->isContinue()) {
    UniqueArray<Int64> uid_particles;
    ENUMERATE_(Particle, ipart, m_particle_family->view()){
      uid_particles.add(ipart->uniqueId());
    }
    std::sort(uid_particles.begin(), uid_particles.end());

    MD5HashAlgorithm hash_algo;
    UniqueArray<Byte> hash_result;
    hash_algo.computeHash64(asBytes(uid_particles.constSpan()), hash_result);
    String hash_str = Convert::toHexaString(hash_result);

    info() << "Hash = " << hash_str;

    if(hash_str != m_hash()) {
      ARCANE_FATAL("Hash different");
    }
    info() << "Même hash que précédemment";

  }
  else {
    UniqueArray<Int64> uid_particles;
    UniqueArray<Int32> lid_cells;
    Int64 uid = 0;

    for(Integer i = 0; i < 10; ++i) {
      ENUMERATE_(Cell, icell, allCells()){
        lid_cells.add(icell.localId());
        uid_particles.add(uid++);
      }
    }

    UniqueArray<Int32> lid_particles(uid_particles.size());

    m_particle_family->toParticleFamily()->addParticles(uid_particles, lid_cells, lid_particles);

    // info() << "uid_particles : " << uid_particles;
    // info() << "lid_cells : " << lid_cells;
    // info() << "lid_particles : " << lid_particles;

    m_particle_family->endUpdate();



    std::sort(uid_particles.begin(), uid_particles.end());

    MD5HashAlgorithm hash_algo;
    UniqueArray<Byte> hash_result;
    hash_algo.computeHash64(asBytes(uid_particles.constSpan()), hash_result);
    String hash_str = Convert::toHexaString(hash_result);

    info() << "Hash = " << hash_str;

    m_hash = hash_str;



    {
      VariableNodeReal3& node_coord = nodesCoordinates();

      ENUMERATE_(Cell, icell, allCells()){
        Real3 coord = {};
        ENUMERATE_(Node, inode, icell->nodes()){
          coord += node_coord[inode];
        }
        coord /= icell->nbNode();

        m_cell_center_coord[icell] = coord;
      }

      ENUMERATE_(Particle, ipart, m_particle_family->view()){
        Real3 coord = m_cell_center_coord[ipart->cell()];
        m_particle_coord[ipart] = coord + Real3{0.1, 0.1, 0.1};
      }
    }
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

  // ENUMERATE_(Cell, icell, allCells()){
  //   info() << "Cell : " << icell->uniqueId() << " -- level : " << icell->level();
  // }

  // ENUMERATE_(Particle, ipart, m_particle_family->view()){
  //   info() << "Particle -- Uid : " << ipart->uniqueId()
  //          << " -- Coord : " << m_particle_coord[ipart];
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
