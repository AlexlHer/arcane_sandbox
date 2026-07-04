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

struct EdgeLite
{
  EdgeLite(const Node aa, const Node bb, const Real3& inter)
  : m_inter(inter)
  {
    const Int32 a = aa.localId();
    const Int32 b = bb.localId();
    if (a < b) {
      m_a = a;
      m_b = b;
    }
    else {
      m_a = b;
      m_b = a;
    }
  }

  EdgeLite() = default;

  bool operator<(const EdgeLite& other) const
  {
    if (m_a != other.m_a) {
      return m_a < other.m_a;
    }
    return m_b < other.m_b;
  }

  bool operator==(const EdgeLite& other) const
  {
    return m_a == other.m_a && m_b == other.m_b;
  }

  Int32 m_a = -1;
  Int32 m_b = -1;
  Real3 m_inter{ -1 };
  Int64 m_uid = -1;
  Int32 d_log = -1;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  // Paramètres du plan : point appartenant au plan et vecteur normal.
  Real3 p0{ 1.6, 0, 0 };
  Real3 normal{ 1, 0.5, 0.1 };

  // Normalise le vecteur normal pour garantir un calcul correct des distances.
  normal = math::normalizeReal3(normal);

  VariableNodeReal3& node_coord = mesh()->nodesCoordinates();
  VariableNodeReal node_dist(VariableBuildInfo(mesh(), "NodeDist"));

  // Calcul de la distance signée de chaque noeud par rapport au plan de coupe.
  ENUMERATE_ (Node, inode, allNodes()) {
    node_dist[inode] = math::dot({ node_coord[inode] - p0 }, normal);
  }

  UniqueArray<Real3> point_coords;
  UniqueArray<EdgeLite> point_coords11;

  UniqueArray<Int64> cells_infos;
  cells_infos.reserve(10000);

  Int32 nb_cell = 0;
  UniqueArray<EdgeLite> tmp11;

  ENUMERATE_ (Cell, icell, allCells()) {
    ItemWithNodes item = *icell;

    // Faut-il continuer ?
    bool cont = true;
    // Si au moins un noeud est sur le plan, peut-être que la face est colinéaire.
    bool has_egal = false;
    {
      // On vérifie que la maille est dans le plan.
      bool has_neg = false, has_pos = false;
      ENUMERATE_ (Node, inode, item.nodes()) {
        Real d = node_dist[inode];
        if (d < 0)
          has_neg = true;
        else if (d > 0)
          has_pos = true;
        else
          has_egal = true;
      }
      if (!(has_neg && has_pos))
        cont = false;
    }
    if (!cont && has_egal) {
      ENUMERATE_ (Face, iface, icell->faces()) {
        has_egal = true;
        ENUMERATE_ (Node, inode, iface->nodes()) {
          Real d = node_dist[inode];
          if (d != 0) {
            has_egal = false;
            break;
          }
        }
        if (has_egal) {
          item = *iface;
          cont = true;
          break;
        }
      }
    }
    if (!cont) {
      continue;
    }

    // Tableaux définissant les arêtes pour chaque type d'élément.
    // hexa_edge: 12 arêtes pour un hexaédre (8 noeuds), chaque arête = [n0, n1] indices locaux.
    static constexpr Integer hexa_edge[12][2] = {
      { 0, 1 },
      { 1, 2 },
      { 2, 3 },
      { 3, 0 },
      { 0, 4 },
      { 1, 5 },
      { 2, 6 },
      { 3, 7 },
      { 4, 5 },
      { 5, 6 },
      { 6, 7 },
      { 7, 4 }
    };
    // tetra_edge: 6 arêtes pour un tétraèdre (4 noeuds).
    static constexpr Integer tetra_edge[6][2] = {
      { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 }
    };
    // quad_edge: 4 arêtes pour un quadrangle (4 noeuds).
    static constexpr Integer quad_edge[4][2] = {
      { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }
    };
    // tri_edge: 3 arêtes pour un triangle (3 noeuds).
    static constexpr Integer tri_edge[3][2] = {
      { 0, 1 }, { 1, 2 }, { 2, 0 }
    };

    Integer nb_edges = 0;
    const Integer(*edge_def)[2] = nullptr;

    Integer nb_node = item.nbNode();

    if (nb_node == 8) {
      nb_edges = 12;
      edge_def = hexa_edge;
    }
    else if (nb_node == 4) {
      nb_edges = 6;
      edge_def = tetra_edge;
    }
    else if (nb_node == 4) {
      nb_edges = 4;
      edge_def = quad_edge;
    }
    else if (nb_node == 3) {
      nb_edges = 3;
      edge_def = tri_edge;
    }
    else {
      ARCANE_FATAL("Type de maille non supporté: nbNode={0}", nb_node);
    }

    // Itère sur toutes les arêtes de la maille.
    for (Integer i = 0; i < nb_edges; ++i) {
      Node node0 = item.node(edge_def[i][0]);
      Node node1 = item.node(edge_def[i][1]);

      info() << "N0 : " << node0.localId() << " -- D = " << std::abs(node_dist[node0]);
      info() << "N1 : " << node0.localId() << " -- D = " << std::abs(node_dist[node1]);

      if (math::isNearlyZero(node_dist[node0])) {
        info() << "AAA";
        for (Integer j = 0; j < 2; ++j) {
          Node n = node0;
          Real3 p = node_coord[n];

          auto aaaa = EdgeLite{ n, n, p };
          aaaa.d_log = 0;

          if (!tmp11.contains(aaaa)) {
            tmp11.add(aaaa);
          }
        }
        continue;
      }
      if (math::isNearlyZero(node_dist[node1])) {
        info() << "BBB";
        for (Integer j = 0; j < 2; ++j) {
          Node n = node1;
          Real3 p = node_coord[n];

          auto aaaa = EdgeLite{ n, n, p };
          aaaa.d_log = 1;

          if (!tmp11.contains(aaaa)) {
            tmp11.add(aaaa);
          }
        }
        continue;
      }

      // Si l'arrête passe à travers le plan.
      if (node_dist[node0] * node_dist[node1] < 0) {
        info() << "CCC";

        // Paramètre d'interpolation t dans [0,1] pour le point d'intersection
        // le long de l'arête de node0 à node1.
        Real t = std::abs(node_dist[node0]) / (std::abs(node_dist[node0]) + std::abs(node_dist[node1]));

        // Calcul du point d'intersection par interpolation linéaire.
        Real3 p;
        p.x = node_coord[node0].x + t * (node_coord[node1].x - node_coord[node0].x);
        p.y = node_coord[node0].y + t * (node_coord[node1].y - node_coord[node0].y);
        p.z = node_coord[node0].z + t * (node_coord[node1].z - node_coord[node0].z);

        auto aaaa = EdgeLite{ node0, node1, p };
        aaaa.d_log = 2;

        if (!tmp11.contains(aaaa)) {
          tmp11.add(aaaa);
        }
        continue;
      }
      info() << "DDD";
    }

    if (tmp11.size() >= 3) {
      // Mappe les indices locaux vers des indices globaux dans point_coords.
      {
        for (auto& edge : tmp11) {
          auto aaa = point_coords11.span().findFirst(edge);
          if (aaa) {
            edge.m_uid = point_coords11[aaa.value()].m_uid;
          }
          else {
            edge.m_uid = point_coords11.size();
            point_coords11.add(edge);
          }
        }
      }

      if (tmp11.size() == 3)
        cells_infos.add(ITI_Triangle3);
      else if (tmp11.size() == 4)
        cells_infos.add(ITI_Quad4);
      else if (tmp11.size() == 5)
        cells_infos.add(ITI_Pentagon5);
      else if (tmp11.size() == 6)
        cells_infos.add(ITI_Hexagon6);
      else
        ARCANE_FATAL("Pas implem : {0}", tmp11.size());

      cells_infos.add(item.uniqueId());

      {
        // Calcul du barycentre de tous les points d'intersection.
        Real3 bary{ 0 };
        for (const auto& node : tmp11) {
          bary += node.m_inter;
        }
        bary /= tmp11.size();

        // On choisit un vecteur de référence arbitraire non parallèle au normal.
        // Si normal.x est grand (proche de l'axe X), utiliser l'axe Y ; sinon utiliser l'axe X.
        Real3 arbitrary = (std::abs(normal.x) > 0.9) ? Real3{ 0.0, 1.0, 0.0 } : Real3{ 1.0, 0.0, 0.0 };
        // u et v forment une base orthonormale du plan de coupe.
        // Ils sont perpendiculaires au normal et l'un à l'autre.
        Real3 u = math::normalizedCrossProduct3(arbitrary, normal);
        Real3 v = math::normalizedCrossProduct3(normal, u);

        // On trie les points d'intersection par angle polaire autour du barycentre.
        // Cela garantit que le polygone résultant est correctement ordonné (sens inverse des aiguilles d'une montre).
        UniqueArray<Int64> indices;
        indices.reserve(tmp11.size());
        for (Int64 i = 0; i < tmp11.size(); ++i) {
          indices.add(i);
        }

        std::sort(indices.begin(), indices.end(),
                  [&](Int64 ia, Int64 ib) {
                    const Real3& pa = tmp11[ia].m_inter;
                    const Real3& pb = tmp11[ib].m_inter;

                    // Vecteurs allant du barycentre vers chaque point.
                    Real3 va{ pa - bary };
                    Real3 vb{ pb - bary };

                    // Projeter sur la base 2D du plan (u, v).
                    Real a_x = math::dot(va, u);
                    Real a_y = math::dot(va, v);

                    Real b_x = math::dot(vb, u);
                    Real b_y = math::dot(vb, v);

                    // Comparer les angles en utilisant atan2.
                    Real angle_a = std::atan2(a_y, a_x);
                    Real angle_b = std::atan2(b_y, b_x);

                    return angle_a < angle_b;
                  });

        for (Int64 idx : indices) {
          ARCANE_FATAL_IF(tmp11[idx].m_uid == -1, "aaa {0}", tmp11[idx].m_uid);
          cells_infos.add(tmp11[idx].m_uid);
        }
      }
      nb_cell++;
    }
    for (auto elem : tmp11) {
      info() << "1UID : " << elem.m_uid << "\t -- Node0 : " << elem.m_a << "\t -- Node1 : " << elem.m_b << "\t -- Log : " << elem.d_log << "\t -- Pos : " << elem.m_inter;
    }
    tmp11.clear();
  }

  for (auto elem : point_coords11) {
    info() << "UID : " << elem.m_uid << "\t -- Node0 : " << elem.m_a << "\t -- Node1 : " << elem.m_b << "\t -- Log : " << elem.d_log << "\t -- Pos : " << elem.m_inter;
  }

  // On crée un nouveau maillage 2D pour stocker la section transversale extraite.
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

  {
    VariableNodeReal3& node_coords(primary_cloned_mesh->nodesCoordinates());
    info() << "Setting " << primary_cloned_mesh->nbNode() << " node coordinates";
    Int64 node_count = 0;
    ENUMERATE_ (Node, inode, primary_cloned_mesh->allNodes()) {
      node_coords[inode] = point_coords11[node_count].m_inter;
      info() << "Node " << node_count << " uid=" << inode->uniqueId() << " coord=" << node_coords[inode];
      node_count++;
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
