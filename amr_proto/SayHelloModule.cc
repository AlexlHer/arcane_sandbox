// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
// Based on AMReX example 06_Advection_Amr available here :
// https://github.com/atmyers/ecp-tutorials/blob/main/06_Advection_Amr

#include "SayHelloModule.h"

#include <arcane/cartesianmesh/ICartesianMesh.h>
#include <arcane/core/IMesh.h>

#include "arcane/core/ICartesianMeshGenerationInfo.h"

#include <arcane/cartesianmesh/CartesianPatch.h>
#include <arcane/cartesianmesh/CellDirectionMng.h>
#include <arcane/cartesianmesh/FaceDirectionMng.h>

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";

  const auto* m_generation_info = ICartesianMeshGenerationInfo::getReference(mesh(), true);
  m_global_deltat = 0.01;

  m_global_length = m_generation_info->globalLength();
  m_origin = m_generation_info->globalOrigin();

  {
    Int64ConstArrayView nb_cells = m_generation_info->globalNbCells();
    m_nb_cells.x = nb_cells[MD_DirX];
    m_nb_cells.y = nb_cells[MD_DirY];
    m_nb_cells.z = (mesh()->dimension() == 2 ? 1 : nb_cells[MD_DirZ]);
  }

  m_cell_size.x = m_global_length.x / static_cast<Real>(m_nb_cells.x);
  m_cell_size.y = m_global_length.y / static_cast<Real>(m_nb_cells.y);
  m_cell_size.z = m_global_length.z / static_cast<Real>(m_nb_cells.z);

  info() << "Global nb_cells : " << m_nb_cells;
  info() << "Global length : " << m_global_length;
  info() << "Global cell_size : " << m_cell_size;
  info() << "Global origin : " << m_origin;

  m_cartesian_mesh = ICartesianMesh::getReference(mesh());

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  ENUMERATE_ (Cell, icell, ownCells()) {
    Real3 n0_coord = node_coords[icell->node(0)];
    Real3 n1_coord = node_coords[(mesh()->dimension() == 2 ? icell->node(2) : icell->node(6))];

    Real3 center = n0_coord + ((n1_coord - n0_coord) / 2);

    Real r2 = (std::pow(center.x - 0.5, 2) + std::pow((center.y - 0.75), 2)) / 0.01;
    m_phi[icell] = 1.0 + std::exp(-r2);
  }
  m_phi.synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  Real dtdx = m_global_deltat() / m_cell_size.x;
  Real dtdy = m_global_deltat() / m_cell_size.y;

  m_cartesian_mesh->computeDirections();

  FaceDirectionMng fdm_x{ m_cartesian_mesh->faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ m_cartesian_mesh->faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ m_cartesian_mesh->cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ m_cartesian_mesh->cellDirection(MD_DirY) };

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  computeVelocity(m_global_time());

  VariableFaceReal phixy = VariableBuildInfo(mesh(), "Phixy");

  VariableCellReal phi = VariableBuildInfo(mesh(), "Phi");
  phi.copy(m_phi);

  {
    VariableCellReal slope2 = VariableBuildInfo(mesh(), "Slope2");
    VariableCellReal slope4 = VariableBuildInfo(mesh(), "Slope4");
    {
      ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
        DirCell dir_cell(cdm_x[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        slope2[icell] = dsgn * std::min(dlim, abs(dcen));
      }
      slope2.synchronize();

      ENUMERATE_ (Cell, icell, cdm_x.innerCells()) {
        DirCell dir_cell(cdm_x[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        Real dq1 = 4.0 / 3.0 * dcen - (1.0 / 6.0) * (slope2[next] + slope2[prev]);
        slope4[icell] = dsgn * std::min(dlim, abs(dq1));
      }
      slope4.synchronize();

      ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
        DirFace cells_of_face(fdm_x[iface]);
        Cell prev = cells_of_face.previousCell();
        Cell next = cells_of_face.nextCell();

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi[next] - slope4[next] * (0.5 + 0.5 * dtdx * m_velocity[iface]) :
          phi[prev] + slope4[prev] * (0.5 - 0.5 * dtdx * m_velocity[iface]));
      }
    }

    {
      ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
        DirCell dir_cell(cdm_y[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        slope2[icell] = dsgn * std::min(dlim, abs(dcen));
      }
      slope2.synchronize();

      ENUMERATE_ (Cell, icell, cdm_y.innerCells()) {
        DirCell dir_cell(cdm_y[icell]);
        Cell next = dir_cell.next();
        Cell prev = dir_cell.previous();

        Real dlft = phi[icell] - phi[prev];
        Real drgt = phi[next] - phi[icell];
        Real dcen = 0.5 * (dlft + drgt);
        Real dsgn = copysign(1.0, dcen);
        Real dslop = 2.0 * ((abs(dlft) < abs(drgt)) ? abs(dlft) : abs(drgt));
        Real dlim = (dlft * drgt >= 0.0) ? dslop : 0.0;
        Real dq1 = 4.0 / 3.0 * dcen - (1.0 / 6.0) * (slope2[next] + slope2[prev]);
        slope4[icell] = dsgn * std::min(dlim, abs(dq1));
      }
      slope4.synchronize();

      ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
        DirFace cells_of_face(fdm_y[iface]);
        Cell prev = cells_of_face.previousCell();
        Cell next = cells_of_face.nextCell();

        phixy[iface] = ((m_velocity[iface] < 0) ?
          phi[next] - slope4[next] * (0.5 + 0.5 * dtdy * m_velocity[iface]) :
          phi[prev] + slope4[prev] * (0.5 - 0.5 * dtdy * m_velocity[iface]));
      }
    }
  }
  phixy.synchronize();

  VariableFaceReal tflux = VariableBuildInfo(mesh(), "Tflux");
  {
    ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
      DirFace cells_of_face(fdm_x[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      Face next_f2 = next.face(2);
      Face next_f0 = next.face(0);

      Face prev_f2 = prev.face(2);
      Face prev_f0 = prev.face(0);

      tflux[iface] = ((m_velocity[iface] < 0) ?
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[next_f2] + m_velocity[next_f0]) * (phixy[next_f2] - phixy[next_f0]))) * m_velocity[iface] :
        (phixy[iface] - 0.5 * dtdy * (0.5 * (m_velocity[prev_f2] + m_velocity[prev_f0]) * (phixy[prev_f2] - phixy[prev_f0]))) * m_velocity[iface]);
    }

    ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
      DirFace cells_of_face(fdm_y[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      Face next_f1 = next.face(1);
      Face next_f3 = next.face(3);

      Face prev_f1 = prev.face(1);
      Face prev_f3 = prev.face(3);

      tflux[iface] = ((m_velocity[iface] < 0) ?
        (phixy[iface] - 0.5 * dtdx * (0.5 * (m_velocity[next_f1] + m_velocity[next_f3]) * (phixy[next_f1] - phixy[next_f3]))) * m_velocity[iface] :
        (phixy[iface] - 0.5 * dtdx * (0.5 * (m_velocity[prev_f1] + m_velocity[prev_f3]) * (phixy[prev_f1] - phixy[prev_f3]))) * m_velocity[iface]);
    }
  }

  {
    ENUMERATE_ (Cell, icell, ownCells()) {
      Face f0 = icell->face(0);
      Face f1 = icell->face(1);
      Face f2 = icell->face(2);
      Face f3 = icell->face(3);

      m_phi[icell] = phi[icell] + ((tflux[f3] - tflux[f1]) * dtdx + (tflux[f0] - tflux[f2]) * dtdy);
    }
  }

  m_phi.synchronize();

  if (globalIteration() > options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computeVelocity(Real time)
{
  constexpr Real PI = 3.1415926535897932384626;

  FaceDirectionMng fdm_x{ m_cartesian_mesh->faceDirection(MD_DirX) };
  FaceDirectionMng fdm_y{ m_cartesian_mesh->faceDirection(MD_DirY) };

  CellDirectionMng cdm_x{ m_cartesian_mesh->cellDirection(MD_DirX) };
  CellDirectionMng cdm_y{ m_cartesian_mesh->cellDirection(MD_DirY) };

  VariableNodeReal3 node_coords = mesh()->nodesCoordinates();

  ENUMERATE_ (Cell, icell, ownCells()) {

    Real3 n0_coord = node_coords[icell->node(0)];
    Real3 n1_coord = node_coords[(mesh()->dimension() == 2 ? icell->node(2) : icell->node(6))];

    Real3 center = n0_coord + ((n1_coord - n0_coord) / 2);

    m_psi[icell] = std::pow(std::sin(PI * center.x), 2) * std::pow(std::sin(PI * center.y), 2) * std::cos(PI * time / 2.0) * 1.0 / PI;
  }

  m_psi.synchronize();

  {
    /*
     * |------|------|
     * |  01  |  11  |
     * |------|------|
     * | prev # next |  (iface = #)
     * |------|------|
     * |  00  |  10  |
     * |------|------|
     */
    ENUMERATE_ (Face, iface, fdm_x.innerFaces()) {
      DirFace cells_of_face(fdm_x[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      DirCell dir_y_prev_cell(cdm_y[prev]);
      Cell c00 = dir_y_prev_cell.previous();
      Cell c01 = dir_y_prev_cell.next();
      if (c00.null() || c01.null())
        continue;

      DirCell dir_y_next_cell(cdm_y[next]);
      Cell c10 = dir_y_next_cell.previous();
      Cell c11 = dir_y_next_cell.next();
      if (c10.null() || c11.null())
        continue;

      m_velocity[iface] = -((m_psi[c11] + m_psi[c01]) - (m_psi[c10] + m_psi[c00])) * (0.25 / m_cell_size.y);
    }

    /*
     * |------|------|------|
     * |  01  | next |  11  |
     * |------|--##--|------|  (iface = --##--)
     * |  00  | prev |  10  |
     * |------|------|------|
     */
    ENUMERATE_ (Face, iface, fdm_y.innerFaces()) {
      DirFace cells_of_face(fdm_y[iface]);
      Cell prev = cells_of_face.previousCell();
      Cell next = cells_of_face.nextCell();

      DirCell dir_x_prev_cell(cdm_x[prev]);
      Cell c00 = dir_x_prev_cell.previous();
      Cell c10 = dir_x_prev_cell.next();
      if (c00.null() || c10.null())
        continue;

      DirCell dir_x_next_cell(cdm_x[next]);
      Cell c01 = dir_x_next_cell.previous();
      Cell c11 = dir_x_next_cell.next();
      if (c01.null() || c11.null())
        continue;

      m_velocity[iface] = ((m_psi[c11] + m_psi[c10]) - (m_psi[c01] + m_psi[c00])) * (0.25 / m_cell_size.x);
    }
  }
  m_velocity.synchronize();

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