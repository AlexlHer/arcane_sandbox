// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include <arcane/core/IMesh.h>
#include <arcane/core/IIndexedIncrementalItemConnectivityMng.h>
#include <arcane/core/IIndexedIncrementalItemConnectivity.h>
#include <arcane/core/IIncrementalItemConnectivity.h>
#include <arcane/core/IDoFFamily.h>
#include <arcane/core/IndexedItemConnectivityView.h>
#include <arcane/core/IParallelMng.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
buildModule()
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";

  // Création de la famille DoF "DoFFamily".
  IItemFamily* dof_family_interface = mesh()->findItemFamily(IK_DoF, "DoFFamily", true);
  IDoFFamily* dof_family = dof_family_interface->toDoFFamily();

  // Création des uniqueIds pour les futurs DoF.
  UniqueArray<Int64> uids;
  ENUMERATE_ (Cell, icell, allCells()) {
    for (Integer i = 0; i < m_nb_dof_per_cell; ++i) {
      uids.add(icell->uniqueId().asInt64() * m_nb_dof_per_cell + i);
    }
  }

  // Création des DoF avec les uniqueIds définis au-dessus.
  UniqueArray<Int32> lids(uids.size());
  dof_family->addDoFs(uids, lids);
  dof_family->endUpdate();

  // Création de la connectivité "AllCells" vers "DoFFamily".
  m_cell_to_dof = mesh()->indexedConnectivityMng()->findOrCreateConnectivity(mesh()->cellFamily(), dof_family_interface, "CellToDoF");

  IIncrementalItemConnectivity* cn = m_cell_to_dof->connectivity();
  ConstArrayView<ItemInternal*> dofs = dof_family_interface->itemsInternal();

  IParallelMng* pm = mesh()->parallelMng();
  Int32 my_rank = pm->commRank();

  // On relie chaque maille à ses DoF et on attribue les proprio.
  Integer index = 0;
  ENUMERATE_ (Cell, icell, allCells()) {
    for (Integer i = 0; i < m_nb_dof_per_cell; ++i) {
      cn->addConnectedItem(icell, DoFLocalId(lids[index]));
      dofs[lids[index]]->setOwner(icell->owner(), my_rank);
      index++;
    }
  }

  // On synchronise les proprio.
  dof_family_interface->notifyItemsOwnerChanged();
  dof_family_interface->computeSynchronizeInfos();

  // Création de la variable.
  m_var_test = makeRef(new VariableDoFArrayReal3x3(VariableBuildInfo(mesh(), "TestVar", "DoFFamily")));

  // On redimensionne la seconde dimension.
  m_var_test->resize(m_var_dim2);

  // On remplit le tableau 2D de "123".
  m_var_test->fill({{0, 1, 2}, {3, 4, 5}, {6, 123., 8}});

  // On explore chaque DoF.
  ENUMERATE_ (DoF, idof, dof_family_interface->allItems()) {

    if ((*m_var_test.get())[idof].size() != m_var_dim2) {
      ARCANE_FATAL("Bad size");
    }

    for (Integer i = 0; i < m_var_dim2; ++i) {
      if ((*m_var_test.get())(idof, i)[2][1] != 123.) {
        ARCANE_FATAL("Bad value 4");
      }

      // On met une valeur dans chaque case de la variable 2D.
      (*m_var_test.get())(idof, i)[2][1] = static_cast<Real>(idof->uniqueId().asInt64()+i);
    }
  }
  m_var_test->synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
compute()
{
  if (globalIteration() % 2 == 0) {
    computeOdd();
  }
  else {
    computeEven();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computeOdd()
{
  info() << "Module SayHello COMPUTEODD";

  IItemFamily* dof_family_interface = mesh()->findItemFamily(IK_DoF, "DoFFamily", false);

  info() << "NbDof : " << dof_family_interface->nbItem();
  if (dof_family_interface->nbItem() != m_nb_dof_per_cell * allCells().size()) {
    ARCANE_FATAL("Bad nbdof");
  }
  info() << "dim1Size : " << m_var_test->asArray().dim1Size();
  if (m_var_test->asArray().dim1Size() != dof_family_interface->nbItem()) {
    ARCANE_FATAL("Bad dim1 var -- Expected: {0} -- Found: {1}", dof_family_interface->nbItem(), m_var_test->asArray().dim1Size());
  }
  info() << "dim2Size : " << m_var_test->asArray().dim2Size();
  if (m_var_test->asArray().dim2Size() != m_var_dim2) {
    ARCANE_FATAL("Bad dim2 var -- Expected: {0} -- Found: {1}", m_var_dim2, m_var_test->asArray().dim2Size());
  }
  info() << "totalNbElement : " << m_var_test->asArray().totalNbElement();
  if (m_var_test->asArray().totalNbElement() != m_var_dim2 * dof_family_interface->nbItem()) {
    ARCANE_FATAL("Bad nbelem var -- Expected: {0} -- Found: {1}", m_var_dim2 * dof_family_interface->nbItem(), m_var_test->asArray().totalNbElement());
  }

  ENUMERATE_ (DoF, idof, dof_family_interface->allItems()) {

    if ((*m_var_test.get())[idof].size() != m_var_dim2) {
      ARCANE_FATAL("Bad size");
    }

    for (Integer i = 0; i < m_var_dim2; ++i) {

      if ((*m_var_test.get())(idof, i)[2][1] != static_cast<Real>(idof->owner())) {
        ARCANE_FATAL("Bad value 3");
      }

      (*m_var_test.get())(idof, i)[2][1] = static_cast<Real>(idof->uniqueId().asInt64() + i);
    }
  }

  m_var_test->synchronize();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
computeEven()
{
  info() << "Module SayHello COMPUTEEVEN";
  IParallelMng* pm = mesh()->parallelMng();

  IItemFamily* dof_family_interface = mesh()->findItemFamily(IK_DoF, "DoFFamily", false);

  info() << "NbDof : " << dof_family_interface->nbItem();
  if (dof_family_interface->nbItem() != m_nb_dof_per_cell * allCells().size()) {
    ARCANE_FATAL("Bad nbdof");
  }
  info() << "dim1Size : " << m_var_test->asArray().dim1Size();
  if (m_var_test->asArray().dim1Size() != dof_family_interface->nbItem()) {
    ARCANE_FATAL("Bad dim1 var -- Expected: {0} -- Found: {1}", dof_family_interface->nbItem(), m_var_test->asArray().dim1Size());
  }
  info() << "dim2Size : " << m_var_test->asArray().dim2Size();
  if (m_var_test->asArray().dim2Size() != m_var_dim2) {
    ARCANE_FATAL("Bad dim2 var -- Expected: {0} -- Found: {1}", m_var_dim2, m_var_test->asArray().dim2Size());
  }
  info() << "totalNbElement : " << m_var_test->asArray().totalNbElement();
  if (m_var_test->asArray().totalNbElement() != m_var_dim2 * dof_family_interface->nbItem()) {
    ARCANE_FATAL("Bad nbelem var -- Expected: {0} -- Found: {1}", m_var_dim2 * dof_family_interface->nbItem(), m_var_test->asArray().totalNbElement());
  }

  ENUMERATE_ (DoF, idof, dof_family_interface->allItems()) {

    if ((*m_var_test.get())[idof].size() != m_var_dim2) {
      ARCANE_FATAL("Bad size");
    }
    for (Integer i = 0; i < m_var_dim2; ++i) {

      if ((*m_var_test.get())(idof, i)[2][1] != static_cast<Real>(idof->uniqueId().asInt64() + i)) {
        ARCANE_FATAL("Bad value 1");
      }

      (*m_var_test.get())(idof, i)[2][1] = static_cast<Real>(pm->commRank());
    }
  }

  m_var_test->synchronize();

  ENUMERATE_ (Cell, icell, allCells()) {
    for (ItemLocalId dof : m_cell_to_dof->view().items(icell)) {
      for (Integer i = 0; i < m_var_dim2; ++i) {
        if ((*m_var_test.get())(DoFLocalId(dof), i)[2][1] != static_cast<Real>(icell->owner())) {
          ARCANE_FATAL("Bad value 2");
        }
      }
    }
  }
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
