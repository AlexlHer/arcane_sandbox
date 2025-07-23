// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"

#include <arcane/core/IParallelMng.h>
#include <arcane/core/MachineMemoryWindow.h>
using namespace Arcane;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";
  m_loop_sum = 0;
}

void SayHelloModule::
compute()
{
  info() << "Module SayHello COMPUTE";

  m_loop_sum = m_loop_sum() + m_global_iteration();

  if (m_global_iteration() > options()->getNSteps())
    subDomain()->timeLoopMng()->stopComputeLoop(true);

  IParallelMng* pm = parallelMng();

  constexpr Integer nb_elem = 10;

  MachineMemoryWindow<Integer> window(pm, nb_elem);

  ConstArrayView<Int32> machine_ranks(window.machineRanks());
  Integer machine_nb_proc = machine_ranks.size();
  Integer my_rank = pm->commRank();

  debug() << "My rank : " << my_rank << " -- Machine ranks : " << machine_ranks;

  {
    Span av_my_segment(window.segmentView());
    Integer iter = 0;
    for (Integer& elem : av_my_segment) {
      elem = iter * (my_rank + 1);
      iter++;
    }
  }

  window.barrier();

  for (Int32 rank : machine_ranks) {
    Span av_segment(window.segmentView(rank));

    for (Integer i = 0; i < nb_elem; ++i) {
      //debug() << "Test " << i << " : " << av_segment[i] << " -- " << rank;
      if (av_segment[i] != i * (rank + 1)) {
        ARCANE_FATAL("Bad element in memory window -- Expected : {0} -- Found : {1}", (i * (rank + 1)), av_segment[i]);
      }
    }
  }

  window.barrier();

  if (my_rank == machine_ranks[0]) {
    Span av_window(window.windowView());

    for (Integer j = 0; j < machine_nb_proc; ++j) {
      for (Integer i = 0; i < nb_elem; ++i) {
        av_window[i + (j * nb_elem)] = machine_ranks[j];
        //info() << "Test3 " << (i + (j * nb_elem)) << " : " << av_window[i + (j * nb_elem)] << " -- " << machine_ranks[j];
      }
    }
  }

  window.barrier();

  Span av_window(window.windowConstView());

  for (Integer j = 0; j < machine_nb_proc; ++j) {
    for (Integer i = 0; i < nb_elem; ++i) {
      //info() << "Test4 " << (i + (j * nb_elem)) << " : " << av_window[i + (j * nb_elem)] << " -- " << machine_ranks[j];
      if (av_window[i + (j * nb_elem)] != machine_ranks[j]) {
        ARCANE_FATAL("Bad element in memory window -- Expected : {0} -- Found : {1}", machine_ranks[j], av_window[i + (j * nb_elem)]);
      }
    }
  }

  window.barrier();
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