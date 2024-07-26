// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-

#include "SayHelloModule.h"
#include <arcane/IMesh.h>

#include <arcane/utils/NumArray.h>
#include <arcane/core/MeshMDVariableRef.h>
#include "arcane/utils/Array.h"

#ifdef ARCANE_HAS_CUDA
#include "arcane/accelerator/cuda/CudaAccelerator.h"
#endif

#ifdef ARCANE_HAS_SYCL
#include "arcane/accelerator/sycl/SyclAccelerator.h"
#endif

using namespace Arcane;
using namespace Arcane::Accelerator;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void SayHelloModule::
startInit()
{
  info() << "Module SayHello INIT";

  Integer nb_blocs = options()->getNbBlocs();
  m_asyncQueue = makeRef(new AsyncRunQueuePool(*(acceleratorMng()->defaultRunner()), nb_blocs));
}

__attribute__((noinline))
void SayHelloModule::
compute1()
{
  Integer nb_sds = options()->getNbSds();
  Integer nb_blocs = options()->getNbBlocs();
  Integer nb_cells = options()->getNbCells();
  Integer nb_nodes = options()->getNbNodes();
  Integer nb_values = options()->getNbValues();

  IMemoryAllocator* allocator = platform::getDefaultDataAllocator();
#if !defined(ARCANE_HAS_CUDA) && !defined(ARCANE_HAS_SYCL)
#endif

#ifdef ARCANE_HAS_CUDA
  //IMemoryAllocator* allocator = Cuda::getCudaDeviceMemoryAllocator();
  //IMemoryAllocator* allocator = Cuda::getCudaUnifiedMemoryAllocator();
  //IMemoryAllocator* allocator = Cuda::getCudaHostPinnedMemoryAllocator();
#endif

#ifdef ARCANE_HAS_SYCL
  //IMemoryAllocator* allocator = Sycl::getSyclDeviceMemoryAllocator();
  //IMemoryAllocator* allocator = Sycl::getSyclUnifiedMemoryAllocator();
  //IMemoryAllocator* allocator = Sycl::getSyclHostPinnedMemoryAllocator();
#endif

  {
    Timer::Action _(subDomain(), "111");


    UniqueArray<Real> data(allocator, nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values);
    data.fill(42.0);//fill sur cpu

    double *dataptr = data.data();

    for (Integer isd = 0; isd < nb_sds; ++isd) {
      for (Integer ibloc = 0; ibloc < nb_blocs; ++ibloc) {

        auto command1 = makeCommand((*m_asyncQueue.get())[ibloc]);

        size_t start = 0;
        start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;
        start += ibloc * nb_cells * nb_nodes * nb_values;

        auto c = makeLoopRanges(nb_cells, nb_nodes, nb_values);
        command1 << RUNCOMMAND_LOOP (iter, c)
        {
          auto [icell, inode, ivalue] = iter();
          size_t i = start;
          i += icell * nb_nodes * nb_values;
          i += inode * nb_values;
          i += ivalue;
          dataptr[i] = sqrt(i);
        };
      }
    }
    m_asyncQueue->waitAll();
  }
}

template <typename IndexType, int rank>
struct MDTruc
{};

template <typename IndexType>
struct MDTruc<IndexType, 1>
{
  constexpr MDTruc(IndexType a)
  : id0(a)
  {}

  template<std::size_t Index>
  ARCCORE_HOST_DEVICE std::tuple_element_t<Index, MDTruc<IndexType, 1>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType> id0 = {};
};

template <typename IndexType>
struct MDTruc<IndexType, 2>
{
  constexpr MDTruc(IndexType a, IndexType b)
  : id0(a, b)
  {}

  template<std::size_t Index>
  ARCCORE_HOST_DEVICE std::tuple_element_t<Index, MDTruc<IndexType, 2>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType> id0 = {};
};

template <typename IndexType>
struct MDTruc<IndexType, 3>
{
  constexpr MDTruc(IndexType a, IndexType b, IndexType c)
  : id0(a, b, c)
  {}

  template<std::size_t Index>
  ARCCORE_HOST_DEVICE std::tuple_element_t<Index, MDTruc<IndexType, 3>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType, IndexType> id0 = {};
};

template <typename IndexType>
struct MDTruc<IndexType, 4>
{
  constexpr MDTruc(IndexType a, IndexType b, IndexType c)
  : id0(a, b, c)
  {}

  template<std::size_t Index>
  ARCCORE_HOST_DEVICE std::tuple_element_t<Index, MDTruc<IndexType, 3>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType, IndexType, IndexType> id0 = {};
};


namespace std
{
template <typename IndexType, int rank>
struct tuple_size<MDTruc<IndexType, rank>>
{
  static constexpr size_t value = rank;
};

template <typename IndexType, int rank>
struct tuple_element<0, MDTruc<IndexType, rank>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<1, MDTruc<IndexType, 2>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<1, MDTruc<IndexType, 3>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<1, MDTruc<IndexType, 4>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<2, MDTruc<IndexType, 3>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<2, MDTruc<IndexType, 4>>
{
  using type = Arccore::Integer;
};

template <typename IndexType>
struct tuple_element<3, MDTruc<IndexType, 4>>
{
  using type = Arccore::Integer;
};
} // namespace std





__attribute__((noinline))
void SayHelloModule::
compute2()
{
  info() << "Module SayHello COMPUTE";

  Integer nb_sds = options()->getNbSds();
  Integer nb_blocs = options()->getNbBlocs();
  Integer nb_cells = options()->getNbCells();
  Integer nb_nodes = options()->getNbNodes();
  Integer nb_values = options()->getNbValues();

  IMemoryAllocator* allocator = platform::getDefaultDataAllocator();
#if !defined(ARCANE_HAS_CUDA) && !defined(ARCANE_HAS_SYCL)
#endif

#ifdef ARCANE_HAS_CUDA
  //IMemoryAllocator* allocator = Cuda::getCudaDeviceMemoryAllocator();
  //IMemoryAllocator* allocator = Cuda::getCudaUnifiedMemoryAllocator();
  //IMemoryAllocator* allocator = Cuda::getCudaHostPinnedMemoryAllocator();
#endif

#ifdef ARCANE_HAS_SYCL
  //IMemoryAllocator* allocator = Sycl::getSyclDeviceMemoryAllocator();
  //IMemoryAllocator* allocator = Sycl::getSyclUnifiedMemoryAllocator();
  //IMemoryAllocator* allocator = Sycl::getSyclHostPinnedMemoryAllocator();
#endif

  {
    Timer::Action _(subDomain(), "222");

    UniqueArray<Real> data(allocator, nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values);
    
    data.fill(42.0);

    double *dataptr = data.data();

    for (Integer isd = 0; isd < nb_sds; ++isd) {
      for (Integer ibloc = 0; ibloc < nb_blocs; ++ibloc) {


        Integer start = 0;
        start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;
        start += ibloc * nb_cells * nb_nodes * nb_values;

        //using paramsarray = std::array<Integer, 3>;
        //using paramsarray = std::tuple<Integer, Integer, Integer>;
        using paramsarray = MDTruc<Integer, 3>;

         //auto truc = [=] __host__ __device__  (Integer icell, Integer inode, Integer ivalue) {
        auto truc = [=] __host__ __device__ (const paramsarray params) {
          const auto [icell, inode, ivalue] = params;
          size_t i = start;
          i += icell * nb_nodes * nb_values;
          i += inode * nb_values;
          i += ivalue;
          dataptr[i] = sqrt(i);
        };

        auto c = makeLoopRanges(nb_cells, nb_nodes, nb_values);
        for (Integer icell = c.lowerBound<0>(); icell < c.upperBound<0>(); ++icell) {
          for (Integer inode = c.lowerBound<1>(); inode < c.upperBound<1>(); ++inode) {
            for (Integer ivalue = c.lowerBound<2>(); ivalue < c.upperBound<2>(); ++ivalue) {
              //truc(icell, inode, ivalue);
              truc({icell, inode, ivalue});
            }
          }
        }
      }
    }
  }
}

void SayHelloModule::
compute()
{
  compute1();
  compute2();
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
