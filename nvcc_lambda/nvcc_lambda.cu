
#include <iostream>
#include <array>
#include <chrono>
#include <iomanip>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


template <typename IndexType, int rank>
struct MyTuple
{};

template <typename IndexType>
struct MyTuple<IndexType, 1>
{
  constexpr MyTuple(IndexType a)
  : id0(a)
  {}

  template<std::size_t Index>
  __host__ __device__ std::tuple_element_t<Index, MyTuple<IndexType, 1>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType> id0 = {};
};

template <typename IndexType>
struct MyTuple<IndexType, 2>
{
  constexpr MyTuple(IndexType a, IndexType b)
  : id0(a, b)
  {}

  template<std::size_t Index>
  __host__ __device__ std::tuple_element_t<Index, MyTuple<IndexType, 2>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType> id0 = {};
};

template <typename IndexType>
struct MyTuple<IndexType, 3>
{
  constexpr MyTuple(IndexType a, IndexType b, IndexType c)
  : id0(a, b, c)
  {}

  template<std::size_t Index>
  __host__ __device__ std::tuple_element_t<Index, MyTuple<IndexType, 3>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType, IndexType> id0 = {};
};

template <typename IndexType>
struct MyTuple<IndexType, 4>
{
  constexpr MyTuple(IndexType a, IndexType b, IndexType c, IndexType d)
  : id0(a, b, c, d)
  {}

  template<std::size_t Index>
  __host__ __device__ std::tuple_element_t<Index, MyTuple<IndexType, 4>>const& get() const&
  {
    return std::get<Index>(id0);
  }

  std::tuple<IndexType, IndexType, IndexType, IndexType> id0 = {};
};


namespace std
{
template <typename IndexType, int rank>
struct tuple_size<MyTuple<IndexType, rank>>
{
  static constexpr size_t value = rank;
};

template <typename IndexType, int rank>
struct tuple_element<0, MyTuple<IndexType, rank>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<1, MyTuple<IndexType, 2>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<1, MyTuple<IndexType, 3>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<1, MyTuple<IndexType, 4>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<2, MyTuple<IndexType, 3>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<2, MyTuple<IndexType, 4>>
{
  using type = int;
};

template <typename IndexType>
struct tuple_element<3, MyTuple<IndexType, 4>>
{
  using type = int;
};
} // namespace std


template<typename ParamArray>
void compute1param(int nb_sds, int nb_blocs, int nb_cells, int nb_nodes, int nb_values)
{
  int size_array = nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values;

  double *data = (double*)malloc(size_array*sizeof(double));

  for(int i = 0; i < size_array; ++i){
    data[i] = 42.0;
  }

  for (int isd = 0; isd < nb_sds; ++isd) {
    for (int ibloc = 0; ibloc < nb_blocs; ++ibloc) {
      for (int icell = 0; icell < nb_cells; ++icell) {
        for (int inode = 0; inode < nb_nodes; ++inode) {

          int start = 0;
          start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;
          start += ibloc * nb_cells * nb_nodes * nb_values;
          start += icell * nb_nodes * nb_values;
          start += inode * nb_values;

          //auto super_lambda = [=] __host__ __device__ (int ivalue) {
          auto super_lambda = [=] __host__ __device__ (const ParamArray params) {
            const auto [ivalue] = params;
            size_t i = start;
            i += ivalue;
            data[i] = i;
          };

          for (int ivalue = 0; ivalue < nb_values; ++ivalue) {
            //super_lambda(ivalue);
            super_lambda({ivalue});
          }
        }
      }
    }
  }

  double moy = 0;
  for(int i = 0; i < size_array; ++i){
    moy += data[i];
  }
  moy /= size_array;

  std::cout << "Moyenne = " << moy << std::endl;

  free(data);
}

template<typename ParamArray>
void compute2params(int nb_sds, int nb_blocs, int nb_cells, int nb_nodes, int nb_values)
{
  int size_array = nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values;

  double *data = (double*)malloc(size_array*sizeof(double));

  for(int i = 0; i < size_array; ++i){
    data[i] = 42.0;
  }

  for (int isd = 0; isd < nb_sds; ++isd) {
    for (int ibloc = 0; ibloc < nb_blocs; ++ibloc) {
      for (int icell = 0; icell < nb_cells; ++icell) {

        int start = 0;
        start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;
        start += ibloc * nb_cells * nb_nodes * nb_values;
        start += icell * nb_nodes * nb_values;

        //auto super_lambda = [=] __host__ __device__ (int inode, int ivalue) {
        auto super_lambda = [=] __host__ __device__ (const ParamArray params) {
          const auto [inode, ivalue] = params;
          size_t i = start;
          i += inode * nb_values;
          i += ivalue;
          data[i] = i;
        };

        for (int inode = 0; inode < nb_nodes; ++inode) {
          for (int ivalue = 0; ivalue < nb_values; ++ivalue) {
            //super_lambda(inode, ivalue);
            super_lambda({inode, ivalue});
          }
        }
      }
    }
  }

  double moy = 0;
  for(int i = 0; i < size_array; ++i){
    moy += data[i];
  }
  moy /= size_array;

  std::cout << "Moyenne = " << moy << std::endl;

  free(data);
}

template<typename ParamArray>
void compute3params(int nb_sds, int nb_blocs, int nb_cells, int nb_nodes, int nb_values)
{
  int size_array = nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values;

  double *data = (double*)malloc(size_array*sizeof(double));

  for(int i = 0; i < size_array; ++i){
    data[i] = 42.0;
  }

  for (int isd = 0; isd < nb_sds; ++isd) {
    for (int ibloc = 0; ibloc < nb_blocs; ++ibloc) {

      int start = 0;
      start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;
      start += ibloc * nb_cells * nb_nodes * nb_values;

      //auto super_lambda = [=] __host__ __device__ (int icell, int inode, int ivalue) {
      auto super_lambda = [=] __host__ __device__ (const ParamArray params) {
        const auto [icell, inode, ivalue] = params;
        size_t i = start;
        i += icell * nb_nodes * nb_values;
        i += inode * nb_values;
        i += ivalue;
        data[i] = i;
      };

      for (int icell = 0; icell < nb_cells; ++icell) {
        for (int inode = 0; inode < nb_nodes; ++inode) {
          for (int ivalue = 0; ivalue < nb_values; ++ivalue) {
            //super_lambda(icell, inode, ivalue);
            super_lambda({icell, inode, ivalue});
          }
        }
      }
    }
  }

  double moy = 0;
  for(int i = 0; i < size_array; ++i){
    moy += data[i];
  }
  moy /= size_array;

  std::cout << "Moyenne = " << moy << std::endl;

  free(data);
}

template<typename ParamArray>
void compute4params(int nb_sds, int nb_blocs, int nb_cells, int nb_nodes, int nb_values)
{
  int size_array = nb_sds * nb_blocs * nb_cells * nb_nodes * nb_values;

  double *data = (double*)malloc(size_array*sizeof(double));

  for(int i = 0; i < size_array; ++i){
    data[i] = 42.0;
  }

  for (int isd = 0; isd < nb_sds; ++isd) {

    int start = 0;
    start += isd * nb_blocs * nb_cells * nb_nodes * nb_values;

    //auto super_lambda = [=] __host__ __device__ (int ibloc, int icell, int inode, int ivalue) {
    auto super_lambda = [=] __host__ __device__ (const ParamArray params) {
      const auto [ibloc, icell, inode, ivalue] = params;
      size_t i = start;
      i += ibloc * nb_cells * nb_nodes * nb_values;
      i += icell * nb_nodes * nb_values;
      i += inode * nb_values;
      i += ivalue;
      data[i] = i;
    };

    for (int ibloc = 0; ibloc < nb_blocs; ++ibloc) {
      for (int icell = 0; icell < nb_cells; ++icell) {
        for (int inode = 0; inode < nb_nodes; ++inode) {
          for (int ivalue = 0; ivalue < nb_values; ++ivalue) {
            //super_lambda(ibloc, icell, inode, ivalue);
            super_lambda({ibloc, icell, inode, ivalue});
          }
        }
      }
    }
  }

  double moy = 0;
  for(int i = 0; i < size_array; ++i){
    moy += data[i];
  }
  moy /= size_array;

  std::cout << "Moyenne = " << moy << std::endl;

  free(data);
}

// __attribute__((noinline))
int main()
{
  int nb_sds = 100;
  int nb_blocs = 20;
  int nb_cells = 30;
  int nb_nodes = 100;
  int nb_values = 100;

  std::cout << std::setprecision (std::numeric_limits<double>::digits10 + 1);

  {
    const auto start{std::chrono::steady_clock::now()};
    compute1param<std::array<int, 1>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_1param_stdarray  : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute2params<std::array<int, 2>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_2params_stdarray : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute3params<std::array<int, 3>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_3params_stdarray : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute4params<std::array<int, 4>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_4params_stdarray : " << std::chrono::duration<double>{end - start} << std::endl;
  }

  std::cout << std::endl << std::endl;

  {
    const auto start{std::chrono::steady_clock::now()};
    compute1param<MyTuple<int, 1>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_1param_MyTuple  : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute2params<MyTuple<int, 2>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_2params_MyTuple : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute3params<MyTuple<int, 3>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_3params_MyTuple : " << std::chrono::duration<double>{end - start} << std::endl;
  }
  std::cout << std::endl;
  {
    const auto start{std::chrono::steady_clock::now()};
    compute4params<MyTuple<int, 4>>(nb_sds, nb_blocs, nb_cells, nb_nodes, nb_values);
    const auto end{std::chrono::steady_clock::now()};
    std::cout << "compute_4params_MyTuple : " << std::chrono::duration<double>{end - start} << std::endl;
  }

  return 0;
}
