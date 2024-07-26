$PZC_NVCC_BIN --extended-lambda --expt-relaxed-constexpr -gencode arch=compute_86,code=sm_86 --std=c++20 -g -O3 nvcc_lambda.cu -o with_nvcc.out

$PZC_SYCL_BIN --acpp-targets=cuda:sm_86 --std=c++20 -g -O3 nvcc_lambda.cu -o with_acpp.out

$PZC_CXX_CLANG_BIN -x cuda --cuda-gpu-arch=sm_86 --std=c++20 -g -O3 nvcc_lambda.cu -o with_clang.out
