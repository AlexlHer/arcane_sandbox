cmake_minimum_required(VERSION 3.21)
project(SharedMemory LANGUAGES CXX C)

find_package(Arcane REQUIRED)
arcane_accelerator_enable()
add_executable(SharedMemory SayHelloModule.cc main.cc SayHello_axl.h)

arcane_generate_axl(SayHello)
arcane_add_arcane_libraries_to_target(SharedMemory)
target_include_directories(SharedMemory PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(HelloWorld.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
