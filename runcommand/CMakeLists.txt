cmake_minimum_required(VERSION 3.21)
project(RunCommand LANGUAGES CXX C)

find_package(Arcane REQUIRED)
arcane_accelerator_enable()
add_executable(RunCommand SayHelloModule.cc main.cc SayHello_axl.h)

arcane_accelerator_add_source_files(SayHelloModule.cc)

arcane_generate_axl(SayHello)
arcane_accelerator_add_to_target(RunCommand)
arcane_add_arcane_libraries_to_target(RunCommand)
target_include_directories(RunCommand PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(HelloWorld.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
