cmake_minimum_required(VERSION 3.21)
project(HelloWorld LANGUAGES CXX C)

find_package(Arcane REQUIRED)
arcane_accelerator_enable()
add_executable(HelloWorld SayHelloModule.cc main.cc SayHello_axl.h)

arcane_generate_axl(SayHello)
arcane_add_arcane_libraries_to_target(HelloWorld)
target_include_directories(HelloWorld PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(HelloWorld.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
