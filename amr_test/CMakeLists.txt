cmake_minimum_required(VERSION 3.21)
project(AMRTest LANGUAGES CXX C)

find_package(Arcane REQUIRED)
arcane_accelerator_enable()
add_executable(AMRTest SayHelloModule.cc main.cc SayHello_axl.h)

arcane_generate_axl(SayHello)
arcane_add_arcane_libraries_to_target(AMRTest)
target_include_directories(AMRTest PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(HelloWorld.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
