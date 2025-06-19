cmake_minimum_required(VERSION 3.21)
project(MeshBackend LANGUAGES CXX C)

find_package(Arcane REQUIRED)
add_executable(MeshBackend main.cc)
arcane_add_arcane_libraries_to_target(MeshBackend)
target_include_directories(MeshBackend PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
