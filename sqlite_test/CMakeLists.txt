cmake_minimum_required(VERSION 3.21)
project(sqlite_test LANGUAGES CXX C)

find_package(Arcane REQUIRED)
find_package(SQLite3 REQUIRED)
add_executable(sqlite_test SayHelloModule.cc main.cc SayHello_axl.h)


target_link_libraries(sqlite_test PUBLIC SQLite::SQLite3)

arcane_generate_axl(SayHello)
arcane_add_arcane_libraries_to_target(sqlite_test)
target_include_directories(sqlite_test PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(HelloWorld.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
