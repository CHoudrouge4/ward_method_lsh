cmake_minimum_required(VERSION 3.10)

project(wards)

set (CMAKE_CXX_STANDARD 11)

include_directories(include)

file(GLOB SOURCES "src/*.cc")

set(CMAKE_CXX_FLAGS " -Wextra -g -O3")

add_executable(wards.exe ${SOURCES})

target_link_libraries(wards.exe)
