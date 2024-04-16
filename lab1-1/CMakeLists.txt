cmake_minimum_required(VERSION 3.15)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_PROJECT_NAME lab0)

project(${CMAKE_PROJECT_NAME})

add_executable(${CMAKE_PROJECT_NAME} src/main.cpp)
