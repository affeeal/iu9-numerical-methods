cmake_minimum_required(VERSION 3.22.1)
project(lab2)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(${PROJECT_NAME}
  ${SOURCES}
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
