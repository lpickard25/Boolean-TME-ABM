cmake_minimum_required(VERSION 3.0)
project(Tumor_Immune_CBM)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++17 -O3")

include_directories(Model/inc)
include_directories(${CMAKE_SOURCE_DIR}/inc/stb_master)
file(GLOB SOURCES "Model/src/*.cpp")

add_executable(runModel ${SOURCES}
        Model/src/modelVisualization.cpp
        Model/inc/Vessel.h
        Model/inc/Vessel.h
        Model/src/Vessel_General.cpp
        Model/src/Vessel_General.cpp)

find_package(glfw3 3.4 REQUIRED)
find_package(OpenGL REQUIRED)

target_include_directories(runModel
PRIVATE
${GLEW_INCLUDE_DIRS}
)

target_link_libraries(runModel
        PRIVATE
        OpenGL::GL
        glfw
        ${stb_master}
        ${GLEW_LIBRARY})
