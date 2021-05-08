cmake_minimum_required(VERSION 3.10)

# set the project name
project(zero_thickness)

find_package(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDE_DIRS})

# add the executable
add_executable(zero_thickness main.cpp)

target_link_libraries(zero_thickness ${ARMADILLO_LIBRARIES})

