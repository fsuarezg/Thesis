cmake_minimum_required(VERSION 3.6)
project(Thesis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp pointcloud.cpp pointcloud.h LUTstack.cpp LUTstack.h angular_spectrum_kernel.cpp angular_spectrum_kernel.h hologram.cpp hologram.h)

FIND_PACKAGE(ArrayFire)
INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})

find_package(Boost 1.57.0 COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(BoostTest main.cpp)
target_link_libraries(BoostTest ${Boost_LIBRARIES})

add_executable(Thesis ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(Thesis ${ArrayFire_LIBRARIES} ${EXTRA_LIBS})