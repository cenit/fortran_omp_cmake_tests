cmake_minimum_required(VERSION 3.25)

set(CMAKE_VERBOSE_MAKEFILE "FALSE" CACHE BOOL "Create verbose makefile")

if(DEFINED ENV{VCPKG_DEFAULT_TRIPLET})
  message(STATUS "Setting default vcpkg target triplet to $ENV{VCPKG_DEFAULT_TRIPLET}")
  set(VCPKG_TARGET_TRIPLET $ENV{VCPKG_DEFAULT_TRIPLET})
endif()

if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "" FORCE)
  message(STATUS "VCPKG found: $ENV{VCPKG_ROOT}")
  message(STATUS "Using VCPKG integration")
  message(STATUS "VCPKG_TARGET_TRIPLET: ${VCPKG_TARGET_TRIPLET}")
elseif(DEFINED CMAKE_TOOLCHAIN_FILE)
  message(STATUS "Using toolchain: ${CMAKE_TOOLCHAIN_FILE}")
else()
  message(WARNING "vcpkg not found, toolchain not defined, build might fail")
endif()

project(openmp_fortran)
enable_language(Fortran)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}" CACHE PATH "Install prefix" FORCE)
endif()

set(INSTALL_BIN_DIR      "bin"                   CACHE PATH "Path where exe and dll will be installed")
set(INSTALL_LIB_DIR      "lib"                   CACHE PATH "Path where lib will be installed")
set(INSTALL_INCLUDE_DIR  "include"               CACHE PATH "Path where headers will be installed (we will explicitly include the subfolder ${PROJECT_NAME}")
set(INSTALL_CMAKE_DIR    "share/${PROJECT_NAME}" CACHE PATH "Path where cmake configs will be installed")
set(CMAKE_DEBUG_POSTFIX d)

find_package(OpenMP REQUIRED)

#add_executable(test1 ${CMAKE_CURRENT_LIST_DIR}/src/test1.f90)
#target_link_libraries(test1 PRIVATE OpenMP::OpenMP_Fortran)
#install(TARGETS test1 DESTINATION ${INSTALL_BIN_DIR})

add_executable(test2 ${CMAKE_CURRENT_LIST_DIR}/src/test2.f90)
target_link_libraries(test2 PRIVATE OpenMP::OpenMP_Fortran)
install(TARGETS test2 DESTINATION ${INSTALL_BIN_DIR})
