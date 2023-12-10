
cmake_minimum_required(VERSION 2.8)
IF(APPLE)
    SET(CMAKE_C_COMPILER "/usr/local/Cellar/gcc/9.1.0/bin/gcc-9")
    SET(CMAKE_CXX_COMPILER "/usr/local/Cellar/gcc/9.1.0/bin/g++-9")
ENDIF()
project(libgrape-lite C CXX)



# ------------------------------------------------------------------------------
# cmake options
# ------------------------------------------------------------------------------
option(USE_JEMALLOC "Whether to use jemalloc." OFF)
option(USE_HUGEPAGES "Whether to use hugepages." OFF)
option(BUILD_SHARED_LIBS "Whether to build libgrape-lite as shared library" ON)
option(PROFILING "Whether to enable profiling" OFF)
option(WITH_ASAN "Build with Address Sanitizer" OFF)

if (USE_HUGEPAGES AND LINUX)
    add_definitions(-DUSE_HUGEPAGES)
endif ()

if (PROFILING)
    message("-- Enable profiling")
    add_definitions(-DPROFILING)
endif ()

include_directories(include/libgrape-lite/thirdparty)
include_directories(include/GUNDAM/include/)
include_directories(include/GUNDAM/)
include_directories(include/libgrape-lite/)
include_directories(include/)

# ------------------------------------------------------------------------------
# setting default cmake type to Release
# ------------------------------------------------------------------------------

SET(DEFAULT_BUILD_TYPE "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
            STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

message("[libgrape-lite] will build in type: ${CMAKE_BUILD_TYPE}")

# ------------------------------------------------------------------------------
# cmake configs
# ------------------------------------------------------------------------------

include(CheckLibraryExists)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (APPLE)
    SET(CMAKE_MACOSX_RPATH ON)
else ()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Werror -Wl,-rpath,$ORIGIN")
endif ()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -fprofile-arcs -ftest-coverage")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -m64 -flto=auto -march=native")

if (WITH_ASAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -O1")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# ------------------------------------------------------------------------------
# find_libraries
# ------------------------------------------------------------------------------
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

# find Threads------------------------------------------------------------------
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)

# find glog---------------------------------------------------------------------
include("cmake/FindGlog.cmake")
if (NOT GLOG_FOUND)
    message(FATAL_ERROR "glog not found, please install the glog library")
else ()
    include_directories(SYSTEM ${GLOG_INCLUDE_DIRS})
endif ()

# find gflags-------------------------------------------------------------------
Include("cmake/FindGFlags.cmake")
#find_package(gflags CONFIG REQUIRED)
if (NOT GFLAGS_FOUND)
    message(STATUS "gflags not found, build without gflags")
else ()
    include_directories(SYSTEM ${GFLAGS_INCLUDE_DIRS})
endif ()

# find yaml-------------------------------------------------------------------
Include("cmake/FindYaml.cmake")
if (NOT YAML_FOUND)
    message(STATUS "yaml not found, build without yaml")
else ()
    include_directories(SYSTEM ${YAML_INCLUDES})
endif ()

#openmp
find_package(OpenMP)
if(OPENMP_FOUND)
  message("OpenMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# find jemalloc-----------------------------------------------------------------
if (USE_JEMALLOC)
    include("cmake/FindJemalloc.cmake")
    if (NOT JEMALLOC_FOUND)
        message(STATUS "jemalloc not found, build without jemalloc")
    else ()
        add_definitions(-DUSE_JEMALLOC)
        include_directories(SYSTEM ${JEMALLOC_INCLUDE_DIRS})
    endif ()
endif ()

# find rdkafka---------------------------------------------------------------------
include("cmake/FindRdkafka.cmake")
if (NOT RDKAFKA_FOUND)
    message(STATUS "rdkafka not found, build without rdkafka")
endif ()

# ------------------------------------------------------------------------------
# generate libgrape-lite
# ------------------------------------------------------------------------------
file(GLOB_RECURSE CORE_SRC_FILES "include/libgrape-lite/grape/*.cc")

add_library(grape-lite ${CORE_SRC_FILES})

target_link_libraries(grape-lite ${MPI_CXX_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})

if (JEMALLOC_FOUND)
    target_link_libraries(grape-lite ${JEMALLOC_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable analytical_apps because gflags not found")
else ()
    add_executable(analytical_apps src/apps/flags.cc src/apps/run_app.cc)
    target_include_directories(analytical_apps PRIVATE src/apps)
    set_target_properties(analytical_apps PROPERTIES OUTPUT_NAME run_app)
    target_link_libraries(analytical_apps grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable gar_discover because gflags not found")
else ()
    add_executable(gar_discover src/apps/flags.cc src/apps/gar_discover.cc)
    target_include_directories(gar_discover PRIVATE src/apps)
    set_target_properties(gar_discover PROPERTIES OUTPUT_NAME gar_discover)
    target_link_libraries(gar_discover grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

if (NOT GFLAGS_FOUND)
    message(WARNING "Disable inc_gar_discover because gflags not found")
else ()
    add_executable(inc_gar_discover src/apps/flags.cc src/apps/inc_gar_discover.cc)
    target_include_directories(inc_gar_discover PRIVATE src/apps)
    set_target_properties(inc_gar_discover PROPERTIES OUTPUT_NAME inc_gar_discover)
    target_link_libraries(inc_gar_discover grape-lite ${MPI_CXX_LIBRARIES}
            ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} ${CMAKE_DL_LIBS} ${YAML_LIBRARIES})
endif ()

