cmake_minimum_required (VERSION 3.5)

project (graphics17a_hw03)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile flags
if(APPLE)
    add_definitions(-Wno-missing-braces)
endif(APPLE)
if(WIN32)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
endif(WIN32)

add_definitions(-DYOBJ_NO_IMAGE -DYGLTF_NO_IMAGE -DYSCN_NO_IMAGE)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_library(hwlib STATIC
    src/image.h src/image.cpp
    src/scene.h src/scene.cpp
    src/vmath.h
    src/ext/yocto_math.h
    src/ext/yocto_utils.h
    src/ext/yocto_scn.h src/ext/yocto_scn.cpp
    src/ext/yocto_obj.h src/ext/yocto_obj.cpp
    src/ext/yocto_gltf.h src/ext/yocto_gltf.cpp
    src/ext/stb_image.h src/ext/stb_image_write.h)

add_executable(raytrace src/raytrace.cpp)
target_link_libraries(raytrace hwlib)
target_link_libraries(raytrace Threads::Threads)

add_executable(model src/model.cpp)
target_link_libraries(model hwlib)
