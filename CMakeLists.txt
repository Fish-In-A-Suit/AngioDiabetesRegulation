cmake_minimum_required(VERSION 3.11-3.18)

# set(CMAKE_FIND_DEBUG_MODE 1)

project(
    AngioDiabetesRegulation
    VERSION 0.1
    DESCRIPTION "AngioDiabetesRegulation Project"
    )

# std::filesystem:: not found error
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CUDA REQUIRED)

message("CMake Version: ${CMAKE_VERSION}")
message("cuda include dirs = ${CUDA_INCLUDE_DIRS}")
message("cuda libraries = ${CUDA_LIBRARIES}")

# need to have include directories so vscode can find vpckg libraries
include_directories(C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/include ${CUDA_INCLUDE_DIRS})

# find_package checks in C:\Program Files\CMake\share\cmake-3.25\Modules and tried to find FindRapidjson.cmake
# (I snatched a custom rapidjson.cmake that sets Rapidjson_INCLUDE_DIRS variable)
# Rapidjson_INCLUDE_DIRS is set to C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/include
# find_package(Rapidjson CONFIG REQUIRED)
# include_directories(${Rapidjson_INCLUDE_DIRS})  
# message("rapidjson include dirs = ${Rapidjson_INCLUDE_DIRS}")

add_subdirectory(cpp-library)   # look in "test-cpp-library" subdirectory for CMakeLists.txt to process
add_subdirectory(apps)    # look in "test-apps" subdirectory for CMakeLists.txt to process

add_executable(main apps/main.cpp)

message("ProjectDir = ${ProjectDir}")
message("SolutionDir = ${SolutionDir}")
message("WorkspaceRoot = ${workspaceRoot}")
message("CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
set_property(TARGET main PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
set(PROJECT_SOURCE_DIR "$(CMAKE_CURRENT_SOURCE_DIR)")

# to enable Visual Studio "Edit and continue" option
if(MSVC)
  target_compile_options(main PUBLIC "/ZI")
  target_link_options(main PUBLIC "/INCREMENTAL")
endif()

add_compile_options("/std:c++17")

target_link_libraries(main PRIVATE cpp-library ${CUDA_LIBRARIES} stdc++fs)

target_compile_features(main PUBLIC cxx_std_20)