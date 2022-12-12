# version 3.11 or later of CMake needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(
    AngioDiabetesRegulation
    VERSION 0.1
    DESCRIPTION "Small CMake example built with VSCode")

# find_package(jsoncpp REQUIRED)
# message("FOUND CONFIG: ${jsoncpp_CONFIG}")
# message("INCLUDE DIRS: ${jsoncpp_INCLUDE_DIRS}")
# message("LIBRARIES: ${jsoncpp_LIBRARIES}")
# message("INTERFACE INCLUDE LIBS: ${jsoncpp_INTERFACE_INCLUDE_DIRECTORIES}")

# need to have include directories so vscode can find vpckg libraries
INCLUDE_DIRECTORIES(C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/include)
find_package(RapidJSON REQUIRED)

add_subdirectory(cpp-library)   # look in "test-cpp-library" subdirectory for CMakeLists.txt to process
add_subdirectory(apps)    # look in "test-apps" subdirectory for CMakeLists.txt to process