cmake_minimum_required(VERSION 3.11-3.18)

project(main)

# add the -g debugging option to the compilation flags, so Dr. Memory can access line information.
# add_compile_options(-g)
# add_executable(main main.cpp)

#target_link_libraries(main PRIVATE cpp-library ${Rapidjson_LIBRARIES} rapidjson)
# target_link_libraries(main PRIVATE cpp-library)

# Tell C++ compiler to use C++20 features. We don't actually use any of them.
# target_compile_features(main PUBLIC cxx_std_20)