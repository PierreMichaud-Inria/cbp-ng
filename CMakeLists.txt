cmake_minimum_required(VERSION 3.10)

project(cbp_ng LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Generate compile_commands.json and create a symlink to it
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
if(NOT ${CMAKE_BINARY_DIR} EQUAL ${PROJECT_SOURCE_DIR})
    file(RELATIVE_PATH relpath ${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/compile_commands.json)
    execute_process(COMMAND ln -sf ${relpath} ${PROJECT_SOURCE_DIR})
endif()

add_executable(cbp cbp.cpp)
target_link_libraries(cbp PRIVATE -lz)

add_executable(test_harcom test_harcom.cpp)
target_link_libraries(test_harcom PRIVATE -lz)
