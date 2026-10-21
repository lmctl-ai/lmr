cmake_minimum_required(VERSION 3.22)
project(lmr VERSION 0.2.0 LANGUAGES CXX)
# The native rcli: a C++17 statistics core (standard library only) and the `_dsl_` bytecode host lmr.
# No R, no Python, no Docker, no third-party libraries. The build is one executable: on Linux it links libstdc++ and
# libgcc statically and glibc (libc, libm) dynamically (-DLMR_FULL_STATIC=ON links glibc too); on Windows (MinGW)
# it links fully static by default, so no MinGW DLLs are needed; macOS links only its system libraries.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(WIN32)
  set(LMR_FULL_STATIC_DEFAULT ON)
else()
  set(LMR_FULL_STATIC_DEFAULT OFF)
endif()
option(LMR_FULL_STATIC "Link lmr fully statically (the C library included)" ${LMR_FULL_STATIC_DEFAULT})
# Sanitizer builds for tests: -DLMR_SANITIZE=address,undefined or -DLMR_SANITIZE=thread (implies a dynamic link).
set(LMR_SANITIZE "" CACHE STRING "Comma-separated -fsanitize= list for every target (empty: none)")
if(LMR_SANITIZE)
  add_compile_options(-fsanitize=${LMR_SANITIZE} -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=${LMR_SANITIZE})
  set(LMR_FULL_STATIC OFF CACHE BOOL "" FORCE)
endif()
find_package(Threads REQUIRED)

# Embed the exact schema bytes: the bytecode is bound to their SHA-256.
set(LMR_SCHEMA "${CMAKE_CURRENT_SOURCE_DIR}/schema/lmr.dsl.json")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${LMR_SCHEMA}")
file(READ "${LMR_SCHEMA}" LMR_SCHEMA_HEX HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," LMR_SCHEMA_INITIALIZER "${LMR_SCHEMA_HEX}")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp.tmp"
  "#pragma once\n#include <cstddef>\nnamespace rcli { inline constexpr unsigned char embedded_schema[] = {${LMR_SCHEMA_INITIALIZER}}; inline constexpr std::size_t embedded_schema_size = sizeof(embedded_schema); }\n")
file(COPY_FILE "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp.tmp"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp" ONLY_IF_DIFFERENT)

set(LMR_WARNINGS -Wall -Wextra -Wpedantic)
# No fused multiply-add contraction: results must not depend on the target CPU (the DSL spec requires the same).
add_compile_options(-ffp-contract=off)
# M_PI and friends are not ISO C++; MinGW hides them under -std=c++17 without this.
add_compile_definitions(_USE_MATH_DEFINES)
file(GLOB LMR_NUM_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/num/*.cpp")
add_library(rcli_num STATIC ${LMR_NUM_SOURCES})
target_include_directories(rcli_num PUBLIC include)
target_compile_options(rcli_num PRIVATE ${LMR_WARNINGS})
target_link_libraries(rcli_num PUBLIC Threads::Threads)

file(GLOB LMR_HOST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/host/*.cpp")
add_executable(lmr ${LMR_HOST_SOURCES})
target_include_directories(lmr PRIVATE cpp src/host "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_compile_options(lmr PRIVATE ${LMR_WARNINGS})
target_compile_definitions(lmr PRIVATE LMR_VERSION="${PROJECT_VERSION}")
target_link_libraries(lmr PRIVATE rcli_num)
set_target_properties(lmr PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# "Clang" also matches AppleClang; macOS has no -static or -static-libstdc++.
if((CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang") AND NOT APPLE AND NOT LMR_SANITIZE)
  target_link_options(lmr PRIVATE -static-libstdc++ -static-libgcc)
  if(LMR_FULL_STATIC)
    target_link_options(lmr PRIVATE -static)
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  file(GLOB LMR_UNIT_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp")
  foreach(test_source ${LMR_UNIT_TESTS})
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable(${test_name} "${test_source}")
    target_include_directories(${test_name} PRIVATE tests)
    target_compile_options(${test_name} PRIVATE ${LMR_WARNINGS})
    target_link_libraries(${test_name} PRIVATE rcli_num)
    add_test(NAME unit_${test_name} COMMAND ${test_name})
  endforeach()
  # Through bash explicitly: Windows cannot execute a .sh file directly.
  find_program(LMR_BASH bash REQUIRED)
  add_test(NAME schema_matches_registry
    COMMAND "${LMR_BASH}" "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check-schema.sh" $<TARGET_FILE:lmr>)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    add_test(NAME host_bindings
      COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/tests/host_test.py" $<TARGET_FILE:lmr>)
    set_tests_properties(host_bindings PROPERTIES TIMEOUT 900)
  endif()
  add_test(NAME dsl_examples
    COMMAND "${LMR_BASH}" "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run-examples.sh" $<TARGET_FILE:lmr>)
  set_tests_properties(dsl_examples PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.." TIMEOUT 600)
endif()
