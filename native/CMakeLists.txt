cmake_minimum_required(VERSION 3.22)
project(rcli_dsl VERSION 0.1.0 LANGUAGES CXX)
# The native rcli: a C++17 statistics core (standard library only) and the `_dsl_` bytecode host rcli-dsl.
# No R, no Python, no Docker, no third-party libraries. The default build is one executable that links libstdc++ and
# libgcc statically and glibc (libc, libm) dynamically; -DRCLI_DSL_FULL_STATIC=ON links glibc statically too.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(RCLI_DSL_FULL_STATIC "Link rcli-dsl fully statically, glibc included (Linux)" OFF)
# Sanitizer builds for tests: -DRCLI_SANITIZE=address,undefined or -DRCLI_SANITIZE=thread (implies a dynamic link).
set(RCLI_SANITIZE "" CACHE STRING "Comma-separated -fsanitize= list for every target (empty: none)")
if(RCLI_SANITIZE)
  add_compile_options(-fsanitize=${RCLI_SANITIZE} -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=${RCLI_SANITIZE})
  set(RCLI_DSL_FULL_STATIC OFF CACHE BOOL "" FORCE)
endif()
find_package(Threads REQUIRED)

# Embed the exact schema bytes: the bytecode is bound to their SHA-256.
set(RCLI_DSL_SCHEMA "${CMAKE_CURRENT_SOURCE_DIR}/schema/rcli.dsl.json")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${RCLI_DSL_SCHEMA}")
file(READ "${RCLI_DSL_SCHEMA}" RCLI_DSL_SCHEMA_HEX HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," RCLI_DSL_SCHEMA_INITIALIZER "${RCLI_DSL_SCHEMA_HEX}")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp.tmp"
  "#pragma once\n#include <cstddef>\nnamespace rcli { inline constexpr unsigned char embedded_schema[] = {${RCLI_DSL_SCHEMA_INITIALIZER}}; inline constexpr std::size_t embedded_schema_size = sizeof(embedded_schema); }\n")
file(COPY_FILE "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp.tmp"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/schema_embedded.hpp" ONLY_IF_DIFFERENT)

set(RCLI_WARNINGS -Wall -Wextra -Wpedantic)
# No fused multiply-add contraction: results must not depend on the target CPU (the DSL spec requires the same).
add_compile_options(-ffp-contract=off)
file(GLOB RCLI_NUM_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/num/*.cpp")
add_library(rcli_num STATIC ${RCLI_NUM_SOURCES})
target_include_directories(rcli_num PUBLIC include)
target_compile_options(rcli_num PRIVATE ${RCLI_WARNINGS})
target_link_libraries(rcli_num PUBLIC Threads::Threads)

file(GLOB RCLI_HOST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/host/*.cpp")
add_executable(rcli-dsl ${RCLI_HOST_SOURCES})
target_include_directories(rcli-dsl PRIVATE cpp src/host "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_compile_options(rcli-dsl PRIVATE ${RCLI_WARNINGS})
target_compile_definitions(rcli-dsl PRIVATE RCLI_DSL_VERSION="${PROJECT_VERSION}")
target_link_libraries(rcli-dsl PRIVATE rcli_num)
set_target_properties(rcli-dsl PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
if((CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang") AND NOT APPLE AND NOT RCLI_SANITIZE)
  target_link_options(rcli-dsl PRIVATE -static-libstdc++ -static-libgcc)
  if(RCLI_DSL_FULL_STATIC)
    target_link_options(rcli-dsl PRIVATE -static)
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  file(GLOB RCLI_UNIT_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp")
  foreach(test_source ${RCLI_UNIT_TESTS})
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable(${test_name} "${test_source}")
    target_include_directories(${test_name} PRIVATE tests)
    target_compile_options(${test_name} PRIVATE ${RCLI_WARNINGS})
    target_link_libraries(${test_name} PRIVATE rcli_num)
    add_test(NAME unit_${test_name} COMMAND ${test_name})
  endforeach()
  add_test(NAME schema_matches_registry
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/scripts/check-schema.sh" $<TARGET_FILE:rcli-dsl>)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    add_test(NAME host_bindings
      COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/tests/host_test.py" $<TARGET_FILE:rcli-dsl>)
    set_tests_properties(host_bindings PROPERTIES TIMEOUT 900)
  endif()
  add_test(NAME dsl_examples
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/scripts/run-examples.sh" $<TARGET_FILE:rcli-dsl>)
  set_tests_properties(dsl_examples PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.." TIMEOUT 600)
endif()
