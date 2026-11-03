# Installs an aerobus build tree into a scratch prefix, then configures, builds and runs
# the consumer project next to this script, which uses find_package(Aerobus).
#
# Usage: cmake -DAEROBUS_BINARY_DIR=... -DWORK_DIR=... -DCXX_COMPILER=... [-DCONFIG=...] -P check_package.cmake
foreach(var AEROBUS_BINARY_DIR WORK_DIR CXX_COMPILER)
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "${var} is not set")
  endif()
endforeach()
if(NOT CONFIG)
  set(CONFIG Release)
endif()

set(prefix "${WORK_DIR}/prefix")
set(consumer_build "${WORK_DIR}/build")
file(REMOVE_RECURSE "${WORK_DIR}")

function(run_step name)
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT rc EQUAL 0)
    message(FATAL_ERROR "${name} failed (${rc})\n${out}\n${err}")
  endif()
endfunction()

run_step(install "${CMAKE_COMMAND}" --install "${AEROBUS_BINARY_DIR}" --prefix "${prefix}" --config "${CONFIG}")
run_step(configure "${CMAKE_COMMAND}" -S "${CMAKE_CURRENT_LIST_DIR}" -B "${consumer_build}"
         "-DCMAKE_PREFIX_PATH=${prefix}" "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" "-DCMAKE_BUILD_TYPE=${CONFIG}")
run_step(build "${CMAKE_COMMAND}" --build "${consumer_build}" --config "${CONFIG}")
run_step(run "${consumer_build}/bin/consumer")
message(STATUS "find_package(Aerobus) consumer built and ran")
