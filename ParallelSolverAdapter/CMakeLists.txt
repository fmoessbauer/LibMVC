set(SOURCES "main")

set(EXEC_NAME "${PROJECT_NAME}-parallel-numvc-standalone")

set(RELEASE_COMPILE_FLAGS "-O3 -Wall -DALGO_NUMVC")
set(DEBUG_COMPILE_FLAGS "-g -Wall -Werror -DALGO_NUMVC")

find_package (Threads REQUIRED)

add_executable(${EXEC_NAME} ${SOURCES})
target_link_libraries(${EXEC_NAME} Threads::Threads)
set_target_properties(${EXEC_NAME} PROPERTIES CXX_STANDARD 14)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS "${DEBUG_COMPILE_FLAGS}")
else()
  set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS "${RELEASE_COMPILE_FLAGS}")
endif()

