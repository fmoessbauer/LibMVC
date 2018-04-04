set(SOURCES "main")
set(EXEC_NAME "${PROJECT_NAME}-fastvc-standalone")

set(RELEASE_COMPILE_FLAGS "-O3 -Wall")
set(DEBUG_COMPILE_FLAGS "-g -Wall -Werror") 

add_executable(${EXEC_NAME} ${SOURCES})
set_target_properties(${EXEC_NAME} PROPERTIES CXX_STANDARD 14)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS "${DEBUG_COMPILE_FLAGS}")
else()
  set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS "${RELEASE_COMPILE_FLAGS}")
endif()

