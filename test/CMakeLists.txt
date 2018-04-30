set(SOURCES main SolverTests InterfaceTests)

add_executable("${PROJECT_NAME}-test" ${SOURCES})
set_target_properties("${PROJECT_NAME}-test" PROPERTIES CXX_STANDARD 14)

target_link_libraries("${PROJECT_NAME}-test" gtest)

if(${LIBMVC_ENABLE_COVERAGE})
  message("Setup coverage tests for ${PROJECT_NAME}-test")
  set(COVERAGE_EXCLUDES "/usr/\\*" "${CMAKE_SOURCE_DIR}/vendor/\\*" "\\*main.cpp")
  setup_target_for_coverage(
    NAME "${PROJECT_NAME}-test-coverage"
    EXECUTABLE ./test/${PROJECT_NAME}-test ${CMAKE_BINARY_DIR}/../data/frb30-15-mis/frb30-15-1.mis 420 450 10
    DEPENDENCIES "${PROJECT_NAME}-test")
endif()

