set(SOURCES main SolverTests InterfaceTests)

add_executable("${PROJECT_NAME}-test" ${SOURCES})
set_target_properties("${PROJECT_NAME}-test" PROPERTIES CXX_STANDARD 14)

target_link_libraries("${PROJECT_NAME}-test" gtest)

