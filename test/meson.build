# create gtest dependency

gtest_dep = dependency('gtest', main: false, required: false)

testing = executable(
  'LibMVC-test',
  'main.cpp', 'SolverTests.cpp', 'InterfaceTests.cpp',
  include_directories: inc,
  dependencies: gtest_dep)

test('all-units', testing, timeout: 120,
  args: ['../data/frb30-15-mis/frb30-15-1.mis', '420', '450', '10'])

