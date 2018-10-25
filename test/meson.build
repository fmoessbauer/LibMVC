# create gtest dependency

gtest_dep = dependency('gtest', main: false, required: false)

if not gtest_dep.found()
  gtest_proj = subproject('gtest')
  gtest_dep  = gtest_proj.get_variable('gtest_dep')
endif

testing = executable(
  'LibMVC-test',
  'main.cpp', 'SolverTests.cpp', 'InterfaceTests.cpp',
  include_directories: inc,
  dependencies: gtest_dep)

test('all-units', testing,
  args: ['../data/frb30-15-mis/frb30-15-1.mis', '420', '450', '10'])

