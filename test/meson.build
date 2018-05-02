# create gtest dependency
gtest_dir = '../vendor/google/googletest/googletest/'
gtest_incdir = include_directories(join_paths(gtest_dir, 'include'))

libgtest = static_library('gtest',
  include_directories : [include_directories(gtest_dir), gtest_incdir],
  sources : [
    join_paths(gtest_dir, 'src', 'gtest-all.cc'),
  ]
)
gtest_dep = declare_dependency(dependencies : dependency('threads'),
include_directories : gtest_incdir,
link_with : libgtest)

testing = executable('LibMVC-test', 'main.cpp', 'SolverTests.cpp', 'InterfaceTests.cpp', dependencies: gtest_dep)

test('all-units', testing,
  args: ['../data/frb30-15-mis/frb30-15-1.mis', '420', '450', '10'])

