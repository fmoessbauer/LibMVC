# create google benchmark dependency
gbench_dir    = '../vendor/google/benchmark'
gbench_incdir = include_directories(join_paths(gbench_dir, 'include'))
gbench_srcdir = join_paths(gbench_dir, 'src')

libgbench = static_library('gbench',
  include_directories : [include_directories(gbench_dir), gbench_incdir],
  sources : [
    join_paths(gbench_srcdir, 'benchmark.cc'),
    join_paths(gbench_srcdir, 'benchmark_register.cc'),
    join_paths(gbench_srcdir, 'colorprint.cc'),
    join_paths(gbench_srcdir, 'commandlineflags.cc'),
    join_paths(gbench_srcdir, 'complexity.cc'),
    join_paths(gbench_srcdir, 'console_reporter.cc'),
    join_paths(gbench_srcdir, 'counter.cc'),
    join_paths(gbench_srcdir, 'csv_reporter.cc'),
    join_paths(gbench_srcdir, 'json_reporter.cc'),
    join_paths(gbench_srcdir, 'reporter.cc'),
    join_paths(gbench_srcdir, 'sleep.cc'),
    join_paths(gbench_srcdir, 'statistics.cc'),
    join_paths(gbench_srcdir, 'string_util.cc'),
    join_paths(gbench_srcdir, 'sysinfo.cc'),
    join_paths(gbench_srcdir, 'timers.cc')
  ]
)
gbench_dep = declare_dependency(dependencies : dependency('threads'),
  include_directories : gbench_incdir,
  link_with : libgbench)

# currently not supported
#openmp  = dependency('openmp')
threads = dependency('threads')

libmvcbench = executable(
  'LibMVC-bench',
  'main.cpp',
  cpp_args: '-fopenmp',
  link_args: '-fopenmp',
  dependencies: [threads, gbench_dep])

benchmark('reg-bench', libmvcbench,
  args: ['../data/frb45-21-mis/frb45-21-1.mis', '900', '10'])

