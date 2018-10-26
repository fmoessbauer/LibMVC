gbench_dep = dependency('benchmark', required: false)

if not gbench_dep.found()
  gbench_proj = subproject('benchmark')
  gbench_dep  = gbench_proj.get_variable('gbenchmark_dep')
endif

# currently not supported
#openmp  = dependency('openmp')
threads = dependency('threads')

libmvcbench = executable(
  'LibMVC-bench',
  'main.cpp',
  include_directories: inc,
  cpp_args: '-fopenmp',
  link_args: '-fopenmp',
  dependencies: [threads, gbench_dep])

benchmark('reg-bench', libmvcbench, timeout: 60,
  args: ['../data/frb45-21-mis/frb45-21-1.mis', '900', '10'])

