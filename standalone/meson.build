solver = ['NuMVC', 'FastVC']

threads_dep = dependency('threads')

foreach s : solver
executable(
  s + '-standalone',
  'main.cpp',
  cpp_args : '-DSOLVER=' + s,
  include_directories: inc,
  install: true,
  install_dir: 'bin')
endforeach

# parallel versions
foreach s : solver
executable(
  s + '-parallel-standalone',
  'main.cpp',
  cpp_args : ['-DSOLVER=' + s, '-DPARALLEL'],
  include_directories: inc,
  install: true,
  install_dir: 'bin',
  dependencies: threads_dep)
endforeach

