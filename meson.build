project('LibMVC', 'cpp',
  version : '1.0.0',
  license : 'MIT',
  default_options: [
    'cpp_std=c++14',
    'buildtype=release'])

inc = include_directories('LibMVC')
install_subdir('LibMVC', install_dir: 'include')
# install tiny graph for testing
install_subdir(
  'data/frb30-15-mis',
  install_dir: 'data')

subdir('standalone')
subdir('test')
subdir('bench')

doxygen = find_program('doxygen', required : false)
if doxygen.found()
  message('build documentation')
  cdata = configuration_data()
  cdata.set('VERSION', meson.project_version())
  subdir('doc')
endif

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  name: 'libmvc',
  version: meson.project_version(),
  filebase: 'libmvc',
  description: 'LibMVC is a collection of fast iterative minimum vertex cover solvers')

