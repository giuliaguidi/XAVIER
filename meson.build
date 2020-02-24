project(
  'xavier',
  ['c','cpp'],
  version : '0.0.1',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.46.0')

# CXXFLAGS
xavier_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
# to remove warnings like
#  '-Wno-unused-command-line-argument',
  ]
  if cpp.has_argument(cflag)
    xavier_flags += cflag
  endif
endforeach

# headers
subdir('include')

# sources + library
subdir('src')

# demo
subdir('examples')

# tests
if (not meson.is_subproject()) and get_option('tests')
  subdir('test')
endif

# dependency info
if (not meson.is_subproject())
  import('pkgconfig').generate(
    xavier_lib,
    version : meson.project_version(),
    name : 'xavier',
    filebase : 'xavier',
    description : 'Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment')
endif

xavier_dep = declare_dependency(
  include_directories : xavier_include_directories,
  link_with : xavier_lib,
  version : meson.project_version())
