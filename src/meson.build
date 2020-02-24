# library
xavier_cpp_sources = files([
  'aligner.cpp',
  'score.cpp',
  'seed.cpp',
  'trace.cpp',
  'xavier.cpp',
])

# install library if
# - either running as a proper project
# - or using shared libraries
xavier_lib_install = (not meson.is_subproject()) or (get_option('default_library') != 'static')

simd_mod = import('unstable_simd')
xavier_simd_libs = []

if get_option('avx2')
  xavier_simd_libs += simd_mod.check(
    'vectors_avx2',
    avx2 : files([
      'vectors.cpp',
    ]),
    include_directories : xavier_include_directories,
    c_args : [xavier_flags],
    compiler : cpp)[0]
endif

if get_option('sse42')
  xavier_simd_libs += simd_mod.check(
    'vectors_sse42',
    sse42 : files([
      'vectors.cpp',
    ]),
    include_directories : xavier_include_directories,
    c_args : [xavier_flags],
    compiler : cpp)[0]
endif

# this is a crutch to unify all compiled objects into
# one library, in order to avoid tons of libraries
# showing up in the pkg-config file
# https://github.com/mesonbuild/meson/issues/2478#issuecomment-351132872
xavier_simd_all_objects = []
foreach i : xavier_simd_libs
  xavier_simd_all_objects += i.extract_all_objects()
endforeach

xavier_lib = library(
  'xavier',
  xavier_cpp_sources,
  install : xavier_lib_install,
  # use boost SONAME practice:
  #   cause ld startup issues before
  #   you even have the chance of running
  #   into ABI issues.
  soversion : meson.project_version(),
  version : meson.project_version(),
  objects : xavier_simd_all_objects,
  include_directories : xavier_include_directories,
  cpp_args : xavier_flags)
