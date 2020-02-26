
# find GoogleTest and GoogleMock
xavier_gmock_main_dep = dependency('gmock_main', fallback : ['gtest', 'gmock_main_dep'])

xavier_test_cpp_sources = files([
  'test.cpp',
])


xavier_test = executable(
  'xavier test',
  xavier_test_cpp_sources,
  dependencies : xavier_gmock_main_dep,
  include_directories : xavier_include_directories,
  link_with : xavier_lib,
  cpp_args : [xavier_flags, '-mavx2'],
  install : false)

test(
  'xavier gtest unittests',
  xavier_test,
  args : [
    '--gtest_output=xml:' + join_paths(meson.build_root(), 'gtest-unittests.xml')],
  env : [
    'REPODIR=' + join_paths([meson.current_source_dir(), '..']),
    'ARGS=-V',
    'VERBOSE=1'],
  timeout : 1800)
