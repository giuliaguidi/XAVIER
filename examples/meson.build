xavier_demo = executable(
  'demo', files([
    'demo.cpp']),
  install : not meson.is_subproject(),
  include_directories : xavier_include_directories,
  link_with : xavier_lib,
  c_args : xavier_flags)
