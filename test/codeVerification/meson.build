# File       : meson.build
# Created    : Wed Mar 24 2021  4:02:52 pm CET (+0100) 
# Author     : Ivan Mihajlovic Milin
# Description: Meson unit test build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

test_data = custom_target('code-verification-plot',
  command: ['rsync', '@INPUT@', '@OUTPUT@'],
  input: files('plot.py'),
  output: 'plot.py',
)
e = executable('codeVerification', 
  [files([
    'codeVerification.cpp',
    ]), test_data], 
  dependencies: [cubismnova_dep, mpi_dep, gtest_dep],
  include_directories: inc, 
  install: false, 
)
test('codeVerificationTest', e, 
  workdir: meson.current_build_dir(), 
  protocol: 'gtest' 
)
