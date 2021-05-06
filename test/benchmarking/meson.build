# File       : meson.build  
# Created    : Wed May 05 2021 03:44:06 PM (+0200)
# Author     : Ivan Mihajlovic Milin
# Description: Meson unit test build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

test_data = custom_target('benchmarking-plot', 
  command: ['rsync', '@INPUT@', '@OUTPUT@'], 
  input: files('plot.py'), 
  output: 'plot.py',
)
e = executable('benchmarking', 
  [files(['benchmarking.cpp',
    ]), test_data],
  cpp_args: '-fopenmp',
  link_args: '-fopenmp',
  dependencies:  [cubismnova_dep, mpi_dep, omp_dep, gtest_dep], 
  include_directories: inc, 
  install: false
)
test('benchmarkingTest', e, 
  workdir: '/tmp',
  protocol: 'gtest'
)
