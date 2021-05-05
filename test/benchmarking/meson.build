# File       : meson.build  
# Created    : Wed May 05 2021 03:44:06 PM (+0200)
# Author     : Ivan Mihajlovic Milin
# Description: Meson unit test build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('benchmarking', 
  files('benchmarking.cpp'),
  dependencies:  [cubismnova_dep, mpi_dep, gtest_dep], 
  include_directories: inc, 
  install: false
)
test('benchmarkingTest', e, 
  workdir: '/tmp',
  protocol: 'gtest'
)
