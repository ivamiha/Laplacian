# File       : meson.build
# Created    : Wed May 26 2021 08:57:12 AM (+0200)
# Author     : Ivan Mihajlovic Milin
# Description: Meson unit test build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('benchmarking',
  files('NeumannTest.cpp'),
  dependencies: [cubismnova_dep, mpi_dep, gtest_dep],
  include_directories: inc, 
  install: false
)
test('NeumannTest', e,
  workdir: '/tmp',
  protocol: 'gtest'
)
