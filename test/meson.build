# File       : meson.build
# Created    : Wed Mar 24 2021  4:02:52 pm CET (+0100) 
# Author     : Ivan Mihajlovic Milin
# Description: Meson unit test build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

e = executable('TestLaplacian', 
  files('TestLaplacian.cpp'), 
  dependencies: [cubismnova_dep, mpi_dep, gtest_dep],
  include_directories: inc, 
  install: false, 
)
test('TestLaplacian', e, 
  workdir: '/tmp', 
  protocol: 'gtest' 
)
