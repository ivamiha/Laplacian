# File       : meson.build  
# Created    : Wed Mar 24 2021  3:36:42 pm CET (+0100)
# Author     : Ivan Mihajlovic Milin
# Description: Meson executable build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

# executable for non-MPI files
executable('main', 
  files('main.cpp'),
  cpp_args: '-fopenmp',
  link_args: '-fopenmp',
  dependencies: [cubismnova_dep, mpi_dep, omp_dep],
  include_directories: inc, 
  install: false, 
)
