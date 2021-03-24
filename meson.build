# File       : meson.build
# Created    : Wed Mar 24 2021  3:30:36 pm CET (+0100)
# Author     : Ivan Mihajlovic Milin
# Description: Main meson build definition
# Copyright 2021 ETH Zurich. All Rights Reserved.

project('Laplacian', ['cpp'], 
  version : '0.0.0',
  meson_version: '>=0.57.0',
  default_options : ['warning_level=3', 'cpp_std=c++14'],
  license: 'BSD-2-Clause',
)

# compiler object
cpp = meson.get_compiler('cpp')

## subproject CubismNova
# mpi dependency for CubismNova
mpi_dep = dependency('mpi', language: 'cpp', required: true)
# CubismNova dependency 
cubismnova_dep = dependency('CubismNova',
  default_options: ['CUBISM_TESTS=false'],
  fallback: ['CubismNova', 'cubismnova_dep'], 
  required: true, 
)
# directly from subprojects without system check
cubismnova_sub = subproject('CubismNova')
cubismnova_dep = cubismnova_sub.get_variable('cubismnova_dep')

## subproject googletest
# gtest subproject & dependencies
gtest_sub = subproject('gtest')
gtest_dep = gtest_sub.get_variable('gtest_dep')
gtest_main_dep = gtest_sub.get_variable('gtest_main_dep')
gmock_dep = gtest_sub.get_variable('gmock_main_dep')
# gtest for MPI: https://github.com/LLNL/gtest-mpi-listener
gtest_mpi_listener_sub = subproject('gtest-mpi-listener') 
gtest_mpi_listener_dep = gtest_mpi_listener_sub.get_variable('gtest_mpi_listener_dep')

# project sources 
subdir('include')
subdir('src')
subdir('test')
