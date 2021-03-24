project('Laplacian', ['cpp'], 
  version : '0.0.0',
  meson_version: '>=0.57.0',
  default_options : ['warning_level=3', 'cpp_std=c++14'],
  license: 'BSD-2-Clause',
)

mpi_dep = dependency('mpi', language: 'cpp', required: true)

cubismnova_dep = dependency('CubismNova',
  default_options: ['CUBISM_TESTS=false'],
  fallback: ['CubismNova', 'cubismnova_dep'], 
  required: true, 
)
# directly from subprojects without system check
cubismnova_sub = subproject('CubismNova')
cubismnova_dep = cubismnova_sub.get_variable('cubismnova_dep')

subdir('include')
subdir('src')
