# Laplacian
Semester thesis focusing on the development of optimized Laplacian compute 
kernels with application to the Gray-Scott reaction-diffusion system. The 
project works with the `CubismNova` framework; a C++ template library used for
solving partial differential equations on structured uniform or stretched 
grids. The developed compute kernels are optimized with the addition of 
data-level parallelism (DLP) with `ISPC` and thread-level parallelism (TLP)
with `OMP`. 

## Project checklist 
The date within brackets `[dd/mm/yy]` denotes the due date originally assigned 
during the planning of the project. A color code has been utilized to denote 
the "relevance" of the task <br /> 
🟢 task necessary for project completion, <br />
🟡 task optional, but would ideally like to implement.

**March**
* [05/03/21] 🟢 literature research on reaction-diffusion systems           ✔️
* [12/03/21] 🟢 familiarization with CubismNova                             ✔️
* [26/03/21] 🟢 implementation of naive FD kernels                          ✔️

**April**
* [02/04/21] 🟢 implementation of order verification study                  ✔️
* [09/04/21] 🟢 literature research on data-level parallelism               ✔️ 
* [09/04/21] 🟢 literature research & implementation of benchmark suite     ✔️
* [16/04/21] 🟢 initial implementation of sliced ISPC kernel 
* [16/04/21] 🟢 performance analysis of naive FD kernels
* [23/04/21] 🟢 optimization of sliced ISPC kernel 
* [23/04/21] 🟢 performance analysis of sliced ISPC kernel
* [30/04/21] 🟡 extension of sliced ISPC kernels to 4th-order FDs
* [30/04/21] 🟢 finalizing serial Gray-Scott solver  

**May**
* [07/05/21] 🟢 extension to non-zero von Neumann boundaries 
* [14/05/21] 🟢 extension of Gray-Scott solver to TLP 
* [14/05/21] 🟢 performance analysis of shared memory code 
* [28/05/21] 🟡 solution verification for Gray-Scott solver
* [28/05/21] 🟡 Gray-Scott solver developmental biology applications
* [28/05/21] 🟢 3D volume rendering og Gray-Scott problem for movie

**June**
* [04/06/21] 🟢 write-up thesis & create presentation of work

Note: Gray-Scott solver-relevant development is undertaken in the 
`Gray-Scott-3D` repository.

## Current focus tasks
Following tasks are being undertaken this week: 
* initial implementation of sliced ISPC kernel 
* performance analysis of naive FD kernels

Further notes: <br />
To be completed as develop ideas throughout the week. 
Use as meeting reference. 
