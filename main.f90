program dynamicsOfOQS
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m
use system_hamiltonian_m
use rdftensor_m
use functions_m
use ordinary_equation_solution
use testes
use verlet_m 

implicit none

call call_ODE_solver(nm_divisoes)

end program DynamicsOfOQS
