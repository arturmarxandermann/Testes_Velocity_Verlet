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
use time_evolution_m
use verlet_m 

implicit none


    call System_Dynamics(nm_divisoes)

end program DynamicsOfOQS
