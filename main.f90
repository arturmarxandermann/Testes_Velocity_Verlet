program dynamicsOfOQS
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m
use system_hamiltonian_m
use system_operators_m
use rdftensor_m
use functions_m
use ordinary_equation_solution
use plot_auto_function_m
use gnugraphs_m
!use edo_solver_m
use testes
use verlet_m 
!use rkf45_m 

implicit none


call create_FT_func

call call_ODE_solver(nm_divisoes)


print*, "EVOLUCAO TERMINOU, FAZENDO ARQUIVOS RELACIONADOS &
COM A FIGURA DOS SITIOS"

!call calculate_fem_autofunction(sites_array)











13 format(3es14.5E3)

end program DynamicsOfOQS
