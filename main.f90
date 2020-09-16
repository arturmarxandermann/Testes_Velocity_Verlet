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


print*, "EVOLUCAO TERMINOU, FAZENDO ARQUIVOS RELACIONADOS &
COM A FIGURA DOS SITIOS"

print*, "ARQUIVOS DE SAÍDA:"
print*, "energiael -> energia do eletron Tr(H RHO);"
print*, "energiazero -> energia do eletron no tempo inicial;"
print*, "energiacinetica -> energia cinetica classica;"
print*, "energiapotencial -> energia potencial classica;"
print*, "energiatotal -> energiazero, energiael + energiacinetica + energiapotencial;"
print*, "popSiteBasis -> populacao na base local;"
print*, "popNonSiteBasis -> populacao na base dos autoestados do hamiltoniano;"
print*, "radius -> raio do sitio 1, raio do sitio 2, ... ;"
print*, "forces -> força no sitio 1, força no sitio2, ... "











13 format(3es14.5E3)

end program DynamicsOfOQS
