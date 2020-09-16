module parameters_m
use types_m
use constants_m 
!modulo que tem todos os parametros que usamos apenas no main.f90
!TODOS OS TEMPOS ESTAO EM UNIDADES DE PS E AS FREQUENCIAS EM UNIDADE DE THZ

real*8, parameter :: tmin = 0.d0  !tempo inicial
real*8, parameter :: tmax = 1.2d0  !tempo final
integer, parameter :: nsites = 2  !numero de sítios que entram na simulacao
integer, parameter :: nmstates_el = 3 !numero de estados para o eletron
integer, parameter :: dim_el = (nsites*nmstates_el)
integer, parameter :: nm_rows = 1 !numero de linhas para os sitios
integer, parameter :: nm_columns = 2 !numero de colunas para os sitios
real*8, parameter :: distance_molecule = 0.25d-9  !distancia entre os sitios
integer, parameter :: neqn_el = dim_el * dim_el !numero de equacoes diferenciais acopladas
integer, parameter :: nm_divisoes = 3000 !numero de passos utilizado no programa
integer, parameter :: initState = 1


real*8, parameter :: siteMass = 1.d-27 !massa dos sitios
real*8, parameter :: raioZero = 0.5d-9 !raio inicial dos sitios
real*8, parameter :: effectiveMassMult = 5.d0 !fator de multiplicacao da massa efetiva do eletron
real*8, parameter :: siteCoupling = 0.6d0 != tc


real*8, parameter :: multFactorDer =  1.d0 !fator que estou multiplicando os elementos 0010 e 1000 da matriz derivada
real*8, parameter :: multFactorOvlp = 1.d0 !fator que estou multiplicando os elementos 0010 e 1000 da matriz de overlap


real*8, parameter :: me = effectiveMassMult * electronMass !massa efetiva do eletron





integer :: neqn
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_ham_in
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: frequency_matrix
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: hamiltoniano, RWMatrix
real*8 :: delta_t

real*8 :: energiazeroel
type(quantum_site), allocatable ::  sites_array(:, :)


!integer, parameter :: ndim = ((dimensao*(dimensao + 1)) / 2) !numero de eq. pra resolver na matriz rho_real


end module parameters_m
