module parameters_m
use types_m
use constants_m 
!modulo que tem todos os parametros que usamos apenas no main.f90
!TODOS OS TEMPOS ESTAO EM UNIDADES DE PS E AS FREQUENCIAS EM UNIDADE DE THZ

    real*8,  parameter :: tmin = 0.d0                  !tempo inicial
    real*8,  parameter :: tmax = 5.d0                  !tempo final
    integer, parameter :: nsites = 2                   !numero de sítios que entram na simulacao
    integer, parameter :: ns_el = 3                    !numero de estados para o eletron
    integer, parameter :: d_el = (nsites*ns_el)        !tamanho da matriz rho, hamiltoniana
    integer, parameter :: nr = 1                       !numero de linhas para os sitios
    integer, parameter :: nc = 2                       !numero de colunas para os sitios
    real*8,  parameter :: distance_molecule = 0.25d-9  !distancia entre os sitios
    integer, parameter :: neqn_el = d_el * d_el        !numero de equacoes diferenciais acopladas
    integer, parameter :: nm_divisoes = 100000         !numero de passos utilizado no programa
    integer, parameter :: initState = 1                !estado inicial do eletron
    real*8,  parameter :: siteMass = 1.d-27            !massa dos sitios
    real*8,  parameter :: raioZero = 0.5d-9            !raio inicial dos sitios
    real*8,  parameter :: velZero  = 0.d0              !2.d3
    real*8,  parameter :: effectiveMassMult = 5.d0     !fator de multiplicacao da massa efetiva do eletron
    real*8,  parameter :: siteCoupling = 0.6d0         != tc


    real*8, parameter :: bathTemp = 300.d0             !temperatura fixa em kelvin
    real*8, parameter :: bathCoup = 1.d-15             !constante de acoplamento Berendsen em segundos
    

    real*8             :: Qn_erg(6)
    INTEGER            :: Qn(6)
    
    real*8,  parameter :: me = effectiveMassMult * electronMass !massa efetiva do eletron
    real*8                          :: energiazeroel
    real*8                          :: dt
    integer                         :: neqn
    real*8,             allocatable :: frequency_matrix(:,:)
    real*8,             allocatable :: h_mtx(:,:), RWMatrix(:,:)
    type(quantum_site), target, allocatable :: site(:,:)
    type(obj_pointer)         , allocatable :: site_point(:)
    type(BasisBuild)          , allocatable :: basis(:,:)

end module parameters_m
