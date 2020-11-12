module parameters_m
use types_m
use constants_m 
!modulo que tem todos os parametros que usamos apenas no main.f90
!TODOS OS TEMPOS ESTAO EM UNIDADES DE PS E AS FREQUENCIAS EM UNIDADE DE THZ

    real*8,  parameter :: tmin = 0.d0                  !tempo inicial
    real*8,  parameter :: tmax = 10.d0                !tempo final
    integer, parameter :: nsites = 2 
    integer, parameter :: nr = 1                       !numero de linhas para os sitios
    integer, parameter :: nc = 2                       !numero de colunas para os sitios
    integer, parameter :: BasisChoice = 2              !{(n=0, l=0)} : 1-> (n=0, l=0) 2-> (n=0, l=-1), (n=0, l=1) 3-> 1 e 2 
    real*8,  parameter :: distance_molecule = 0.9d-9   !distancia entre os sitios
    integer, parameter :: nm_divisoes = 25000          !numero de passos utilizado no programa
    integer, parameter :: initState = 2                !estado inicial do eletron
    real*8,  parameter :: siteMass = 1.9936d-26        !massa dos sitios
    real*8,  parameter :: raioZero = 0.35d-9           !raio inicial dos sitios
    real*8,  parameter :: velZero  = 0.d0              !2.d3
    real*8,  parameter :: effectiveMassMult = 5.d0     !fator de multiplicacao da massa efetiva do eletron
    real*8,  parameter :: siteCoupling = 0.6d0         != tc

    integer, parameter :: BRD_on = 1                   !1 para berendsen desligado, 0 para ligado 

    real*8, parameter :: freqMode = 8.154356d12        !frequencia do modo vibracional em Hz 
    real*8, parameter :: bathTemp = 150.d0             !temperatura fixa em kelvin
    real*8, parameter :: bathCoup = 1.d-14             !constante de acoplamento Berendsen em segundos
    
    integer, parameter :: ns_el = BasisChoice          !numero de estados para o eletron 
    integer, parameter :: d_el = (nsites*ns_el)        !tamanho da matriz rho, hamiltoniana
    integer, parameter :: neqn_el = d_el * d_el        !numero de equacoes diferenciais acopladas

    real*8, allocatable             :: Qn_erg(:)
    integer, allocatable            :: Qn(:)
   
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
