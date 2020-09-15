module parameters_m
use types_m
use constants_m 
!modulo que tem todos os parametros que usamos apenas no main.f90
!TODOS OS TEMPOS ESTAO EM UNIDADES DE PS E AS FREQUENCIAS EM UNIDADE DE THZ

real*8, parameter :: tmin = 0.d0  !tempo inicial que utilizamos na integral para calcular o tensor de Redfield
real*8, parameter :: tmax = 1.2d0  !tempo maximo que utilizamos na evolucao
real*8, parameter :: temp = 347.d0        !temperatura em kelvin
real*8, parameter :: coupling = 1.d0 !teste para modificar o acoplamento eletronico 
integer, parameter :: nsites = 2  !numero de sítios que entram na simulacao
integer, parameter :: nsites_rec = 0 !numero de sítios de recombinacao
integer, parameter :: nqdots = nsites + nsites_rec !numero total 
integer, parameter :: nmstates_el = 3 !numero de estados para o eletron
integer, parameter :: nmstates_hl = 3 !numero de estados para o buraco 
integer, parameter :: dim_el = (nsites*nmstates_el) + nsites_rec
integer, parameter :: dim_hl = (nsites*nmstates_hl) + nsites_rec
integer, parameter :: nm_rows = 1 !numero de linhas para os sitios
integer, parameter :: nm_columns = 2 !numero de colunas para os sitios
real*8, parameter :: v_dephasing_el = 0.d0  !0.15d0 !0.15d0 !v de dephasing para o eletron
real*8, parameter  :: v_dephasing_hl = 0.d0 !0.15d0 !0.15d0 !v de dephasing para o buraco 
!real*8, parameter :: taxa_rec = 0.5d0 
real*8, parameter :: distance_molecule = 0.25d-9  !distancia entre os sitios
real*8, parameter :: distance_layer = 0.25d-9 !0.35d-9 !distancia na heterojuncao
real*8, parameter :: distance_eletrodes_sites = 0.25d-9 !0.35d-9 !!0.25 !distancia dos eletrodos
real*8, parameter :: distance_eletrodes = 0.25d-9 !0.35d-9 !distancia dos eletrodos
integer, parameter :: neqn_el = dim_el * dim_el !numero de equacoes diferenciais acopladas que queremos resolver ########### resolve equacao nao homogenea tbm (testado) ###########
integer, parameter :: neqn_hl = dim_hl * dim_hl
real*8, parameter :: lim_values = 0.d0 !1.d-2
integer, parameter :: nm_divisoes = 3000 !numero de passos utilizado no programa
!integer, parameter :: half_ndim = (nsites *(nsites-1)) / 2 !numero de eq. pra resolver na matriz - modo 1
!integer, parameter :: half_ndim = 2*nsites-1 !numero de eq. pra resolver na matriz - modo 2 (relaxacao separada por sitio)
integer, parameter :: half_ndim = nsites + 1 !nsites !numero de eq. pra resolver na matriz - modo 2 (relaxacao separada por sitio)
real*8, parameter :: Dterm = 0.d0 !1.d-16
integer, parameter :: initState = 1



real*8, parameter :: siteMass = 1.d-27
real*8, parameter :: raioZero = 0.5d-9
real*8, parameter :: effectiveMassMult = 5.d0
real*8, parameter :: siteCoupling = 0.6d0
real*8, parameter :: elhlCouplingCte = 0.1d0


real*8, parameter :: multFactorDer =  1.d0 !fator que estou multiplicando os elementos nao diagonais da matriz derivada
real*8, parameter :: multFactorOvlp = 1.d0 !fator que estou multiplicando os elementos nao diagonais da matriz de overlap
real*8, parameter :: eletricCte = 1.d0     !fator que estou multiplicando a força eletrica 
real*8, parameter :: cteForce = 1.d0       !fator que estou multiplicando a força da mola 


!nm de sitios para cada material. nqdots é a soma de todos eles
integer, parameter :: nm_anodes_sites  = 2
integer, parameter :: nm_donnors_sites = 1 
integer, parameter :: nm_acceptors_sites = 1 
integer, parameter :: nm_cathodes_sites  = 1
!---------------------------------------------------------------
!nqdots deve ser a soma de nm_anodes + nm_donnors +  nm_acceptors + nm_cathode


!nm de colunas para cada material
integer, parameter :: nm_anode_columns  = 2
integer, parameter :: nm_donnor_columns = 1 !5  !12
integer, parameter :: nm_acceptor_columns = 1   !5 !8
integer, parameter :: nm_cathode_columns = 1
!--------------------------

!parameterização de cada coluna para cada material
integer, parameter :: anode_layer_b = 1 !anode layer begin
integer, parameter :: anode_layer_e = 1 !anode layer end
integer, parameter :: donnor_layer_b  = 3
integer, parameter :: donnor_layer_e  = 3 
integer, parameter :: acceptor_layer_b = 0  
integer, parameter :: acceptor_layer_e = 0  
integer, parameter :: cathode_layer_b = 0 
integer, parameter :: cathode_layer_e = 0 


real*8, parameter :: me = effectiveMassMult * electronMass !massa efetiva total



!V_DEPHASING ERA O ANTIGO DR E NAO TEM MAIS DIMENSAO POIS MULTIPLICO POR HBAR AS CONTAS JA NA DEFINIÇÃO DO H_INT
!DESSA FORMA  A FUNCAO DE CORRELACAO TEM UNIDADE DE FREQUENCIA
!COUPLING É DEFINIDO DE FORMA QUE PARA T=300K A TAXA SEJA A DO PAPER = 5 PS^-1 = 26.52 CM
!ISSO FORNECCE 26.52 = 2 * COUPLING * Kb * T => COUPLING = 26.52/Kb*T = 26.52/416.988 = 0.0635





integer, parameter :: xypoints = 500
integer, parameter :: max_hmt = 6
real*8, parameter :: xmin = -5.d0
real*8, parameter :: xmax = 20.d0
real*8, parameter :: ymin = -9.d0
real*8, parameter :: ymax = 9.d0


real*8, parameter :: EF = 0.d0 



real*8 :: trace_rho_el, trace_rho_hl 
integer :: dims, neqn, nvezes
real*8 :: trace_rho
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_ham_in
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: frequency_matrix
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: hamiltoniano, RWMatrix
real*8 :: delta_t
real*8, allocatable :: spectral_density(:)

integer, allocatable :: basis(:, :)
real*8, allocatable :: OPS_site(:, :, :)
real*8, allocatable :: system_operators(:, :, :), system_operators_el(:, :, :), system_operators_hl(:, :, :)
real*8, allocatable, dimension(:, :) :: ovlpm_el, ovlpm_hl 
real*8 :: energiazeroel, energiazerohl
type(quantum_site), allocatable ::  sites_array(:, :)


!integer, parameter :: ndim = ((dimensao*(dimensao + 1)) / 2) !numero de eq. pra resolver na matriz rho_real


end module parameters_m
