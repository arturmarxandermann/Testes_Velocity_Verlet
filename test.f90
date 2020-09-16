module testes
use f95_precision
use blas95
use lapack95
use types_m
use parameters_m
use constants_m
use functions_m
use system_hamiltonian_m
use rdftensor_m
use ordinary_equation_solution
use omp_lib 
use verlet_m 

private

public :: test, f0, call_ODE_solver


contains


!se as matrizes são definidas no parameters.f90, não preciso criar elas como intent(in)
!pois todas as subrotinas utilizam parameters.f90


subroutine call_ODE_solver(nm_divisoes)
implicit none
integer, intent(in) :: nm_divisoes
integer :: pl
real*8 :: ti, tf
complex*16, DIMENSION(:, :), ALLOCATABLE :: rho_el_sites_in !matriz densidade inicial
REAL*8, ALLOCATABLE, DIMENSION(:, :):: pop_el_in, pop_el !populacoes eletron
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_el !matriz densidade calculada em t + delta t

real*8, dimension(:, :), allocatable :: BWEletricForce !aloco na subrotina eletric_force
real*8, dimension(:, :), allocatable :: BWSpringForce !aloco na subrotina spring_force
real*8 :: EnergiaEl !energia do eletron
real*8, dimension(:, :), allocatable :: hamiltoniana_el_in !hamiltoniano inicial
real*8, dimension(:, :), allocatable :: hamiltoniana_el  !hamitoniano calculado em t + delta t


real*8, dimension(:, :), allocatable :: FWRadius, FWRadialVel, FWEletricForce, FWSpringForce !Forward variables
real*8, dimension(nm_rows, nm_columns) :: BWRadius, BWRadialVel !Backward variables
real*8, dimension(:, :), allocatable :: rho_el_real !parte real da matriz densidade
real*8 :: SpringEnergy, KinectEnergy !energia da mola e energia cinetica


ALLOCATE(rho_el_sites_in(dim_el, dim_el), source = (0.d0, 0.d0))
ALLOCATE(pop_el_in(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(rho_el_real(dim_el, dim_el ), source = 0.d0 ) 




call define_sitios(sites_array) 

!============= CONDICOES INICIAIS NA BASE DOS SITIOS ==============
  rho_el_sites_in(initState, initState) = 1.d0 + 0.d0 * zi  !matriz densidade
  call rho_matrix_to_pop(rho_el_sites_in, pop_el_in) !populacao de cada sitio
  BWRadius = sites_array(:, :)%radius !raio inicial
  BWRadialVel = sites_array(:, :)%radial_vel !velocidade inicial
  rho_el_real = real(rho_el_sites_in) !parte real da matriz densidade
  call build_hamiltonian(hamiltoniana_el_in) !hamiltoniano inicial
  call eletric_force(pl, rho_el_real, hamiltoniana_el_in, BWEletricForce) !Backward EletricForce
  call spring_force(BWSpringForce) !Backward SpringForce
!================================================================

!======= SUBROTINA QUE ABRE OS ARQUIVOS DE ESCRITA =============
call open_write_files
!==============================================================

ti = 0.d0


do pl = 1, nm_divisoes !NÃO POSSO USAR i NESSSE LOOP PQ ELE JÁ É UTILIZADO EM TODOS OS OUTROS LOOPS NA SUBROTINA TEST01..
   if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
      print*, "==== PASSO",pl,"===="
    endif
  
    !============= DEFINE O TEMPO FINAL E O DELTA_t ===========
    tf = float(pl)*(tmax)/float(nm_divisoes)
    delta_t = tf - ti 
    !==========================================================



    !=============== CALCULO DO HAMILTONIANO  ================
    call build_hamiltonian(hamiltoniana_el)
    !==========================================================


    !============== EVOLUCAO PARA O ELETRON DE ti ATÉ tf =======
    call test(pl, ti, tf, hamiltoniana_el, rho_el_sites_in, pop_el_in, rho_el, EnergiaEl, pop_el)
    !===========================================================
  
    !============= ATUALIZO AS CONDICOES INICIAIS QUANTICAS ==============
    rho_el_sites_in = rho_el !atualizo cond. inicial para o eletron
    pop_el_in = pop_el
    rho_el_real = real(rho_el_sites_in) 
    !===========================================================

    
    !============ CALCULO A EVOLUCAO CLASSICA COM O ALGORITMO VELOCITY VERLET ===========================

    call velocity_verlet(pl, BWRadius, BWRadialVel, BWEletricForce, BWSpringForce, delta_t, rho_el_real, &
                      hamiltoniana_el, FWRadius, FWRadialVel, FWEletricForce, FWSpringForce, SpringEnergy, KinectEnergy) 
              !calculamos a energia cinetica e da mola no algoritmo de verlet em um tempo ti para ser igual ao resto do programa
    !====================================================================================================

    !============ ESCREVE AS ENERGIAS CLASSICAS E QUANTICAS ============================================
    write(84, "(60F20.8)") ti, KinectEnergy
    write(85, "(60F20.8)") ti, SpringEnergy
    write(86, "(60F20.8)") ti, energiazeroel
    write(87, "(60F20.8)") ti, EnergiaEl
    write(103, "(60F20.8)") ti, energiazeroel, EnergiaEl + SpringEnergy + KinectEnergy
    !===================================================================================================
  

    !=========== ATUALIZA A PARTE CLASSICA =============================================================
    BWRadius(:, :) = FWRadius(:, :) !aloco em velocity_verlet
    BWRadialVel(:, :) = FWRadialVel(:, :) !aloco em velocity_verlet
    BWEletricForce(:, :) = FWEletricForce(:, :) !aloco em eletric_force
    BWSpringForce(:, :) = FWSpringForce(:, :) !aloco em spring_force
    deallocate(FWEletricForce, FWRadialVel, FWRadius, FWSpringForce) 
    !===================================================================================================
  

    !========== ESCREVE OS RAIOS E AS FORÇAS ===========================================================
    write(100, "(60F20.8)"), tf, sites_array(1, 1)%radius*1.d9, sites_array(1, 2)%radius*1.d9 
    write(101, "(60F20.8)"), tf, BWEletricForce(1, 1)*1.d11, BWSpringForce(1, 1)*1.d11, BWEletricForce(1, 2)*1.d11 &
          , BWSpringForce(1, 2)*1.d11
    !===================================================================================================
  
    !========= ATUALIZO O TEMPO ===========
    ti = tf
    !=====================================
enddo

!========= FECHO OS ARQUIVOS DE ESCRITA =====
call close_write_files
!============================================


return

DEALLOCATE(pop_el_in, rho_el_sites_in)

end subroutine call_ODE_solver



subroutine test(step, ti, tf, hamiltoniana, rho_el_sites_in, pop_el_in, rho_sites, ParticleEnergy, pop)
implicit none

integer, intent(in) :: step 
REAL*8, INTENT(IN) :: ti !tempo inicial que vamos utilizar no ODE e que não  vai ser atualizado uma vez chamada a rotina
REAL*8, INTENT(IN) :: tf ! tempo final que nao vai ser atualizado uma vez chamado a rotina
real*8, dimension(dim_el, dim_el),  intent(in) :: hamiltoniana
COMPLEX*16, INTENT(IN) :: rho_el_sites_in(dim_el, dim_el)
REAL*8, INTENT(IN) :: pop_el_in(nm_rows, nm_columns)
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: rho_sites
real*8, intent(out) :: ParticleEnergy
REAL*8, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: pop


real*8, ALLOCATABLE, DIMENSION(:) :: energias, y, yp, work
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: phi, phi_transpose
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_sites_in, matriz_temporaria
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: pop_in
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_ham
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: rhomtx
INTEGER :: number_file_site, number_file_ham

real ( kind = 8 ) abserr
integer ( kind = 4 ) iflag
integer :: flag
integer ( kind = 4 ) iwork(5)
real ( kind = 8 ) relerr
real ( kind = 8 ) t !esse t é a variável independente do ODE
real ( kind = 8 ) tout 
real*8 :: energ(dim_el) 



 ALLOCATE(rho_sites_in(dim_el, dim_el), source = (0.d0, 0.d0))
 ALLOCATE(pop_in(nm_rows, nm_columns), source = 0.d0 )
 rho_sites_in = rho_el_sites_in
 pop_in = pop_el_in
 number_file_ham = 14
 number_file_site = 15

  

ALLOCATE(y(neqn_el), source = 0.d0)
ALLOCATE(yp(neqn_el), source = 0.d0)
ALLOCATE(matriz_temporaria(dim_el, dim_el), source = (0.d0, 0.d0) )

ALLOCATE(rho_sites(dim_el, dim_el), source = (0.d0, 0.d0) )

ALLOCATE(pop(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(work(100+21*neqn_el))

ALLOCATE(rho_ham_in(dim_el, dim_el), source = (0.d0, 0.d0) ) !populaca inicial na base da hamiltoniana
ALLOCATE(rho_ham(dim_el, dim_el), source = (0.d0, 0.d0) )


ALLOCATE(rhomtx(dim_el, dim_el), source =  0.d0 ) 






abserr = 1.d-9
relerr = 1.d-9
flag = +1 !FLAG -1 NAO FUNCIONA
iflag = 1 
t = ti !construo o tempo inicial de uma forma que o ODE nao atualiza!
tout = tf !TOUT É O TEMPO DE SAÍDA, OU SEJA, DELTA T = tout(2) - tout(1)



!==== CALCULA OS AUTOESTADOS E AUTOVETORES =======
call calculate_eigenvectors(step, hamiltoniana, energias, phi, phi_transpose, frequency_matrix) 
!===============================================================================================


!==== PRINTA A MATRIZ DENSIDADE, AUTOVETORES E HAMILTONIANO =====
call print_matrices(step, rho_sites_in, phi, hamiltoniana) 
!===============================================================================================


!==== DEFINE A ENERGIA ZERO ===========
if (t == 0.d0 ) then
  energiazeroel = hamiltoniana(initState, initState)
endif
!=====================================


!=========== CONDICOES INICIAIS NA BASE EXC ======================
call rhosite_TO_rhoham(phi, phi_transpose, rho_sites_in, rho_ham_in)
!=================================================================



!========== CALCULO DA ENERGIA DA PARTICULA ====================
call particle_energy(hamiltoniana, rho_sites_in, ParticleEnergy)
!==============================================================


!======== TRANSFORMO AS CONDICOES INICIAIS PARA OS VETORES LINHAS DO ODE ========
call monta_y(rho_ham_in, y) !condição inicial para resolver as ODE'S
!================================================================================


!====== CRIO A MATRIZ RW QUE É UTILIZADA PARA CALCULAR O ODE =================
call createRwMatrix(frequency_matrix, RWMatrix) 
!=============================================================================



!y esta na forma (para d=2):  y(1) = rho(1, 1)  y(3) = real(rho(2, 1))
!                             y(2) = rho(2, 2)  y(4) = aimag(rho(2, 1))
!----------------------------------------------------------------------------------
  

!====== ESCREVE AS POPULACOES =====================
call printa_resultado(number_file_ham, t, rho_ham_in)
call printa_resultado(number_file_site, t, rho_sites_in)
!==================================================



!======== SUBROTINA QUE CALCULA AS EQUACOES DIFERENCIAIS DE T ATÉ TOUT ========
call ode ( f0, neqn_el, y, t, tout, relerr, abserr, iflag, work, iwork )

if ( iflag /= 2 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST - Fatal error!'
  write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
  stop
end if
!==============================================================================



!========== TRANSFORMO O RESULTADO DO VETOR LINHA Y DO ODE PARA AS MATRIZES RHO ===========
call monta_rho(y, rho_ham)
!==========================================================================================


!============== ESCREVO O OPERADOR DENSIDADE NA BASE DO SITIO =============================
call rhoham_TO_rhosite(phi, phi_transpose, rho_ham, rho_sites)
!==========================================================================================



!=========== CALCULO AS POPULACOES DE CARGA EM CADA SITIO ================================
call rho_matrix_to_pop(rho_sites, pop)
!=========================================================================================


    

DEALLOCATE(phi, phi_transpose, rhomtx, rho_sites_in, rho_ham_in, matriz_temporaria, pop_in, RWMatrix)
13 format(3es14.3E3)
return
end subroutine test



subroutine f0 ( t, y, yp )
!não definimos as condicoes iniciais aqui, apenas a funcao que vamos calcular a derivada. A condição inicial vem antes do call
implicit none
real*8 ::  t !tempo que vai entrar na equação diferencial
real*8 ::  y(neqn_el)  !funcao que estamos querendo, neste caso x
real*8 ::  yp(neqn_el) !derivada que estamos querendo calcular, por exemplo dx/dt = -x => sol: x = exp(-t)


call gemv(RWMatrix, y, yp)  !gemv calcula o produto da matriz RWMatrix com
                            !o vetor y resultando em yp => muito mais rapido

return
13 format(3es14.3E3)
end subroutine f0



end module testes
