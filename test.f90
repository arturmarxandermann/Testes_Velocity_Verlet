module testes
use f95_precision
use blas95
use lapack95
use types_m
use parameters_m
use constants_m
use functions_m
use system_hamiltonian_m
use system_operators_m
use rdftensor_m
use ordinary_equation_solution
use gnugraphs_m
!use rkf45_m
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
complex*16, DIMENSION(:, :), ALLOCATABLE :: rho_el_sites_in, rho_hl_sites_in
REAL*8, ALLOCATABLE, DIMENSION(:, :):: pop_el_in, pop_el, pop_hl_in, pop_hl
real*8, DIMENSION(nm_divisoes) :: cathode_pop, anode_pop
real*8, DIMENSION(nm_divisoes) :: time_grid
real*8, DIMENSION(nm_divisoes) :: eletron_current, hole_current
integer :: counter, image_number, ver, k
real*8, parameter :: ponto_simulacao = tmax/1.d0 !a partir deste ponto pego menos pontos na simulacao
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_el, rho_hl
real*8 :: tb, ta, tmax1 
real*8 :: taxarec_el, taxarec_hl 
complex*16, dimension(:, :), allocatable :: rho_rec_site_out_el, rho_rec_site_out_hl

real*8, dimension(:, :), allocatable :: BWEletricForce !aloco na subrotina eletric_force
real*8, dimension(:, :), allocatable :: BWSpringForce
real*8 :: EnergiaEl, EnergiaHl
real*8, dimension(:, :), allocatable :: hamiltoniana_el_in, hamiltoniana_hl_in
real*8, dimension(:, :), allocatable :: hamiltoniana_el, hamiltoniana_hl


real*8, dimension(:, :), allocatable :: FWRadius, FWRadialVel, FWEletricForce, FWSpringForce
real*8, dimension(nm_rows, nm_columns) :: BWRadius, BWRadialVel
real*8, dimension(:, :), allocatable :: rho_el_real, rho_hl_real
!real*8, dimension(:, :), allocatable :: SpringEnergy
real*8 :: SpringEnergy, KinectEnergy
real*8, dimension(:, :), allocatable :: DerOvlp


ALLOCATE(rho_el_sites_in(dim_el, dim_el), source = (0.d0, 0.d0))
ALLOCATE(rho_hl_sites_in(dim_hl, dim_hl), source = (0.d0, 0.d0))
ALLOCATE(pop_el_in(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(rho_el_real(dim_el, dim_el ), source = 0.d0 ) 
ALLOCATE(rho_hl_real(dim_hl, dim_hl ), source = 0.d0 ) 
ALLOCATE(pop_hl_in(nm_rows, nm_columns), source = 0.d0)

allocate(hamiltoniana_hl_in(dim_hl, dim_hl), source = 0.d0) 
allocate(hamiltoniana_hl(dim_hl, dim_hl), source = 0.d0) 


!allocate( rho_hl_sites_in(dim_hl, dim_hl), source = (0.d0, 0.d0) ) 

call define_sitios(sites_array) 

!condicoes iniciais na base dos sitios 
rho_el_sites_in(initState, initState) = 1.d0 + 0.d0 * zi  !condicao inicial para o eletron
rho_hl_sites_in(initState, initState) = 0.d0 + 0.d0 * zi !condicao inical para o buraco


call rho_matrix_to_pop(nmstates_el, dim_el, rho_el_sites_in, pop_el_in)
call rho_matrix_to_pop(nmstates_hl, dim_hl, rho_hl_sites_in, pop_hl_in)


BWRadius = sites_array(:, :)%radius
BWRadialVel = sites_array(:, :)%radial_vel

rho_el_real = real(rho_el_sites_in)
rho_hl_real = real(rho_hl_sites_in)

  call build_hamiltonian(1, pop_el_in, pop_hl_in, ti, hamiltoniana_el_in)
  call build_hamiltonian(2, pop_el_in, pop_hl_in, ti, hamiltoniana_hl_in)


  call eletric_force(1, dim_el, rho_el_real, rho_hl_real, hamiltoniana_el_in, hamiltoniana_hl_in, BWEletricForce) !Backward EletricForce
  call spring_force(1, dim_el, BWSpringForce) !, SpringEnergy) !Backward SpringForce




!---------------------------------------------------
open(20, file = 'plot.gnu', status = 'replace') !arquivos para os filmes
open(17, file = 'fort.17', status = 'replace')  !buraco base local
open(16, file = 'fort.16', status = 'replace')  !buraco base desloc.
open(15, file = 'fort.15', status = 'replace')  !eletron base local
open(13, file = 'fort.13', status = 'replace')  !eletron base desloc.
open(14, file = 'fort.14', status = 'replace')  !eletron base desloc.
open(10, file = 'fort.10', status = 'replace')  !evol. el com (RHO * S + S * RHO) / 2
open(11, file = 'fort.11', status = 'replace')  !evol. hl com (RHO * S + S * RHO) /2
open(90, file = 'fort.90', status = 'replace') 
open(100, file = 'radius.txt', status = 'replace') 
open(101, file = 'forces.txt', status = 'replace') 
open(103, file = 'energiatotal.txt', status = 'replace') 
open(79, file = 'fort.79', status = 'replace') 
open(84, file = 'fort.84', status = 'replace') 
open(70, file = 'fort.70', status = 'replace')
open(86, file = 'fort.86', status = 'replace')
open(87, file = 'fort.87', status = 'replace')


counter = 1
image_number = 1
ti = 0.d0



!call GnuFiles(sites_array, pop_el_in, pop_hl_in, image_number, ti)

do pl = 1, nm_divisoes !NÃO POSSO USAR i NESSSE LOOP PQ ELE JÁ É UTILIZADO EM TODOS OS OUTROS LOOPS NA SUBROTINA TEST01..
  ver = int(mod(pl, 500))
   if (ver == 0 ) then
      print*, "==== PASSO",pl,"===="
    endif
  ta = omp_get_wtime()
  
  time_grid(pl) = ti  !transforma o tempo para ps com objetivo de calcular a corrente
  
  tf = float(pl)*(tmax)/float(nm_divisoes)

  delta_t = tf - ti 



  !========== PARTE REC
  !!trace_rho_el = 0.d0 
  !!trace_rho_hl = 0.d0 
  
  !!do k  = 1, dim_el
  !!  trace_rho_el = trace_rho_el + real(rho_el_sites_in(k, k))
  !!  trace_rho_hl = trace_rho_hl + real(rho_hl_sites_in(k, k))  
  !!enddo


  !!taxarec_el = 0.d0
  !!taxarec_hl = 0.d0

  !!taxarec_el = trace_rho_el * taxa_rec
  !!taxarec_hl = trace_rho_hl * taxa_rec 


  !!do i = 1, dims
  !!rho_el_sites_in(i, i) =  real(rho_el_sites_in(i, i)) * exp(-taxarec_el*real(rho_hl_sites_in(i, i)) * (tf - ti))
  !!rho_hl_sites_in(i, i) =  real(rho_hl_sites_in(i, i)) * exp(-taxarec_hl*real(rho_el_sites_in(i, i)) * (tf - ti))
  !!enddo

  !=============
  

  ! =============== calculando a variação nos raios dos sítios ========


  print*, "=========== PASSO ===========", pl
  call build_hamiltonian(1, pop_el_in, pop_hl_in, ti, hamiltoniana_el)
  call build_hamiltonian(2, pop_el_in, pop_hl_in, ti, hamiltoniana_hl)

  hamiltoniana_hl = 0.d0
  pop_hl_in = 0.d0
  !======================================================




  !EVOLUCAO PARA O ELETRON DE ti ATÉ tf
  call test(1, ti, tf, hamiltoniana_el, rho_el_sites_in, rho_hl_sites_in, pop_el_in, pop_hl_in, rho_el, EnergiaEl, pop_el)

  !EVOLUCAO PARA O BURACO DE ti ATÉ tf
  !call test(2, ti, tf, hamiltoniana_hl, rho_el_sites_in, rho_hl_sites_in, pop_el_in, pop_hl_in, rho_hl, EnergiaHl, pop_hl)
  
  EnergiaHl = 0.d0

  rho_el_sites_in = rho_el !atualizo cond. inicial para o eletron
  pop_el_in = pop_el
  !rho_hl_sites_in = rho_hl !atualizo cond. inicial para o buraco 
  !pop_hl_in = pop_hl
  tb = omp_get_wtime()
  !!print*, "tempo para o passo", pl, "é", tb - ta
  
  rho_el_real = real(rho_el_sites_in) 
  !rho_hl_real = real(rho_hl_sites_in) 
  

  hamiltoniana_hl = 0.d0 

   call velocity_verlet(BWRadius, BWRadialVel, BWEletricForce, BWSpringForce, ti, delta_t, rho_el_real, &
                      rho_hl_real, hamiltoniana_el, hamiltoniana_hl, FWRadius, FWRadialVel, FWEletricForce, &
                      FWSpringForce, SpringEnergy, KinectEnergy) 
              !calculamos a energia cinetica e da mola no algoritmo de verlet em um tempo ti para ser igual ao resto do programa
  
  write(85, "(60F20.8)") ti, SpringEnergy, KinectEnergy
  write(86, "(60F20.8)") ti, energiazeroel + energiazerohl
  write(87, "(60F20.8)") ti, EnergiaEl
  write(103, "(60F20.8)") ti, energiazeroel + energiazerohl, EnergiaHl + EnergiaEl + SpringEnergy + KinectEnergy

  
   do j = 1, nm_columns
    do i = 1, nm_rows
    BWRadius(i, j) = FWRadius(i, j) !aloco em velocity_verlet
    BWRadialVel(i, j) = FWRadialVel(i, j) !aloco em velocity_verlet
    BWEletricForce(i, j) = FWEletricForce(i, j) !aloco em eletric_force
    BWSpringForce(i, j) = FWSpringForce(i, j) !aloco em spring_force
    enddo
  enddo
  deallocate(FWEletricForce, FWRadialVel, FWRadius, FWSpringForce) !, SpringEnergy) 
  
  

 
  write(100, "(60F20.8)"), tf, sites_array(1, 1)%radius*1.d9, sites_array(1, 2)%radius*1.d9 
  write(101, "(60F20.8)"), tf, BWEletricForce(1, 1)*1.d11, BWSpringForce(1, 1)*1.d11, BWEletricForce(1, 2)*1.d11 &
          , BWSpringForce(1, 2)*1.d11

  image_number = image_number + 1
!  call GnuFiles(sites_array, pop_el, pop_hl, image_number, tf)

  ti = tf
enddo

close(13)
close(14)
close(15)
close(10)
close(11)
close(16)
close(17)
close(20)
close(90)
close(101)
close(100)
close(79)
close(84)
close(70) 
close(86) 
close(87) 
close(103) 

return

DEALLOCATE(pop_el_in, pop_hl_in, rho_el_sites_in, rho_hl_sites_in, rho_rec_site_out_el, rho_rec_site_out_hl)

end subroutine call_ODE_solver



subroutine test(el_or_hl, ti, tf, hamiltoniana, rho_el_sites_in, rho_hl_sites_in, pop_el_in, pop_hl_in, rho_sites, ParticleEnergy, pop)
implicit none

integer, INTENT(IN) :: el_or_hl
REAL*8, INTENT(IN) :: ti !tempo inicial que vamos utilizar no ODE e que não  vai ser atualizado uma vez chamada a rotina
REAL*8, INTENT(IN) :: tf ! tempo final que nao vai ser atualizado uma vez chamado a rotina
real*8, dimension(dim_el, dim_el),  intent(in) :: hamiltoniana
COMPLEX*16, INTENT(IN) :: rho_el_sites_in(dim_el, dim_el), rho_hl_sites_in(dim_hl, dim_hl)
REAL*8, INTENT(IN) :: pop_el_in(nm_rows, nm_columns), pop_hl_in(nm_rows, nm_columns)
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: rho_sites
real*8, intent(out) :: ParticleEnergy
REAL*8, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: pop

REAL*8, ALLOCATABLE, DIMENSION(:, :, :, :) :: RDtensor

real*8, ALLOCATABLE, DIMENSION(:) :: energias, y, yp, work
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: phi, phi_transpose
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_sites_in, matriz_temporaria
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_sites_in_squared, energymtc
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: pop_in
COMPLEX*16, ALLOCATABLE, DIMENSION(:, :) :: rho_ham, rho_rec 
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: multmtx
real*8, allocatable, dimension(:, :) :: ovlpm  !, ovlpmcopy, INovlpm
INTEGER :: number_file_site, number_file_ham, nstates
real*8 :: t1, t2, t3, t4, t5, t6
real*8, dimension(:, :), allocatable :: results!, IS_ovlpm, ovlpm_copy

real ( kind = 8 ) abserr
integer ( kind = 4 ) iflag
integer :: flag, number_file_pop_ovl
integer ( kind = 4 ) iwork(5)
real ( kind = 8 ) relerr
real ( kind = 8 ) t !esse t é a variável independente do ODE
real ( kind = 8 ) tout 
real*8 :: ta, tb
real*8 :: purity, energ(dim_el) 


if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if




if (el_or_hl == 1 ) then
  !print*, "calculando Hs eletron" 

  dims = dim_el
  neqn = neqn_el
  nstates = nmstates_el
  ALLOCATE(rho_sites_in(dims, dims), source = (0.d0, 0.d0))
  ALLOCATE(pop_in(nm_rows, nm_columns), source = 0.d0 )
  rho_sites_in = rho_el_sites_in
  pop_in = pop_el_in
  ALLOCATE(rho_rec(dims, dims), source = (0.d0, 0.d0) )
  rho_rec = rho_hl_sites_in
  number_file_pop_ovl = 10
  number_file_ham = 14
  number_file_site = 15

  
end if

if (el_or_hl == 2 ) then
  !print*, "calculando Hs buraco"
  dims = dim_hl
  neqn = neqn_hl
  nstates = nmstates_hl
  ALLOCATE(rho_sites_in(dims, dims), source = (0.d0, 0.d0))
  ALLOCATE(pop_in(nm_rows, nm_columns), source = 0.d0 )
  rho_sites_in = rho_hl_sites_in
  pop_in = pop_hl_in
  ALLOCATE(rho_rec(dims, dims), source = (0.d0, 0.d0) )
  rho_rec = rho_el_sites_in
  number_file_pop_ovl = 11
  number_file_ham = 16
  number_file_site = 17


end if


ALLOCATE(y(neqn), source = 0.d0)
ALLOCATE(yp(neqn), source = 0.d0)
ALLOCATE(matriz_temporaria(dims, dims), source = (0.d0, 0.d0) )

ALLOCATE(rho_sites(dims, dims), source = (0.d0, 0.d0) )

ALLOCATE(pop(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(work(100+21*neqn))

ALLOCATE(rho_ham_in(dims, dims), source = (0.d0, 0.d0) ) !populaca inicial na base da hamiltoniana
ALLOCATE(rho_ham(dims, dims), source = (0.d0, 0.d0) )
ALLOCATE(energymtc(dims, dims), source = (0.d0, 0.d0) ) 
ALLOCATE(rho_sites_in_squared(dims, dims), source = (0.d0, 0.d0) ) 

ALLOCATE(ovlpm(dims, dims), source = 0.d0 ) 
!ALLOCATE(INovlpm(dims, dims), source = 0.d0 ) 
!ALLOCATE(ovlpmcopy(dims, dims), source = 0.d0 ) 



ALLOCATE(multmtx(dims, dims), source =  0.d0 ) 
ALLOCATE(results(dims, dims), source =  0.d0 ) 






abserr = 1.d-9
relerr = 1.d-9
flag = +1 !FLAG -1 NAO FUNCIONA
iflag = 1 
t = ti !construo o tempo inicial de uma forma que o ODE nao atualiza!
tout = tf !TOUT É O TEMPO DE SAÍDA, OU SEJA, DELTA T = tout(2) - tout(1)



!===== CALCULA A MATRIZ DE OVERLAP ENTRE OS SITIOS ======
call build_overlap_matrix(nstates, dims, ovlpm)


!==== CALCULA OS AUTOESTADOS E AUTOVETORES EM UMA BASE NÃO ORTOGONAL (depende de IS_ovlpm)  =======
call calculate_eigenvectors(dims, hamiltoniana, energias, phi, phi_transpose, frequency_matrix) !, Ztrf)
!===============================================================================================

if (el_or_hl == 1) then
multmtx = 0.d0
multmtx = real(rho_sites_in)
print*, "PARTE REAL COMECO DO PROGRAMA"
call print_mat2(multmtx, dims, dims) 
multmtx = 0.d0 
multmtx = aimag(rho_sites_in)
print*, "PARTE IMAGINARIA COMECO DO PROGRAMA"
call print_mat2(multmtx, dims, dims) 
print*, "autovetores"
call print_mat2(phi, dims, dims)
print*, "hamiltoniano"
call print_mat2(hamiltoniana, dims, dims) 
endif 



if (el_or_hl == 1 .AND. t == 0.d0 ) then
  energiazeroel = hamiltoniana(initState, initState)
endif

if (el_or_hl == 2 .AND. t == 0.d0 ) then
  energiazerohl = hamiltoniana(initState, initState)
endif


!=== CALCULA O TENSOR DE REDFIELD PARA TODOS OS ESTADOS ====
t1 = omp_get_wtime()
call create_redfield_tensor(el_or_hl, dims, phi, frequency_matrix, ovlpm, RDtensor)
t2 = omp_get_wtime()


!if (ti .EQ. 0.d0) then
!    if (el_or_hl == 1 ) then
!    call printaletters2('el_ham', 1, dims, dims, dims, hamiltoniana)
!    call printaletters2('el_phi', 2, dims, dims, dims, phi)
!    call printaletters2('el_ovm', 3, dims, dims, dims, ovlpm)
!!    call printaletters2('el_Ztr', 4, dims, dims, (nsites/2)*nstates, Ztrf)
!   
! 
!    open( unit = 4, file = 'el_erg', status = 'replace') 
!    do j = 1, dims
!        write(4, *) energias(j)
!    enddo
!    close(4)
!
!
!   endif 
!
!    if (el_or_hl == 2 ) then
!    call printaletters2('hl_ha1', 5, dims, dims, dims, hamiltoniana)
!    call printaletters2('hl_phi', 6, dims, dims, dims, phi)
!    call printaletters2('hl_ovm', 7, dims, dims, dims, ovlpm) 
!    
!    open( unit = 8, file = 'hl_erg', status = 'replace') 
!    do j = 1, dims
!        write(8, *) energias(j)
!    enddo
!    close(8)
!    
!    endif
!
!endif 


! --------------------------------------------------------

!------------------- CONDICOES INICIAIS NA BASE EXC ---------------
call rhosite_TO_rhoham(dims, phi, phi_transpose, rho_sites_in, rho_ham_in)
!------------------------------------------------------------------



!--------------- CALCULO DA PUREZA E ENERGIA DA PARTICULA ---------
!purity = 0.d0
!call gemm(rho_sites_in, rho_sites_in, rho_sites_in_squared)
!do i = 1, dims
!purity = purity + real(rho_sites_in_squared(i, i))
!enddo  

call particle_energy(dims, hamiltoniana, rho_sites_in, ParticleEnergy)

!-------- ENERGIA DO ELETRON ---------
if (el_or_hl == 1) then
write(79, "(60F20.8)") t, energiazeroel
write(80, "(60F20.8)") t, ParticleEnergy
!write(81, "(60F20.8)") t, purity
endif

!-------- ENERGIA DO BURACO ---------
if (el_or_hl == 2) then
write(82, "(60F20.8)") t, energiazerohl
write(83, "(60F20.8)") t, ParticleEnergy
!write(84, "(60F20.8)") t, purity
!print*, hamiltoniana
!print*, "omega é", sites_array(1, 1)%omega
endif
!-------------------------------------------------------------------



!------------------ TRANSFORMO AS CONDICOES INICIAIS PARA OS VETORES LINHAS DO ODE -----
call monta_y(dims, rho_ham_in, y) !condição inicial para resolver as ODE'S
!--------------------------------------------------------------------------------------


!--------------- CRIO A MATRIZ RW QUE É UTILIZADA PARA CALCULAR O ODE --------
t3 = omp_get_wtime()
call createRwMatrix(dims, frequency_matrix, RDtensor, RWMatrix) 
t4 = omp_get_wtime()
! ---------------------------------------------------------------------------



!y esta na forma (para d=2):  y(1) = rho(1, 1)  y(3) = real(rho(2, 1))
!                             y(2) = rho(2, 2)  y(4) = aimag(rho(2, 1))
!----------------------------------------------------------------------------------
  

!---------------- ESCREVO O  RESULTADO  ----------------


call printa_resultado(nstates, dims, number_file_ham, t, rho_ham_in)
call printa_resultado(nstates, dims, number_file_site, t, rho_sites_in)



if (el_or_hl == 1) then
  do i = 1, dim_el
    energ(i) = real(rho_sites_in(i, i)) * hamiltoniana(i, i)
  enddo

open(90, position = 'append')
open(70, position = 'append')

write(75, "(60F12.5)") t, sum(energ(:))
write(70, "(60F12.5)") t, HB_ev * sites_array(1, 1)%omega * thz_to_hz

!write(90, "(60F12.5)") t, real(rho_sites_in(3, 1)), real(rho_sites_in(4, 1)), real(rho_sites_in(3, 2)) &
!        , real(rho_sites_in(4, 2))

close(90)
close(70)
endif 




!--------------------------------------------------------------------------

! ----------------- CHAMO A SUBROTINA QUE CALCULA PARA CADA TEMPO TOUT AS EQUACOES DIFERENCIAIS

t5 = omp_get_wtime() 
call ode ( f0, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )
t6 = omp_get_wtime() 

!!print*, "Tempo para montar o tensor de redfield", t2 - t1
!!print*, "Tempo para montar a matrix Rw", t4 - t3
!!print*, "Tempo para calcular a ODE", t6 - t5


if ( iflag /= 2 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST - Fatal error!'
  write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
  stop
end if
!-----------------------------------------------------------------------------------------------



!------- TRANSFORMO O RESULTADO DO VETOR LINHA Y DO ODE PARA AS MATRIZES RHO --- 
call monta_rho(dims, y, rho_ham)
!----------------------- ESCREVO O OPERADOR DENSIDADE NA BASE DO SITIO----------------------
call rhoham_TO_rhosite(dims, phi, phi_transpose, rho_ham, rho_sites)
!-------------------------------------------------------------------------------------------



!-------- CALCULO AS POPULACOES DE CARGA EM CADA SITIO ----------
call rho_matrix_to_pop(nstates, dims, rho_sites, pop)
!-------------------------------------------------------------------------------------------------------------------


    

DEALLOCATE(phi, phi_transpose, multmtx, rho_sites_in, rho_ham_in, matriz_temporaria, pop_in, RDtensor, RWMatrix, energymtc, &
rho_sites_in_squared, ovlpm, rho_rec)
13 format(3es14.3E3)
return
end subroutine test



subroutine f0 ( t, y, yp )
!não definimos as condicoes iniciais aqui, apenas a funcao que vamos calcular a derivada. A condição inicial vem antes do call
implicit none
real*8 ::  t !tempo que vai entrar na equação diferencial
real*8 ::  y(neqn)  !funcao que estamos querendo, neste caso x
real*8 ::  yp(neqn) !derivada que estamos querendo calcular, por exemplo dx/dt = -x => sol: x = exp(-t)


call gemv(RWMatrix, y, yp)  !gemv calcula o produto da matriz RWMatrix com
                            !o vetor y resultando em yp => muito mais rapido

return
13 format(3es14.3E3)
end subroutine f0



end module testes
