module verlet_m
use f95_precision
use lapack95
use blas95
use parameters_m
use constants_m 
use types_m 
use system_hamiltonian_m
use functions_m 

contains 

subroutine particle_energy(hamiltoniano, rho, ParticleEnergy)
implicit none
real*8, dimension(dim_el, dim_el), intent(in) :: hamiltoniano
complex*16, dimension(dim_el, dim_el), intent(in) :: rho
real*8, intent(out) :: ParticleEnergy

integer :: i 
complex*16, allocatable, dimension(:, :) :: energymtx


allocate(energymtx(dim_el, dim_el), source = (0.d0, 0.d0) ) 

call gemm(hamiltoniano, rho, energymtx)


ParticleEnergy = 0.d0
do i = 1, dim_el
ParticleEnergy = ParticleEnergy + real(energymtx(i, i))
enddo



deallocate(energymtx) 
return
end subroutine particle_energy        



subroutine kinect_energy(RadialVel, KinectEnergy)
implicit none
real*8, intent(in), dimension(nm_rows, nm_columns) :: RadialVel
real*8, intent(out) :: KinectEnergy
real*8, dimension(:, :), allocatable :: kinect_energy_sites


allocate(kinect_energy_sites(nm_rows, nm_columns), source = 0.d0)


do j = 1, nm_columns
  do i = 1, nm_rows
  kinect_energy_sites(i, j) = HALF * sites_array(i, j)%mass * RadialVel(i, j)**2.0
  enddo
enddo

KinectEnergy = sum(kinect_energy_sites(:, :)) * joule_to_ev

!print*, KinectEnergy


deallocate(kinect_energy_sites)
return
end subroutine kinect_energy

 
subroutine spring_energy(Radius, SpringEnergy)
implicit none
real*8, dimension(nm_rows, nm_columns), intent(in) :: Radius
real*8, intent(out) :: SpringEnergy


real*8, dimension(:, :), allocatable :: spring_energy_sites


ALLOCATE( spring_energy_sites(nm_rows, nm_columns), source = 0.d0 ) 


do j = 1, nm_columns
  do i = 1, nm_rows
spring_energy_sites(i, j) = HALF*sites_array(i, j)%mass*((sites_array(i, j)%omegazero)**2.0) * &
                      (Radius(i, j) - sites_array(i, j)%radiuszero)**2.0 
  enddo
enddo


SpringEnergy = sum(spring_energy_sites(:, :)) * joule_to_ev


deallocate(spring_energy_sites)
return
end subroutine spring_energy        



subroutine spring_force(SpringForce)
implicit none
real*8, dimension(:, :), allocatable, intent(out) :: SpringForce

integer :: i, j 


ALLOCATE(SpringForce(nm_rows, nm_columns), source = 0.d0)

do j = 1, nm_columns
  do i = 1, nm_rows
SpringForce(i, j) = sites_array(i, j)%mass*((sites_array(i, j)%omegazero)**2.0)*& 
        (sites_array(i, j)%radiuszero-sites_array(i, j)%radius)
  enddo
enddo


return
end subroutine spring_force



subroutine eletric_force(pl, rho_el_real, elHam, EletricForce) 
implicit none
integer, intent(in) :: pl
real*8, dimension(dim_el, dim_el), intent(in) :: rho_el_real
real*8, dimension(dim_el, dim_el), intent(in) :: elHam !hamiltoniano eletron
real*8, dimension(:, :), allocatable, intent(out) :: EletricForce


real*8, dimension(:, :), allocatable :: EletricForceDiagonal, EletricForceNonDiagOne
real*8, dimension(:, :), allocatable :: EletricForceNonDiagTwo
real*8, dimension(6) :: Qnumbers_energy
real*8, dimension(:), allocatable :: ForceStates
real*8, dimension(:, :), allocatable :: derivativeTerm
real*8, dimension(:, :), allocatable :: hamTimesRhoEl
real*8, dimension(:, :), allocatable :: dermtx, ElRhoDerTerm

real*8, dimension(:, :), allocatable :: MtxElForceNonDiagOne, MtxElForceNonDiagTwo


real*8, dimension(:, :), allocatable :: ElHamRhoDer


real*8 :: ElForceNonDiagOne, ElForceNonDiagTwo



Qnumbers_energy = (/ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 /)




ALLOCATE( EletricForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( EletricForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 ) 

ALLOCATE( ElHamRhoDer(nm_rows, nm_columns), source = 0.d0 ) 

ALLOCATE( MtxElForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( MtxElForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 )



ALLOCATE( ElRhoDerTerm(nm_rows, nm_columns), source = 0.d0 )

ALLOCATE( hamTimesRhoEl(dim_el, dim_el), source = 0.d0 ) 

ALLOCATE( EletricForce(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( EletricForceDiagonal(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( derivativeTerm(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( ForceStates(nmstates_el), source = 0.d0 )


!========== CALCULO DAS DERIVADAS DW/DR ===========================
derivativeTerm(:, :) = - (( 4.d0 * hbar ) / ( me * (sites_array(:, :)%radius)**3.0 ))
!==================================================================

!=========== CALCULO DA MATRIZ DERIVA =========
call calculate_derivative_matrix(dermtx)
!==============================================

if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
  print*, "Matriz derivada"
  call print_mat2(dermtx, dim_el, dim_el) 
endif


!======= PARTE DIAGONAL DA FORCA - DEPENDENTE DAS POPULACOES =============
l = 1
do j = 1, nm_columns
  do i = 1, nm_rows

  do k = 1, nmstates_el
  ForceStates(k) = rho_el_real(l, l) * Qnumbers_energy(k) * hbar
  l = l + 1
  enddo

   EletricForceDiagonal(i, j) = derivativeTerm(i, j) * sum(ForceStates(:))
  enddo
enddo
!=========================================================================


!===== PRIMEIRA PARTE NAO DIAGONAL - DEPENDENTE APENAS DAS COERÊNCIAS =====

ElRhoDerTerm = matmul( rho_el_real, dermtx)   


ElForceNonDiagOne = 0.d0

l = 1
do j = 1, nm_columns
  do i = 1, nm_rows

    do k = 1, nmstates_el
    ElForceNonDiagOne = ElForceNonDiagOne + ElRhoDerTerm(l, l)
    l = l + 1
    enddo
 
  MtxElForceNonDiagOne(i, j) = ElForceNonDiagOne
  ElForceNonDiagOne = 0.d0

  enddo
enddo

!==========================================================================

!========== SEGUNDA PARTE NÃO DIAGONAL - DEPENDENTE DOS ACOPLAMENTOS E COERÊNCIAS ========


!hamTimesRhoEl = matmul( rho_el_real, elHam * ev_to_joule)
!hamTimesRhoHl = matmul( rho_hl_real, hlHam * ev_to_joule)



!ElHamRhoDer = matmul( hamTimesRhoEl, dermtx) 
!HlHamRhoDer = matmul( hamTimesRhoHl, dermtx) 


!ElForceNonDiagTwo = 0.d0
!HlForceNonDiagTwo = 0.d0

!l = 1
!do j = 1, nm_columns
!  do i = 1, nm_rows

!    do k = 1, nmstates_el
!    ElForceNonDiagTwo = ElForceNonDiagTwo + ElHamRhoDer(l, l)
!    HlForceNonDiagTwo = HlForceNonDiagTwo + HlHamRhoDer(l, l)
!    l = l + 1
!    enddo
 
!  MtxElForceNonDiagTwo(i, j) = ElForceNonDiagTwo
!  MtxHlForceNonDiagTwo(i, j) = HlForceNonDiagTwo
!  ElForceNonDiagTwo = 0.d0
!  HlForceNonDiagTwo = 0.d0 

!  enddo
!enddo

!====================================================================


EletricForceNonDiagOne(:, :) =  2.d0 * ev_to_joule * sites_array(:, :)%t * MtxElForceNonDiagOne(:, :)  
!EletricForceNonDiagTwo(i, j) = 2.d0 * MtxElForceNonDiagTwo(i, j) + 2.d0 * MtxHlForceNonDiagTwo(i, j)


! ================== FORCA ELETRICA TOTAL =====================
EletricForce(:, :) = - EletricForceDiagonal(:, :) - EletricForceNonDiagOne(:, :) !+ EletricForceNonDiagTwo(i, j)
!==============================================================


if ( pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
  print*, "-Eletric force diagonal SITIO 1", -EletricForceDiagonal(1, 1)
  print*, "-Eletric force non diagonal 1 SITIO 1", -EletricForceNonDiagOne(1, 1)
  !print*, "+Eletric force non diagonal2 part1", +EletricForceNonDiagTwo(1, 1)


  print*, "-Eletric force diagonal SITIO 2", -EletricForceDiagonal(1, 2)
  print*, "-Eletric force non diagonal 1 SITIO 2", -EletricForceNonDiagOne(1, 2)
  !print*, "+Eletric force non diagonal2 part2", +EletricForceNonDiagTwo(1, 2)
endif 


!-------------------------------




deallocate(derivativeTerm, ForceStates, EletricForceDiagonal, EletricForceNonDiagOne, EletricForceNonDiagTwo, ElRhoDerTerm &
          , MtxElForceNonDiagOne, MtxElForceNonDiagTwo, ElHamRhoDer)
return
end subroutine eletric_force

subroutine velocity_verlet(pl, BWRadius, BWRadialVel, BWEletricForce, BWSpringForce, delta_t, rho_el_real, & 
                elHam, FWRadius, FWRadialVel, FWEletricForce, FWSpringForce, SpringEnergy, KinectEnergy )
implicit none
integer, intent(in) :: pl 
real*8, intent(in), dimension(nm_rows, nm_columns) :: BWRadius, BWRadialVel, BWEletricForce, BWSpringForce
real*8, intent(in) :: delta_t
real*8, intent(in), dimension(dim_el, dim_el) :: rho_el_real
real*8, intent(in), dimension(dim_el, dim_el) :: elHam
real*8, allocatable, dimension(:, :), intent(out) :: FWRadius, FWRadialVel, FWEletricForce, FWSpringForce
real*8, intent(out) :: SpringEnergy, KinectEnergy


real*8, allocatable, dimension(:, :) :: BWForce, FWForce

!BWRadius -> BackWardRadius (raio anterior).    FWRadius -> ForWardRadius (raio posterior).
integer :: k, i, j 


ALLOCATE( FWRadius(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( FWRadialVel(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( BWForce(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( FWForce(nm_rows, nm_columns), source = 0.d0 )
!ALOCO FWEletricForce na subrotina eletric_force

     
forall (i = 1:nm_rows, j = 1:nm_columns) BWForce(i, j) =  BWEletricForce(i, j) + BWSpringForce(i, j) 



!=============== CALCULO DA ENERGIA CINETICA E POTENCIAL BASEADO NA FORCA EM t =====
call kinect_energy(BWRadialVel, KinectEnergy) 
call spring_energy(BWRadius, SpringEnergy)
!===================================================================================

    !========== CALCULO DOS NOVOS RAIOS BASEDOS NA FORCA EM t ======================

    FWRadius(:, :) = BWRadius(:, :) + ( BWRadialVel(:, :) *  (delta_t * 1.d-12)) + (BWForce(:, :) * &
                              ( ( (delta_t * 1.d-12))**2.d0 / (2.d0 * sites_array(:, :)%mass ) ) )

    !===============================================================================


    !========= ATUALIZO O RAIO E OMEGA ==============================================
     sites_array(:, :)%radius = FWRadius(:, :)
     sites_array(:, :)%omega = ( 2.d0 * hbar / ( me  * (sites_array(:, :)%radius)**2.0 ) ) * hz_to_thz
    !================================================================================

    !======== CALCULO DAS NOVAS FORCAS EM t + delta t ===============================
     call eletric_force(pl, rho_el_real, elHam, FWEletricForce)
     call spring_force(FWSpringForce) 
    !===============================================================================
    
    !======= TOTAL FORWARD FORCE ==================================================
     FWForce(:, :) = FWEletricForce(:, :) + FWSpringForce(:, :)
    !==============================================================================

    
     !FWRadialVel(i, j) = ( FWRadius(i, j) - BWRadius(i, j) ) / (2.d0 * delta_t * 1.d-12)   
  
     !====== NOVA VELOCIDADE RADIAL ===================================================
     FWRadialVel(:, :) = BWRadialVel(:, :) + ( BWForce(:, :) + FWForce(:, :) ) * ( (delta_t * 1.d-12) / (2.d0 &
             * sites_array(:, :)%mass ) ) 

     !=================================================================================

    
    
     !===== ATUALIZO A VELOCIDADE RADIAL  ===========================================
      sites_array(:, :)%radial_vel = FWRadialVel(:, :) 
     !==============================================================================


 


deallocate(BWForce, FWForce) 
return
end subroutine velocity_verlet




end module verlet_m
