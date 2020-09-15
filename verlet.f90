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

subroutine particle_energy(hamsize, hamiltoniano, rho, ParticleEnergy)
implicit none
integer, intent(in) :: hamsize
real*8, dimension(hamsize, hamsize), intent(in) :: hamiltoniano
complex*16, dimension(hamsize, hamsize), intent(in) :: rho
real*8, intent(out) :: ParticleEnergy

integer :: i 
complex*16, allocatable, dimension(:, :) :: energymtx


allocate(energymtx(hamsize, hamsize), source = (0.d0, 0.d0) ) 

call gemm(hamiltoniano, rho, energymtx)


ParticleEnergy = 0.d0
do i = 1, hamsize
ParticleEnergy = ParticleEnergy + real(energymtx(i, i))
enddo



deallocate(energymtx) 
return
end subroutine particle_energy        



subroutine kinect_energy(el_or_hl, rhosize, RadialVel, KinectEnergy)
implicit none
integer, intent(in) :: el_or_hl, rhosize
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

 
subroutine spring_energy(el_or_hl, rhosize, Radius, SpringEnergy)
implicit none
integer, intent(in) :: el_or_hl, rhosize
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



subroutine spring_force(el_or_hl, rhosize, SpringForce )!, SpringEnergy)
implicit none
integer, intent(in) :: el_or_hl
integer, intent(in) :: rhosize
real*8, dimension(:, :), allocatable, intent(out) :: SpringForce
!real*8, dimension(:, :), allocatable, intent(out) :: SpringEnergy

integer :: i, j 
integer :: nmstates, dims

if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if

if (el_or_hl == 1) then
nmstates = nmstates_el
dims = dim_el
endif


if (el_or_hl == 2) then
nmstates = nmstates_hl
dims = dim_hl
endif

ALLOCATE(SpringForce(nm_rows, nm_columns), source = 0.d0)
!ALLOCATE(SpringEnergy(nm_rows, nm_columns), source = 0.d0 ) 

do j = 1, nm_columns
  do i = 1, nm_rows
SpringForce(i, j) = sites_array(i, j)%mass*((sites_array(i, j)%omegazero)**2.0)*& 
        (sites_array(i, j)%radiuszero-sites_array(i, j)%radius)
  enddo
enddo


return
end subroutine spring_force



subroutine eletric_force(el_or_hl, rhosize, rho_el_real, rho_hl_real, elHam, hlHam, EletricForce) 
implicit none
integer, intent(in) :: el_or_hl
integer, intent(in) :: rhosize
real*8, dimension(dim_el, dim_el), intent(in) :: rho_el_real
real*8, dimension(dim_hl, dim_hl), intent(in) :: rho_hl_real
real*8, dimension(dim_el, dim_el), intent(in) :: elHam
real*8, dimension(dim_hl, dim_hl), intent(in) :: hlHam
real*8, dimension(:, :), allocatable, intent(out) :: EletricForce


real*8, dimension(:, :), allocatable :: EletricForceDiagonal, EletricForceNonDiagOne
real*8, dimension(:, :), allocatable :: EletricForceNonDiagTwo
real*8, dimension(6) :: Qnumbers_energy
real*8, dimension(:), allocatable :: ForceStates
integer :: nmstates
real*8 :: dims 
real*8, dimension(:, :), allocatable :: derivativeTerm
real*8, dimension(:, :), allocatable :: hamTimesRhoEl, hamTimesRhoHl
real*8, dimension(:, :), allocatable :: dermtx, ElRhoDerTerm, HlRhoDerTerm

real*8, dimension(:, :), allocatable :: MtxElForceNonDiagOne, MtxElForceNonDiagTwo
real*8, dimension(:, :), allocatable :: MtxHlForceNonDiagOne, MtxHlForceNonDiagTwo

real*8, dimension(:, :), allocatable :: ElHamRhoDer, HlHamRhoDer


real*8 :: ElForceNonDiagOne, ElForceNonDiagTwo
real*8 :: HlForceNonDiagOne, HlForceNonDiagTwo


if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if

if (el_or_hl == 1) then
nmstates = nmstates_el
dims = dim_el
endif


if (el_or_hl == 2) then
nmstates = nmstates_hl
dims = dim_hl
endif


Qnumbers_energy = (/ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 /)




ALLOCATE( EletricForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( EletricForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 ) 

ALLOCATE( ElHamRhoDer(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( HlHamRhoDer(nm_rows, nm_columns), source = 0.d0 ) 

ALLOCATE( MtxElForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( MtxElForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 )

ALLOCATE( MtxHlForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( MtxHlForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 )


ALLOCATE( ElRhoDerTerm(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( HlRhoDerTerm(nm_rows, nm_columns), source = 0.d0 )

ALLOCATE( hamTimesRhoEl(dim_el, dim_el), source = 0.d0 ) 
ALLOCATE( hamTimesRhoHl(dim_hl, dim_hl), source = 0.d0 ) 

ALLOCATE( EletricForce(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( EletricForceDiagonal(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( derivativeTerm(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( ForceStates(nmstates), source = 0.d0 )


do j = 1, nm_columns
  do i = 1, nm_rows
  derivativeTerm(i, j) = - (( 4.d0 * hbar ) / ( me * (sites_array(i, j)%radius)**3.0 ))
  enddo
enddo


call calculate_derivative_matrix(nmstates_el, dim_el, dermtx)

print*, "Matriz derivada"
call print_mat2(dermtx, dim_el, dim_el) 

! -------- parte diagonal --------- 
l = 1
do j = 1, nm_columns
  do i = 1, nm_rows

  do k = 1, nmstates
  ForceStates(k) = rho_el_real(l, l) * Qnumbers_energy(k) * hbar + rho_hl_real(l, l) * Qnumbers_energy(k) * hbar
  l = l + 1
  enddo



   EletricForceDiagonal(i, j) = derivativeTerm(i, j) * sum(ForceStates(:))
  enddo
enddo
! ------------------------------


! ----- primeira parte não diagonal --------

ElRhoDerTerm = matmul( rho_el_real, dermtx)   
HlRhoDerTerm = matmul( rho_hl_real, dermtx)


ElForceNonDiagOne = 0.d0
HlForceNonDiagOne = 0.d0

l = 1
do j = 1, nm_columns
  do i = 1, nm_rows

    do k = 1, nmstates_el
    ElForceNonDiagOne = ElForceNonDiagOne + ElRhoDerTerm(l, l)
    HlForceNonDiagOne = HlForceNonDiagOne + HlRhoDerTerm(l, l)
    l = l + 1
    enddo
 
  MtxElForceNonDiagOne(i, j) = ElForceNonDiagOne
  MtxHlForceNonDiagOne(i, j) = HlForceNonDiagOne
  ElForceNonDiagOne = 0.d0
  HlForceNonDiagOne = 0.d0 

  enddo
enddo

! --------------------------------------------

! ----- segunda parte não diagonal ----------


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




!------------------------------------------


do j = 1, nm_columns
  do i = 1, nm_rows
EletricForceNonDiagOne(i, j) =  2.d0 * ev_to_joule * sites_array(i, j)%t * MtxElForceNonDiagOne(i, j) & 
                             +  2.d0 * ev_to_joule * sites_array(i, j)%t * MtxHlForceNonDiagOne(i, j) 


!EletricForceNonDiagTwo(i, j) = 2.d0 * MtxElForceNonDiagTwo(i, j) + 2.d0 * MtxHlForceNonDiagTwo(i, j)
  enddo
enddo


do j = 1, nm_columns
  do i = 1, nm_rows
EletricForce(i, j) = -EletricForceDiagonal(i, j) -EletricForceNonDiagOne(i, j) !+ EletricForceNonDiagTwo(i, j)
  enddo
enddo


EletricForce = eletricCte * EletricForce

print*, "-Eletric force diagonal SITIO 1", -EletricForceDiagonal(1, 1)
print*, "-Eletric force non diagonal 1 SITIO 1", -EletricForceNonDiagOne(1, 1)
!print*, "+Eletric force non diagonal2 part1", +EletricForceNonDiagTwo(1, 1)


print*, "-Eletric force diagonal SITIO 2", -EletricForceDiagonal(1, 2)
print*, "-Eletric force non diagonal 1 SITIO 2", -EletricForceNonDiagOne(1, 2)
!print*, "+Eletric force non diagonal2 part2", +EletricForceNonDiagTwo(1, 2)



!-------------------------------




deallocate(derivativeTerm, ForceStates, EletricForceDiagonal, EletricForceNonDiagOne, EletricForceNonDiagTwo, ElRhoDerTerm, HlRhoDerTerm &
        , MtxElForceNonDiagOne, MtxElForceNonDiagTwo, MtxHlForceNonDiagOne, MtxHlForceNonDiagTwo, ElHamRhoDer, HlHamRhoDer)
return
end subroutine eletric_force

subroutine velocity_verlet(BWRadius, BWRadialVel, BWEletricForce, BWSpringForce, time,  delta_t, rho_el_real, & 
                rho_hl_real, elHam, hlHam, FWRadius, FWRadialVel, FWEletricForce, FWSpringForce, SpringEnergy, KinectEnergy )
implicit none
real*8, intent(in), dimension(nm_rows, nm_columns) :: BWRadius, BWRadialVel, BWEletricForce, BWSpringForce
real*8, intent(in) :: time, delta_t
real*8, intent(in), dimension(dim_el, dim_el) :: rho_el_real
real*8, intent(in), dimension(dim_hl, dim_hl) :: rho_hl_real
real*8, intent(in), dimension(dim_el, dim_el) :: elHam
real*8, intent(in), dimension(dim_hl, dim_hl) :: hlHam 
real*8, allocatable, dimension(:, :), intent(out) :: FWRadius, FWRadialVel, FWEletricForce, FWSpringForce
real*8, intent(out) :: SpringEnergy, KinectEnergy
real*8 :: rn
real*8 :: BWDamp(nm_rows, nm_columns) 


real*8, allocatable, dimension(:, :) :: BWForce, FWForce

!BWRadius -> BackWardRadius (raio anterior).    FWRadius -> ForWardRadius (raio posterior).
integer :: k, i, j 


ALLOCATE( FWRadius(nm_rows, nm_columns), source = 0.d0 ) 
ALLOCATE( FWRadialVel(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( BWForce(nm_rows, nm_columns), source = 0.d0 )
ALLOCATE( FWForce(nm_rows, nm_columns), source = 0.d0 )
!ALOCO FWEletricForce na subrotina eletric_force

     

!forall (i = 1 : nm_rows, j = 1 : nm_columns) BWDamp(i, j) = - Dterm * BWRadialVel(i, j) 

forall (i = 1:nm_rows, j = 1:nm_columns) BWForce(i, j) =  BWEletricForce(i, j) + cteForce*BWSpringForce(i, j) !+ BWDamp(i, j) 







call kinect_energy(1, dim_el, BWRadialVel, KinectEnergy) 
call spring_energy(1, dim_el, BWRadius, SpringEnergy)


! ==============================================


    do j = 1, nm_columns
      do i = 1, nm_rows


      FWRadius(i, j) = BWRadius(i, j) + ( BWRadialVel(i, j) *  (delta_t * 1.d-12)) + (BWForce(i, j) * &
                              ( ( (delta_t * 1.d-12))**2.d0 / (2.d0 * sites_array(i,j)%mass ) ) )



 !    !ATUALIZO O RAIO E OMEGA
     sites_array(i, j)%radius = FWRadius(i, j)
     sites_array(i, j)%omega = ( 2.d0 * hbar / ( me  * (sites_array(i, j)%radius)**2.0 ) ) * hz_to_thz

     call eletric_force(1, dim_el, rho_el_real, rho_hl_real, elHam, hlHam, FWEletricForce)
     call spring_force(1, dim_el, FWSpringForce) 
    
     FWForce(i, j) = FWEletricForce(i, j) + cteForce*FWSpringForce(i, j)

    
     !FWRadialVel(i, j) = ( FWRadius(i, j) - BWRadius(i, j) ) / (2.d0 * delta_t * 1.d-12)   
  
     
     FWRadialVel(i, j) = BWRadialVel(i, j) + ( BWForce(i, j) + FWForce(i, j) ) * ( (delta_t * 1.d-12) / (2.d0 &
             * sites_array(i, j)%mass ) ) 


     enddo
   enddo


  
    
    
     !ATUALIZO A VELOCIDADE RADIAL 
    

  
   sites_array(:, :)%radial_vel = FWRadialVel(:, :) 



 


deallocate(BWForce, FWForce) 
return
end subroutine velocity_verlet




end module verlet_m
