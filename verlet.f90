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

    ! args
    real*8,     intent(in) :: hamiltoniano(dim_el, dim_el)
    complex*16, intent(in) :: rho(dim_el, dim_el)
    real*8,     intent(out):: ParticleEnergy
   
    ! local  
    integer :: i 
    complex*16, allocatable :: energymtx(:,:)

    allocate(energymtx(dim_el, dim_el), source = (0.d0, 0.d0) ) 
    
    call gemm(hamiltoniano, rho, energymtx)
    
    ParticleEnergy = 0.d0
    do i = 1, dim_el
       ParticleEnergy = ParticleEnergy + real(energymtx(i, i))
    enddo
    
    deallocate(energymtx) 
end subroutine particle_energy        
    
    
    
subroutine kinect_energy(RadialVel, KinectEnergy)
    implicit none

    ! args
    real*8, intent(in)  :: RadialVel(nm_rows, nm_columns)
    real*8, intent(out) :: KinectEnergy

    ! local
    real*8,  allocatable :: kinect_energy_sites(:,:)
    
    allocate(kinect_energy_sites(nm_rows, nm_columns), source = 0.d0)

    do j = 1, nm_columns
       do i = 1, nm_rows
          kinect_energy_sites(i, j) = HALF * sites_array(i, j)%mass * RadialVel(i, j)**2.0
       enddo
    enddo
    
    KinectEnergy = sum(kinect_energy_sites(:, :)) * joule_to_ev
    
    deallocate(kinect_energy_sites)
end subroutine kinect_energy

 
subroutine spring_energy(Radius, SpringEnergy)
    implicit none

    ! args
    real*8, intent(in)  :: Radius(nm_rows, nm_columns)
    real*8, intent(out) :: SpringEnergy
    
    ! local  
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
end subroutine spring_energy        



subroutine spring_force(SpringForce)
    implicit none

    ! args
    real*8, allocatable, intent(out) :: SpringForce(:,:)
    
    ! local 
    integer :: i, j 
    
    ALLOCATE(SpringForce(nm_rows, nm_columns), source = 0.d0)
    
    do j = 1, nm_columns
       do i = 1, nm_rows
          SpringForce(i, j) = sites_array(i, j)%mass*((sites_array(i, j)%omegazero)**2.0)*& 
            (sites_array(i, j)%radiuszero-sites_array(i, j)%radius)
       enddo
    enddo
end subroutine spring_force



subroutine eletric_force(pl, rho_el_real, elHam, EletricForce) 
    implicit none

    ! args
    integer,             intent(in)  :: pl
    real*8,              intent(in)  :: rho_el_real(dim_el, dim_el)
    real*8,              intent(in)  :: elHam (dim_el, dim_el)
    real*8, allocatable, intent(out) :: EletricForce(:,:)
    
    ! local 
    real*8,  allocatable :: EletricForceDiagonal(:,:), EletricForceNonDiagOne(:,:)
    real*8,  allocatable :: EletricForceNonDiagTwo(:,:)
    real*8,              :: Qnumbers_energy(6)
    real*8,  allocatable :: ForceStates(:)
    real*8,  allocatable :: derivativeTerm(:,:)
    real*8,  allocatable :: hamTimesRhoEl(:,:)
    real*8,  allocatable :: dermtx, ElRhoDerTerm(:,:)
    real*8,  allocatable :: MtxElForceNonDiagOne(:,:), MtxElForceNonDiagTwo(:,:)
    real*8,  allocatable :: ElHamRhoDer(:,:)
    real*8               :: ElForceNonDiagOne, ElForceNonDiagTwo
    
    Qnumbers_energy = [ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 ]

    allocate( EletricForceNonDiagOne(nm_rows, nm_columns), source = 0.d0 )
    allocate( EletricForceNonDiagTwo(nm_rows, nm_columns), source = 0.d0 ) 
    allocate( ElHamRhoDer(nm_rows, nm_columns)           , source = 0.d0 ) 
    allocate( MtxElForceNonDiagOne(nm_rows, nm_columns)  , source = 0.d0 )
    allocate( MtxElForceNonDiagTwo(nm_rows, nm_columns)  , source = 0.d0 )
    allocate( ElRhoDerTerm(nm_rows, nm_columns)          , source = 0.d0 )
    allocate( hamTimesRhoEl(dim_el, dim_el)              , source = 0.d0 ) 
    allocate( EletricForce(nm_rows, nm_columns)          , source = 0.d0 ) 
    allocate( EletricForceDiagonal(nm_rows, nm_columns)  , source = 0.d0 )
    allocate( derivativeTerm(nm_rows, nm_columns)        , source = 0.d0 ) 
    allocate( ForceStates(nmstates_el)                   , source = 0.d0 )
    
    !========== CALCULO DAS DERIVADAS DW/DR ===========================
    derivativeTerm(:,:) = - (( 4.d0 * hbar ) / ( me * (sites_array(:, :)%radius)**3.0 ))
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
   EletricForceNonDiagOne(:, :) =  2.d0 * ev_to_joule * sites_array(:, :)%t * MtxElForceNonDiagOne(:, :)  
   
   ! ================== FORCA ELETRICA TOTAL =====================
   EletricForce(:, :) = - EletricForceDiagonal(:, :) - EletricForceNonDiagOne(:, :) !+ EletricForceNonDiagTwo(i, j)
   !==============================================================
   
   
   if ( pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
     print*, "-Eletric force diagonal SITIO 1", -EletricForceDiagonal(1, 1)
     print*, "-Eletric force non diagonal 1 SITIO 1", -EletricForceNonDiagOne(1, 1)
     print*, "-Eletric force diagonal SITIO 2", -EletricForceDiagonal(1, 2)
     print*, "-Eletric force non diagonal 1 SITIO 2", -EletricForceNonDiagOne(1, 2)
   endif 
   
   deallocate(derivativeTerm, ForceStates, EletricForceDiagonal, EletricForceNonDiagOne, EletricForceNonDiagTwo, ElRhoDerTerm &
             , MtxElForceNonDiagOne, MtxElForceNonDiagTwo, ElHamRhoDer)
end subroutine eletric_force



subroutine velocity_verlet(pl, BWRadius, BWRadialVel, BWEletricForce, BWSpringForce, delta_t, rho_el_real, & 
                elHam, FWRadius, FWRadialVel, FWEletricForce, FWSpringForce, SpringEnergy, KinectEnergy )
    implicit none

    ! args
    integer,               intent(in)  :: pl 
    real*8,                intent(in), dimension(nm_rows, nm_columns) :: BWRadius, BWRadialVel, BWEletricForce, BWSpringForce
    real*8,                intent(in)  :: delta_t
    real*8,                intent(in)  :: rho_el_real(dim_el, dim_el)
    real*8,                intent(in)  :: elHam(dim_el, dim_el)
    real*8,  allocatable,  intent(out) :: FWRadius(:,:), FWRadialVel(:,:), FWEletricForce(:,:), FWSpringForce(:,:)
    real*8,                intent(out) :: SpringEnergy, KinectEnergy
    
    ! local  
    real*8, allocatable :: BWForce(:,:), FWForce(:,:)
    integer             :: k, i, j 
    !BWRadius -> BackWardRadius (raio anterior).    FWRadius -> ForWardRadius (raio posterior).
    
    allocate( FWRadius(nm_rows, nm_columns)    , source = 0.d0 ) 
    allocate( FWRadialVel(nm_rows, nm_columns) , source = 0.d0 )
    allocate( BWForce(nm_rows, nm_columns)     , source = 0.d0 )
    allocate( FWForce(nm_rows, nm_columns)     , source = 0.d0 )
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
    
    !====== NOVA VELOCIDADE RADIAL ===================================================
    FWRadialVel(:, :) = BWRadialVel(:, :) + ( BWForce(:, :) + FWForce(:, :) ) * ( (delta_t * 1.d-12) / (2.d0 &
            * sites_array(:, :)%mass ) ) 
    !=================================================================================
    
    !===== ATUALIZO A VELOCIDADE RADIAL  ===========================================
     sites_array(:, :)%radial_vel = FWRadialVel(:, :) 
    !==============================================================================

    deallocate(BWForce, FWForce) 
end subroutine velocity_verlet

end module verlet_m
