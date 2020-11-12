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

subroutine particle_energy(h_mtx, rho, ParticleEnergy)
    implicit none

    ! args
    real*8,     intent(in) :: h_mtx(d_el, d_el)
    complex*16, intent(in) :: rho(d_el, d_el)
    real*8,     intent(out):: ParticleEnergy
   
    ! local  
    integer :: i 
    complex*16, allocatable :: MtxEnergy(:,:)

    allocate(MtxEnergy(d_el, d_el), source = (0.d0, 0.d0) ) 
    
    call gemm(h_mtx, rho, MtxEnergy)
    
    ParticleEnergy = 0.d0
    do i = 1, d_el
       ParticleEnergy = ParticleEnergy + real(MtxEnergy(i, i))
    enddo
    
    deallocate(MtxEnergy) 
end subroutine particle_energy        
    
    
    
subroutine kinect_energy(VelSite, K_Energy)
    implicit none

    ! args
    real*8, intent(in)  :: VelSite(nr, nc)
    real*8, intent(out) :: K_Energy

    ! local
    real*8,  allocatable :: K_EnergySites(:,:)
    
    allocate(K_EnergySites(nr, nc), source = 0.d0)

    K_EnergySites = HALF * site%mass * VelSite * VelSite
    
    K_Energy = sum(K_EnergySites(:,:)) * joule_to_ev
    
    deallocate(K_EnergySites)
end subroutine kinect_energy

 

subroutine spring_energy(RadiusSite, V_Energy)
    implicit none

    ! args
    real*8, intent(in)  :: RadiusSite(nr, nc)
    real*8, intent(out) :: V_Energy
    
    ! local  
    real*8, allocatable :: V_EnergySites(:,:)
    real*8, allocatable :: omega0Squared(:,:)
    real*8, allocatable :: raDiff(:,:)
    real*8, allocatable :: raDiffSquared(:,:)

    
    
    allocate( V_EnergySites(nr, nc),         source = 0.d0 ) 
    allocate( omega0Squared(nr, nc),         source = 0.d0 ) 
    allocate( raDiff(nr, nc),                source = 0.d0 ) 
    allocate( raDiffSquared(nr, nc),         source = 0.d0 ) 
    
    
    omega0Squared = freqMode * freqMode
    raDiff        = RadiusSite - site%radius0
    raDiffSquared      = raDiff * raDiff


    V_EnergySites = HALF * site%mass * omega0Squared * raDiffSquared
    
    V_Energy = sum(V_EnergySites(:, :)) * joule_to_ev
    
    deallocate(V_EnergySites, omega0Squared, raDiff, raDiffSquared)
end subroutine spring_energy        



subroutine spring_force(Vforce)
    implicit none

    ! args
    real*8, allocatable, intent(out) :: Vforce(:,:)
    
    ! local 
    integer :: i, j 
    real*8, allocatable :: omega0Squared(:,:), raDiff(:,:)
    
    
    ALLOCATE( Vforce(nr, nc),        source = 0.d0) !Potential Force
    ALLOCATE( omega0Squared(nr, nc), source = 0.d0)
    ALLOCATE( raDiff(nr, nc),        source = 0.d0) !Radius Difference

    omega0Squared = freqMode * freqMode
    raDiff        = site%radius0 - site%radius

    Vforce = site%mass * omega0Squared * raDiff

    deallocate(omega0Squared, raDiff)
end subroutine spring_force



subroutine eletric_force(pl, rhoElReal, time, elHam, Eforce) 
    implicit none

    ! args
    integer,             intent(in)  :: pl
    real*8,              intent(in)  :: rhoElReal(d_el, d_el)
    real*8,              intent(in)  :: time 
    real*8,              intent(in)  :: elHam (d_el, d_el)
    real*8, allocatable, intent(out) :: Eforce(:,:)
    
    ! local 
    real*8,  allocatable :: MtxEforceDiagonal(:,:)
    real*8,  allocatable :: ForceStates(:)
    real*8,  allocatable :: derTerm(:,:)
    real*8,  allocatable :: DerMtx(:,:), RhoTimesDer(:,:)
    real*8,  allocatable :: MtxEforceND1(:,:)
    real*8               :: EforceND1
    


    allocate( MtxEforceND1(nr, nc)      , source = 0.d0 )
    allocate( RhoTimesDer(nr, nc)       , source = 0.d0 )
    allocate( Eforce(nr, nc)            , source = 0.d0 ) 
    allocate( MtxEforceDiagonal(nr, nc) , source = 0.d0 )
    allocate( derTerm(nr, nc)           , source = 0.d0 ) 
    allocate( ForceStates(ns_el)        , source = 0.d0 )
    
    !========== CALCULO DAS DERIVADAS DW/DR ===========================
    derTerm = - (( 4.d0 * hbar ) / ( me * (site%radius)**3.0 ))
    !==================================================================
    
    !=========== CALCULO DA MATRIZ DERIVA =========
    call build_derivative_matrix(DerMtx)
    !==============================================
    
    if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
      print*, "Matriz derivada"
      call print_mat2(DerMtx, d_el, d_el) 
    endif
    
    !======= PARTE DIAGONAL DA FORCA - DEPENDENTE DAS POPULACOES =============
    l = 1
    do j = 1, nc
       do i = 1, nr
    
       do k = 1, ns_el
          ForceStates(k) = rhoElReal(l, l) * hbar * Qn_erg(k)
          l = l + 1
       enddo
    
       MtxEforceDiagonal(i, j) = derTerm(i, j) * sum(ForceStates(:))
       enddo
    enddo
    !=========================================================================
    

   !===== PRIMEIRA PARTE NAO DIAGONAL - DEPENDENTE APENAS DAS COERÊNCIAS =====
    
   RhoTimesDer = matmul(rhoElReal, DerMtx)   
    
   EforceND1 = 0.d0 !Eletric Force Non Diagonal One
    
   l = 1
   do j = 1, nc
     do i = 1, nr
    
         do k = 1, ns_el
            EforceND1 = EforceND1 + RhoTimesDer(l, l)
            l = l + 1
         enddo
     
       MtxEforceND1(i, j) = EforceND1
       EforceND1 = 0.d0
    
       enddo
    enddo
    
   !==========================================================================
   MtxEforceND1 =  2.d0 * ev_to_joule * site%t * MtxEforceND1
    

   ! ================== FORCA ELETRICA TOTAL =====================
   Eforce = - MtxEforceDiagonal - MtxEforceND1
   !==============================================================

 
   deallocate(derTerm, ForceStates, MtxEforceDiagonal, RhoTimesDer, MtxEforceND1)
end subroutine eletric_force


subroutine calculate_temperature(K_erg, insTemp)
    implicit none
    !args
    real*8,                intent(in)  :: K_erg   !energia cinetica
    real*8,                intent(out) :: insTemp !temperatura instantanea


    insTemp = ( (2.d0/3.d0) * K_erg) / (float(nsites) * kboltz)   !temp instantanea variavel global


end subroutine calculate_temperature



subroutine calculate_BRDscale(dt, insTemp, BRDscale) 
    implicit none

    !args
    real*8,                intent(in)  :: dt, insTemp
    real*8,                intent(out) :: BRDscale

    !local
    real*8                             :: timeScale, tempScale, BRDsq, temp


    temp = insTemp

    if (temp == 0.d0 ) temp = 1.d0 

    tempScale = bathTemp/temp 
    timeScale = dt/bathCoup  


    BRDsq = ( (tempScale - 1.d0) * timeScale ) + 1.d0 

    BRDscale = sqrt(BRDsq)  

end subroutine calculate_BRDscale


subroutine velocity_verlet(pl, BWRadius, BWVel, BWEforce, BWVforce, time, dt, rhoElReal, & 
                elHam, FWRadius, FWVel, FWEforce, FWVforce, FWTemp, V_Energy, K_Energy )
    implicit none

    ! args
    integer,               intent(in)  :: pl 
    real*8,                intent(in), dimension(nr, nc) :: BWRadius, BWVel, BWEforce, BWVforce
    real*8,                intent(in)  :: time, dt
    real*8,                intent(in)  :: rhoElReal(d_el, d_el)
    real*8,                intent(in)  :: elHam(d_el, d_el)
    real*8,  allocatable,  intent(out) :: FWRadius(:,:), FWVel(:,:), FWEforce(:,:), FWVforce(:,:)
    real*8,                intent(out) :: FWTemp, V_Energy, K_Energy
    
    ! local  
    real*8, allocatable :: BWForce(:,:), FWForce(:,:)
    real*8, allocatable :: BWAcc(:,:), FWAcc(:,:)
    integer             :: k, i, j 
    real*8              :: dt_ps
    real*8              :: BRDscale
    !BWRadius -> BackWardRadius (raio anterior).    FWRadius -> ForWardRadius (raio posterior).
    
    allocate( FWRadius(nr, nc)    , source = 0.d0 ) 
    allocate( FWVel(nr, nc)       , source = 0.d0 )
    allocate( BWForce(nr, nc)     , source = 0.d0 )
    allocate( FWForce(nr, nc)     , source = 0.d0 )
    allocate( BWAcc(nr, nc)       , source = 0.d0 )
    allocate( FWAcc(nr, nc)       , source = 0.d0 ) 


    !ALOCO FWEforce na subrotina eletric_force
    dt_ps = dt * 1.d-12   

    !====== FORCA TOTAL EM t   
    BWForce =  BWEforce + BWVforce
    
    !======= ACELERACAO COM A FORCA EM t
    BWAcc = BWForce / site%mass
    
    !====== CALCULO DA ENERGIA CINETICA E POTENCIAL COM A FORCA EM t 
    call kinect_energy(BWVel, K_Energy) 
    call spring_energy(BWRadius, V_Energy)

    
    !===== CALCULO A TEMPERATURA DO SISTEMA COM A NOVA ENERGIA CINETICA
    call calculate_temperature(K_energy, FWTemp) 

    !===== BERENDSEN SCALE 
    call calculate_BRDscale(dt_ps, FWTemp, BRDscale)

    if (BRD_on == 1 ) BRDscale = 1.d0
    
    !======= CALCULO DOS NOVOS RAIOS EM t
    FWRadius = BWRadius + (BWVel * dt_ps) + (HALF * BWAcc * dt_ps * dt_ps)

    !======= ATUALIZO O RAIO E OMEGA 
    site%radius = FWRadius
    site%omega = ( ( 2.d0 * hbar) / ( me  * site%radius * site%radius ) ) * hz_to_thz

    !======= CALCULO DAS NOVAS FORCAS EM t + delta t
    call eletric_force(pl, rhoElReal, time, elHam, FWEforce)
    call spring_force(FWVforce) 
    
    !======= TOTAL FORWARD FORCE 
    FWForce = FWEforce + FWVforce

    !======== NOVA ACELERACAO
    FWAcc = FWForce / site%mass
    
    !======== NOVA VELOCIDADE RADIAL
    FWVel = BWVel + HALF * dt_ps * (BWAcc + FWAcc) 
    
    !======== BERENDSEN THERMOSTAT 
    FWVel = FWVel * BRDscale
 
    !======== ATUALIZO A VELOCIDADE RADIAL 
    site%vel = FWVel


    deallocate(BWForce, FWForce, BWAcc, FWAcc) 
end subroutine velocity_verlet

end module verlet_m
