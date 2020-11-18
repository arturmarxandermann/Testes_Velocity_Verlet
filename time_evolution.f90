module time_evolution_m
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
use gnugraphs_m 

private

public :: PropOfWavePacket, f0, System_Dynamics


contains


!se as matrizes são definidas no parameters.f90, não preciso criar elas como intent(in)
!pois todas as subrotinas utilizam parameters.f90


subroutine System_Dynamics(InitSitesRadius, InitSitesVel, nm_divisoes)
    implicit none

    ! args 
    real*8,  intent(in) :: InitSitesRadius(nr, nc), InitSitesVel(nr, nc) 
    integer, intent(in) :: nm_divisoes

    ! local 
    integer                     :: pl
    real*8                      :: ti, tf
    real*8                      :: EnergiaEl                       !energia do eletron
    real*8                      :: BWTemp, FWTemp                  !temperatura instantanea do sistema
    real*8                      :: K_energ_init
    real*8                      :: BWRadius(nr, nc), BWVel(nr, nc) !Backward variables
    real*8                      :: V_Energy, K_Energy              !energia da mola e energia cinetica
    real*8,     allocatable     :: BWEforce(:,:)                   !aloco na subrotina eletric_force
    real*8,     allocatable     :: BWVforce(:,:)                   !aloco na subrotina spring_force
    real*8,     allocatable     :: hMtx(:,:)                       !hamitoniano calculado em t + delta t
    real*8,     allocatable     :: FWRadius(:,:), FWVel(:,:)
    real*8,     allocatable     :: FWEforce(:,:), FWVforce(:,:)    !Forward variables
    real*8,     allocatable     :: rhoReal(:,:)                    !parte real da matriz densidade
    complex*16, allocatable     :: rhoSites_in(:,:)                !matriz densidade inicial
    complex*16, allocatable     :: rhoSitesPolar_in(:,:)                !matriz densidade inicial
    complex*16, allocatable     :: rhoSites(:,:)                   !matriz densidade calculada em t + delta t
    real*8                      :: pop_el(nr, nc)          

     
    integer :: nn, ll, points 
    real*8 :: EQForce 
    complex*16, allocatable :: TMtx(:,:), TMtxDg(:,:)

    allocate( rhoSites_in(d_el, d_el), source = (0.d0, 0.d0))
    allocate( rhoSitesPolar_in(d_el, d_el), source = (0.d0, 0.d0))
    allocate( rhoReal(d_el, d_el ),    source = 0.d0        ) 
    allocate( TMtxDg(d_el, d_el),      source = (0.d0, 0.d0)) 
    allocate( basis(nsites, nsites)                         ) 
    do j = 1, nsites
     do i = 1, nsites
        allocate ( basis(i, j)%hMtx(ns_el, ns_el),   source = 0.d0         )
        allocate ( basis(i, j)%DerMtx(ns_el, ns_el), source = 0.d0         )
        allocate ( basis(i, j)%TMtx(ns_el, ns_el),   source = (0.d0, 0.d0) ) 
     enddo
    enddo

  
    call Basis_Builder_Blocks
    call define_sitios(InitSitesRadius, InitSitesVel, site, site_point) 


    call build_TMatrix(TMtx)
    

    TMtxDg = transpose( conjg(TMtx) )


    call build_EnergyVectors(Qn, Qn_erg)  



    !============= CONDICOES INICIAIS NA BASE(x,y) DOS SITIOS ==============
    rhoSitesPolar_in(initState, initState) = 1.d0 + 0.d0 * zi  !matriz densidade
    BWRadius = site(:, :)%radius                          !raio inicial
    BWVel = site(:, :)%vel                                !velocidade inicial
    rhoReal = real(rhoSites_in)                           !parte real da matriz densidade
    call build_hamiltonian(hMtx)
    call eletric_force(pl, rhoReal, 0.d0, hMtx, BWEforce) !Backward Eforce
    call spring_force(BWVforce)                           !Backward SpringForce
    call kinect_energy(BWVel, K_energ_init)               !Energia cinetica inicial
    !================================================================
    

    !============ CONDICOES INICIAIS NA BASE(r, theta) DOS SITIOS
    rhoSites_in = matmul(TMtx, matmul(rhoSitesPolar_in, TMtxDg) ) 

    !======= SUBROTINA QUE ABRE OS ARQUIVOS DE ESCRITA =============
    call open_write_files
    !==============================================================
    
    ti = 0.d0

    time = 0.d0 
   
    do pl = 1, nm_divisoes !NÃO POSSO USAR i NESSSE LOOP PQ ELE JÁ É UTILIZADO EM TODOS OS OUTROS LOOPS NA SUBROTINA PropOfWavePacket01..
      if (pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
        print*, "==== PASSO",pl,"===="
      endif


      call rhoMtx_to_pop(rhoSites_in, pop_el)
      
      call GnuFiles(pop_el, pl) 

      !============= DEFINE O TEMPO FINAL E O DELTA_t ===========
      tf = float(pl)*(tmax)/float(nm_divisoes)
      dt = tf - ti 
      !==========================================================
  
      time = time + dt    
    
    
      !=============== CALCULO DO HAMILTONIANO  ================
      call build_hamiltonian(hMtx)
      !==========================================================
     
      write(80, "(60F20.8)") ti, hMtx(1, 1) !, hMtx(4, 4), hMtx(4, 4) - hMtx(1, 1)       
    
      !============== EVOLUCAO PARA O ELETRON DE ti ATÉ tf =======
      call PropOfWavePacket(pl, ti, tf, hMtx, rhoSites_in, rhoSites, EnergiaEl) !propago função de onda
      !===========================================================
      
      write(201, "(60F20.8)") ti,  ( ( real(rhoSites_in(i, i)) ), i = 1, d_el ) 
      
      !============= ATUALIZO AS CONDICOES INICIAIS QUANTICAS ==============
      rhoSites_in = rhoSites !atualizo cond. inicial para o eletron
      rhoReal = real(rhoSites_in) 
      !===========================================================
    

      !============ CALCULO A EVOLUCAO CLASSICA COM O ALGORITMO VELOCITY VERLET ===========================
    
      call velocity_verlet(pl, BWRadius, BWVel, BWEforce, BWVforce, ti, dt, rhoReal, &
                           hMtx, FWRadius, FWVel, FWEforce, FWVforce, FWTemp, V_Energy, K_Energy) 
                  !calculamos a energia cinetica e da mola no algoritmo de verlet em um tempo ti para ser igual ao resto do programa
      !====================================================================================================

      !============ ESCREVE AS ENERGIAS CLASSICAS E QUANTICAS ============================================
      write(84, "(60F20.8)")  ti, K_Energy
      write(85, "(60F20.8)")  ti, V_Energy
      write(86, "(60F20.8)")  ti, energiazeroel 
      write(87, "(60F20.8)")  ti, EnergiaEl
      write(103, "(60F20.8)") ti, EnergiaEl + V_Energy + K_Energy
      write(105, "(60F20.8)") ti, FWTemp       
      !===================================================================================================
 
    
      !=========== ATUALIZA A PARTE CLASSICA =============================================================
      BWRadius(:, :) = FWRadius(:, :) !aloco em velocity_verlet
      BWVel(:, :)    = FWVel(:, :)    !aloco em velocity_verlet
      BWEforce(:, :) = FWEforce(:, :) !aloco em eletric_force
      BWVforce(:, :) = FWVForce(:, :) !aloco em spring_force
      deallocate(FWEforce, FWVel, FWRadius, FWVforce) 
      !===================================================================================================

      
      !========== ESCREVE OS RAIOS E AS FORÇAS ===========================================================
      write(100,  "(60F20.8)") tf, ( ( site(1, i)%radius*1.d9 ), i = 1, nsites ) 
      write(101,  "(60F20.8)") tf, ( ( BWVforce(1, i)*1.d10   ), i = 1, nsites ) 
      write(102,  "(60F20.8)") tf, ( ( BWEforce(1, i)*1.d10   ), i = 1, nsites ) 
      !===================================================================================================
      
      !========= ATUALIZO O TEMPO ===========
      ti = tf
      !=====================================
    enddo

    if ( therm_on == .true. ) then 

        print*, "Termalização acabou, escrevendo os arquivos de saída"

        open( unit = 301, file = "therm_radius.dat", status = "replace" )
        write(301, "(60F20.8)") ( ( site_point(i)%np%radius*1.d9 ),  i = 1, nsites ) 
        close(301) 


        open( unit = 302, file = "therm_vel.dat", status = "replace" )
        write(302, "(60F20.8)") ( ( site_point(i)%np%vel ),  i = 1, nsites ) 
        close(302) 

    endif 

    !========= FECHO OS ARQUIVOS DE ESCRITA =====
    call close_write_files
    !============================================
    
    
    
    DEALLOCATE(rhoSites_in, rhoSitesPolar_in)
end subroutine System_Dynamics



subroutine PropOfWavePacket(step, ti, tf, hMtx, rhoSites_in, rhoSites, ParticleEnergy)
    implicit none
   
    ! args 
    integer,                 intent(in)  :: step 
    real*8,                  intent(in)  :: ti 
    real*8,                  intent(in)  :: tf 
    real*8,                  intent(in)  :: hMtx(d_el, d_el) 
    complex*16,              intent(in)  :: rhoSites_in(d_el, d_el)
    complex*16, allocatable, intent(out) :: rhoSites(:,:)
    real*8,                  intent(out) :: ParticleEnergy
    
    ! local 
    real*8,     allocatable  :: energias(:), y(:), yp(:), work(:)
    REAL*8,     allocatable  :: phi(:,:), phi_transpose(:,:)
    COMPLEX*16, allocatable  :: rhoHam_in(:,:), rhoHam(:,:)
    INTEGER                  :: number_file_site, number_file_ham
    real*8                   :: abserr
    integer                  :: iflag
    integer                  :: flag
    integer                  :: iwork(5)
    real*8                   :: relerr
    real*8                   :: t !esse t é a variável independente do ODE
    real*8                   :: tout 
    real*8                   :: energ(d_el) 
    
    
    
     allocate( rhoSites(d_el, d_el)           , source = (0.d0, 0.d0) )
     allocate( y(neqn_el)                     , source = 0.d0         )
     allocate( yp(neqn_el)                    , source = 0.d0         )
     allocate( work(100+21*neqn_el))
     allocate( rhoHam_in(d_el, d_el)          , source = (0.d0, 0.d0) ) !populaca inicial na base da hMtx
     allocate( rhoHam(d_el, d_el)             , source = (0.d0, 0.d0) )
    
    
     number_file_ham = 14
     number_file_site = 15
     abserr = 1.d-9
     relerr = 1.d-9
     flag = +1 !FLAG -1 NAO FUNCIONA
     iflag = 1 
     t = ti !construo o tempo inicial de uma forma que o ODE nao atualiza!
     tout = tf !TOUT É O TEMPO DE SAÍDA, OU SEJA, DELTA T = tout(2) - tout(1)
    
    
    
    !==== CALCULA OS AUTOESTADOS E AUTOVETORES =======
    call calculate_eigenvectors(step, hMtx, energias, phi, phi_transpose, frequency_matrix) 
    !===============================================================================================
    
    
    !==== PRINTA A MATRIZ DENSIDADE, AUTOVETORES E HAMILTONIANO =====
    call print_matrices(step, rhoSites_in, phi, hMtx) 
    !===============================================================================================
    

    !==== DEFINE A ENERGIA ZERO ===========
    if (t == 0.d0 ) then
      energiazeroel = hMtx(initState, initState)
    endif
    !=====================================
    
    
    !=========== CONDICOES INICIAIS NA BASE EXC ======================
    call rhosite_TO_rhoham(phi, phi_transpose, rhoSites_in, rhoHam_in)
    !=================================================================
    
    
    
    !========== CALCULO DA ENERGIA DA PARTICULA ====================
    call particle_energy(hMtx, rhoSites_in, ParticleEnergy)
    !==============================================================
    
    
    !======== TRANSFORMO AS CONDICOES INICIAIS PARA OS VETORES LINHAS DO ODE ========
    call monta_y(rhoHam_in, y) !condição inicial para resolver as ODE'S
    !================================================================================
    
    
    !====== CRIO A MATRIZ RW QUE É UTILIZADA PARA CALCULAR O ODE =================
    call createRwMatrix(frequency_matrix, RWMatrix) 
    !=============================================================================
    
    
    
    !y esta na forma (para d=2):  y(1) = rho(1, 1)  y(3) = real(rho(2, 1))
    !                             y(2) = rho(2, 2)  y(4) = aimag(rho(2, 1))
    !----------------------------------------------------------------------------------
      
    
    !====== ESCREVE AS POPULACOES =====================
    call printa_resultado(number_file_ham,  t,  rhoHam_in)
    call printa_resultado(number_file_site, t, rhoSites_in)
    !==================================================
    
    
    
    !======== SUBROTINA QUE CALCULA AS EQUACOES DIFERENCIAIS DE T ATÉ TOUT ========
    call ode ( f0, neqn_el, y, t, tout, relerr, abserr, iflag, work, iwork )
    
    if ( iflag /= 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PropOfWavePacket - Fatal error!'
      write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
      stop
    end if
    !==============================================================================
    
    
    
    !========== TRANSFORMO O RESULTADO DO VETOR LINHA Y DO ODE PARA AS MATRIZES RHO ===========
    call monta_rho(y, rhoHam)
    !==========================================================================================
    
    
    !============== ESCREVO O OPERADOR DENSIDADE NA BASE DO SITIO =============================
    call rhoham_TO_rhosite(phi, phi_transpose, rhoHam, rhoSites)
    !==========================================================================================
    
     
        
    
    DEALLOCATE(phi, phi_transpose, rhoHam_in, rhoHam, RWMatrix)
    13 format(3es14.3E3)
end subroutine PropOfWavePacket



subroutine f0 ( t, y, yp )
!não definimos as condicoes iniciais aqui, apenas a funcao que vamos calcular a derivada. A condição inicial vem antes do call
    implicit none

    real*8 ::  t !tempo que vai entrar na equação diferencial
    real*8 ::  y(neqn_el)  !funcao que estamos querendo, neste caso x
    real*8 ::  yp(neqn_el) !derivada que estamos querendo calcular, por exemplo dx/dt = -x => sol: x = exp(-t)
    
    
    call gemv(RWMatrix, y, yp)  !gemv calcula o produto da matriz RWMatrix com
                                !o vetor y resultando em yp => muito mais rapido
    
    13 format(3es14.3E3)
end subroutine f0



end module time_evolution_m
