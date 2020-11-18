program dynamicsOfOQS
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m
use system_hamiltonian_m
use rdftensor_m
use functions_m
use ordinary_equation_solution
use time_evolution_m
use verlet_m 

implicit none

    real*8  :: InitSitesRadius(nr, nc), InitSitesVel(nr, nc) 

    if ( therm_on == .true. ) then
        ElCoup = .false.  
        InitSitesRadius = raioZero
        InitSitesVel = velZero
    else
        ElCoup = .true. 
        open ( unit=301, file="therm_radius.dat", status="old", access="sequential" )
        read(301, *) ( ( InitSitesRadius(1, i) ), i = 1, nsites )
        close(301)

        open ( unit=302, file="therm_vel.dat"  ,  status="old", access="sequential" )
        read(302, *) ( ( InitSitesVel(1, i) ), i = 1, nsites )
        close(302)

        InitSitesRadius = InitSitesRadius * 1.d-9 
    endif




    call System_Dynamics(InitSitesRadius, InitSitesVel, nm_divisoes)





end program DynamicsOfOQS
