module system_hamiltonian_m
!modulo com subrotinas para criar a hamiltoniana do sistema
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use overlap_m

contains

subroutine print_mat3(aa, nn, mm)
    implicit none

    ! args
    integer, intent(in) :: nn, mm
    real*8 , intent(in) :: aa(nn, mm)

    !local 
    integer :: i, j

    do, i=1,mm
        write(*,'(100g12.4)') ( aa(i,j), j=1,nn )
    enddo
end subroutine print_mat3


subroutine build_hamiltonian(hamiltonian)
    implicit none

    !args
    real*8, INTENT(OUT), allocatable :: hamiltonian(:,:)

    !local
    integer             :: ns, c1, r1, s1, s2, c2, r2
    real*8, allocatable :: S(:,:)
    real*8, allocatable :: coup_temp(:,:)
    real*8              :: site_type(nr, nc), omegas(nr, nc)
    real*8              :: Qn_erg(6)
    INTEGER             :: Qn(6)
    
    !criando o vetor associado aos valores dos numeros quanticos dos estados do osciladores harmonicos
    Qn = [ 00, 01, 10, 02, 11, 20 ]
    
    !criando o vetor com as energias associado a cada numero quantico
    Qn_erg = [ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 ]
    
    ns       = ns_el
    site_type(:,:) = site(:,:)%V0
    omegas(:,:)    = site(:,:)%omega
    
    allocate(hamiltonian(d_el , d_el) , source = 0.d0 )
    allocate(coup_temp(ns, ns), source = 0.d0 )
    allocate(S(ns, ns)        , source = 0.d0 ) 

    !---- PARTE DIAGONAL -----
    k = 1
    do j = 1, nc
      do i = 1, nr
        do s1 = 1, ns
           hamiltonian(k, k) = site_type(i, j) + HB_ev_ps * omegas(i, j) *  Qn_erg(s1) 
           k = k + 1
         enddo
       enddo
    enddo
    !-------------------------
    
    
    !---- PARTE NÂO DIAGONAL ---------
    l = 0
    k = 0
    do c1 = 1, nc
      do r1 = 1, nr
        k = 0
        do c2 = 1, nc
          do r2 = 1, nr 
            if ( k > l ) then 
             call NewOverlap(site, ns, r2, c2, r1, c1, Qn, S)
              do s2 = 1, ns
               do s1 = 1, ns
                 
                  coup_temp(s1, s2) = S(s1, s2) * site(r2, c2)%t
                 
                enddo
              enddo
              hamiltonian(k*ns+1:k*ns+ns, l*ns+1:l*ns+ns) = &
              RESHAPE(coup_temp, [ns, ns] ) !triangulo de baixo

              hamiltonian(l*ns+1:l*ns+ns, k*ns+1:k*ns+ns) = &
              RESHAPE(transpose(coup_temp), [ns, ns] ) !triangulo de cima
            endif
          k = k + 1
          enddo
        enddo
        l = l + 1
      enddo
    enddo
    !------------------------------
    
DEALLOCATE(S)
end subroutine build_hamiltonian



subroutine calculate_eigenvectors(pl, hamiltoniana, energias, phi, phi_transpose, omega_matrix)
!calcula os autovetores (phi) e autovalores (energias) da hamiltoniana. omega_matriz é a matriz de frequencias angulares
    implicit none

    ! args
    integer             , intent(in)  :: pl
    real*8              , intent(in)  :: hamiltoniana(d_el,d_el)
    real*8 , allocatable, intent(out) :: energias(:)
    real*8 , allocatable, intent(out) :: phi(:,:), phi_transpose(:,:), omega_matrix(:,:)

    ! local 
    real*8, allocatable ::  H_temp(:,:)
    character(len=1)    :: jobz, uplo
    integer             :: info, is, js 

    allocate(phi(d_el, d_el)               , source = 0.d0)
    allocate(phi_transpose(d_el, d_el)     , source = 0.d0)
    allocate(omega_matrix(d_el, d_el)      , source = 0.d0)
    allocate(energias(d_el)                  , source = 0.d0)
    allocate(H_temp(d_el, d_el) , source = 0.d0)
   
    info = 0 
    H_temp = hamiltoniana
    
    call syevd(H_temp, energias, "V", "L", info)
    
    !If info = 0, contains the n orthonormal eigenvectors, stored by columns.
    !(The i-th column corresponds to the ith eigenvalue.)
    !Como o iézimo autovetor é disposto na iézima coluna, C_n^N = <n|N> é o coeficiente da expansão do autovetor
    
    if ( info /= 0 ) write(*,*) info, "Error to found eigenvalues"
    
    if ( pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
      do i = 1, d_el
      print*, "ENERGIA", i, energias(i)
      enddo
    endif
    
    !========== CONSTRUCAO DA MATRIZ DE AUTOVETORES ========
    do n = 1, d_el
      phi(n,:) = H_temp(n,:)
    enddo
    !=======================================================
    
    phi_transpose = transpose(phi)
    
    !============== CONSTRUÇÃO DA MATRIZ DE FREQUÊNCIAS ========================
    do j = 1, d_el
      do i = 1, d_el
        omega_matrix(i, j) = (energias(i) - energias(j)) / HB_ev_ps !omegas em 1/s
      enddo
    enddo
    !============================================================================
    
    13 format(3es14.5E3)
    
    deallocate(H_temp) 

end subroutine calculate_eigenvectors


subroutine build_overlap_matrix(ns, d_el, ovlpmatrix)
    implicit none

    !args
    integer            , intent(in)  :: ns, d_el
    real*8, allocatable, intent(out) :: ovlpmatrix(:,:)
   
    !local 
    real*8, allocatable :: S(:,:)
    integer             :: c1, c2, r1, r2
    INTEGER             :: Qn(6)

    allocate(ovlpmatrix(d_el, d_el), source = 0.d0)
    
    Qn = [ 00, 01, 10, 02, 11, 20 ]

    l = 0
    k = 0
    do c1 = 1, nc
      do r1 = 1, nr
          k = 0
        do c2 = 1, nc
          do r2 = 1, nr
            if ( k > l ) then
            call Overlap(site, ns, r2, c2, r1, c1, Qn, S)

            ovlpmatrix(k*ns+1:k*ns+ns, l*ns+1:l*ns+ns) =  &
            RESHAPE(S, [ns, ns] )

            ovlpmatrix(l*ns+1:l*ns+ns, k*ns+1:k*ns+ns) =  &
            RESHAPE(transpose(S), [ns, ns] )
    
            endif
            k = k + 1
          enddo
        enddo
          l = l + 1
      enddo
    enddo
    
    do j = 1, d_el
        ovlpmatrix(j,j) = 1.d0 
    enddo 
    
end subroutine build_overlap_matrix


subroutine calculate_derivative_matrix(dermtx)
    implicit none

    !args
    real*8, allocatable, intent(out) :: dermtx(:,:)

    !local 
    real*8, allocatable :: derivativeMtx(:, :)
    integer             :: c1, r1, c2, r2
    
    allocate(dermtx(d_el, d_el), source = 0.d0 )

    l = 0
    k = 0
    do c1 = 1, nc
      do r1 = 1, nr
          k = 0
        do c2 = 1, nc
          do r2 = 1, nr
            if ( k > l ) then

            call DerivativeOverlap(site, ns_el, r2, c2, r1, c1, derivativeMtx)

            dermtx(k*ns_el+1:k*ns_el+ns_el, l*ns_el+1:l*ns_el+ns_el) =  &
            RESHAPE(derivativeMtx, (/ns_el, ns_el/) ) !triangulo de baixo
    
            dermtx(l*ns_el+1:l*ns_el+ns_el, k*ns_el+1:k*ns_el+ns_el) =  &
            RESHAPE(transpose(derivativeMtx), (/ns_el, ns_el/) ) !triangulo de cima
                                    
            endif 
            k = k + 1
          enddo
        enddo
          l = l + 1
      enddo
    enddo

end subroutine calculate_derivative_matrix


end module system_hamiltonian_m
