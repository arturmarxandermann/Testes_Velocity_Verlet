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
!se el_or_hl = 1 a ham. output é para eletron
    implicit none

    !args
    real*8, INTENT(OUT), allocatable :: hamiltonian(:,:)

    !local
    integer             :: nmstates, c1, r1, s1, s2, c2, r2
    real*8, allocatable :: S(:,:)
    real*8, allocatable :: coup_temp(:,:)
    real*8,             :: site_type, omegas(nm_rows, nm_columns)
    real*8,             :: Qnumbers_energy(6)
    INTEGER,            :: Qnumbers_row(6)
    
    !criando o vetor associado aos valores dos numeros quanticos dos estados do osciladores harmonicos
    Qnumbers_row = [ 00, 01, 10, 02, 11, 20 ]
    
    !criando o vetor com as energias associado a cada numero quantico
    Qnumbers_energy = [ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 ]
    
    nmstates       = nmstates_el
    site_type(:,:) = sites_array(:,:)%site_type_el
    omegas(:,:)    = sites_array(:,:)%omega
    
    allocate(hamiltonian(dim_el , dim_el) , source = 0.d0 )
    allocate(coup_temp(nmstates, nmstates), source = 0.d0 )
    allocate(S(nmstates, nmstates)        , source = 0.d0 ) 

    !---- PARTE DIAGONAL -----
    k = 1
    do j = 1, nm_columns
      do i = 1, nm_rows
        do s1 = 1, nmstates
           hamiltonian(k, k) = site_type(i, j) + HB_ev_ps * omegas(i, j) *  Qnumbers_energy(s1) 
           k = k + 1
         enddo
       enddo
    enddo
    !-------------------------
    
    
    !---- PARTE NÂO DIAGONAL ---------
    l = 0
    k = 0
    do c1 = 1, nm_columns
      do r1 = 1, nm_rows
        k = 0
        do c2 = 1, nm_columns
          do r2 = 1, nm_rows 
            if ( k > l ) then 
             call NewOverlap(sites_array, nmstates, r2, c2, r1, c1, Qnumbers_row, S)
              do s2 = 1, nmstates
               do s1 = 1, nmstates
                 
                  coup_temp(s1, s2) = S(s1, s2) * sites_array(r2, c2)%t
                 
                enddo
              enddo
              hamiltonian(k*nmstates+1:k*nmstates+nmstates, l*nmstates+1:l*nmstates+nmstates) = &
              RESHAPE(coup_temp, (/nmstates, nmstates/) ) !triangulo de baixo

              hamiltonian(l*nmstates+1:l*nmstates+nmstates, k*nmstates+1:k*nmstates+nmstates) = &
              RESHAPE(transpose(coup_temp), (/nmstates, nmstates/) ) !triangulo de cima
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
    real*8              , intent(in)  :: hamiltoniana(dim_el,dim_el)
    real*8 , allocatable, intent(out) :: energias(:)
    real*8 , allocatable, intent(out) :: phi, phi_transpose(:,:), omega_matrix(:,:)

    ! local 
    real*8, allocatable ::  hamiltoniana_diag(:,:)
    character(len=1)    :: jobz, uplo
    integer             :: info, is, js 

    allocate(phi(dim_el, dim_el)               , source = 0.d0)
    allocate(phi_transpose(dim_el, dim_el)     , source = 0.d0)
    allocate(omega_matrix(dim_el, dim_el)      , source = 0.d0)
    allocate(energias(dim_el)                  , source = 0.d0)
    allocate(hamiltoniana_diag(dim_el, dim_el) , source = 0.d0)
    
    hamiltoniana_diag = hamiltoniana
    
    call syevd(hamiltoniana_diag, energias, "V", "L", 0)
    
    !If info = 0, contains the n orthonormal eigenvectors, stored by columns.
    !(The i-th column corresponds to the ith eigenvalue.)
    !Como o iézimo autovetor é disposto na iézima coluna, C_n^N = <n|N> é o coeficiente da expansão do autovetor
    
    if ( info /= 0 ) write(*,*) info, "Error to found eigenvalues"
    
    if ( pl == 1 .OR. pl == nm_divisoes/2 .OR. pl == nm_divisoes ) then
      do i = 1, dim_el
      print*, "ENERGIA", i, energias(i)
      enddo
    endif
    
    !========== CONSTRUCAO DA MATRIZ DE AUTOVETORES ========
    do n = 1, dim_el
      phi(n,:) = hamiltoniana_diag(n,:)
    enddo
    !=======================================================
    
    phi_transpose = transpose(phi)
    
    !============== CONSTRUÇÃO DA MATRIZ DE FREQUÊNCIAS ========================
    do j = 1, dim_el
      do i = 1, dim_el
        omega_matrix(i, j) = (energias(i) - energias(j)) / HB_ev_ps !omegas em 1/s
      enddo
    enddo
    !============================================================================
    
    13 format(3es14.5E3)
    
    deallocate(hamiltoniana_diag) 

end subroutine calculate_eigenvectors


subroutine build_overlap_matrix(nmstates, dim_el, ovlpmatrix)
    implicit none

    !args
    integer            , intent(in)  :: nmstates, dim_el
    real*8, allocatable, intent(out) :: ovlpmatrix(:,:)
   
    !local 
    real*8, allocatable :: S(:,:)
    integer             :: c1, c2, r1, r2
    INTEGER             :: Qnumbers_row(6)

    allocate(ovlpmatrix(dim_el, dim_el), source = 0.d0)
    
    Qnumbers_row = (/ 00, 01, 10, 02, 11, 20 /)

    l = 0
    k = 0
    do c1 = 1, nm_columns
      do r1 = 1, nm_rows
          k = 0
        do c2 = 1, nm_columns
          do r2 = 1, nm_rows
            if ( k > l ) then
            call Overlap(sites_array, nmstates, r2, c2, r1, c1, Qnumbers_row, S)

            ovlpmatrix(k*nmstates+1:k*nmstates+nmstates, l*nmstates+1:l*nmstates+nmstates) =  &
            RESHAPE(S, (/nmstates, nmstates/) )

            ovlpmatrix(l*nmstates+1:l*nmstates+nmstates, k*nmstates+1:k*nmstates+nmstates) =  &
            RESHAPE(transpose(S), (/nmstates, nmstates/) )
    
            endif
            k = k + 1
          enddo
        enddo
          l = l + 1
      enddo
    enddo
    
    do j = 1, dim_el
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
    
    allocate(dermtx(dim_el, dim_el), source = 0.d0 )

    l = 0
    k = 0
    do c1 = 1, nm_columns
      do r1 = 1, nm_rows
          k = 0
        do c2 = 1, nm_columns
          do r2 = 1, nm_rows
            if ( k > l ) then

            call DerivativeOverlap(sites_array, nmstates_el, r2, c2, r1, c1, derivativeMtx)

            dermtx(k*nmstates_el+1:k*nmstates_el+nmstates_el, l*nmstates_el+1:l*nmstates_el+nmstates_el) =  &
            RESHAPE(derivativeMtx, (/nmstates_el, nmstates_el/) ) !triangulo de baixo
    
            dermtx(l*nmstates_el+1:l*nmstates_el+nmstates_el, k*nmstates_el+1:k*nmstates_el+nmstates_el) =  &
            RESHAPE(transpose(derivativeMtx), (/nmstates_el, nmstates_el/) ) !triangulo de cima
                                    
            endif 
            k = k + 1
          enddo
        enddo
          l = l + 1
      enddo
    enddo

end subroutine calculate_derivative_matrix


subroutine printaletters3(nomearquivo, nmfile, matrizsize, nmlines, nmcol, matriz)
    implicit none

    !args
    character(len=6), intent(in) :: nomearquivo
    integer,          intent(in) :: nmfile, matrizsize, nmlines, nmcol
    real*8,           intent(in) :: matriz(matrizsize,matrizsize)

    open(file = nomearquivo, status = "replace", unit = nmfile)
    do i = 1, nmlines
       do j = 1, nmcol
          write(nmfile, 13, advance = "no") matriz(i, j)
       enddo
       write(nmfile, 13, advance = "yes")
    enddo
    
    close(nmfile)
    
    13 format (8F6.2)
end subroutine printaletters3


end module system_hamiltonian_m
