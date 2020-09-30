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


subroutine build_hamiltonian(hMtx)
    implicit none

    !args
    real*8, INTENT(OUT), allocatable :: hMtx(:,:)

    !local
    integer             :: c1, r1, s1, s2, c2, r2

    
    integer :: row_min, row_max, col_min, col_max
    
    
    allocate(hMtx(d_el , d_el)        , source = 0.d0 )


    do j = 1, nsites-1
      do i = j+1, nsites
        row_min = basis(i, j)%rmin
        row_max = basis(i, j)%rmax
        col_min = basis(i, j)%cmin
        col_max = basis(i, j)%cmax
          
        hMtx( row_min : row_max, col_min : col_max ) = basis(i, j)%hMtx  
      enddo
    enddo
  
    hMtx = hMtx + transpose(hMtx) 

    do i = 1, nsites
      row_min = basis(i, i)%rmin
      row_max = basis(i, i)%rmax
      col_min = basis(i, i)%cmin
      col_max = basis(i, i)%cmax
      
      hMtx( row_min : row_max, col_min : col_max ) = basis(i, i)%hMtx
    enddo



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
    allocate(energias(d_el)                , source = 0.d0)
    allocate(H_temp(d_el, d_el)            , source = 0.d0)
   
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
    do j = 1, d_el-1
      do i = j + 1, d_el
        omega_matrix(i, j) = (energias(i) - energias(j)) / HB_ev_ps 
      enddo
    enddo

    omega_matrix = omega_matrix + transpose(omega_matrix) 
    !============================================================================
    
    13 format(3es14.5E3)
    
    deallocate(H_temp) 

end subroutine calculate_eigenvectors


subroutine build_derivative_matrix(DerMtx)
    implicit none

    !args
    real*8, allocatable, intent(out) :: DerMtx(:,:)

    !local 
    integer :: row_min, row_max, col_min, col_max
    
    allocate( DerMtx(d_el, d_el), source = 0.d0 )


    do j = 1, nsites-1
      do i = j + 1, nsites

        row_min = basis(i, j)%rmin
        row_max = basis(i, j)%rmax
        col_min = basis(i, j)%cmin
        col_max = basis(i, j)%cmax
        
        DerMtx( row_min : row_max, col_min : col_max ) = basis(i, j)%DerMtx
      enddo
    enddo
    
    DerMtx = DerMtx + transpose(DerMtx) 

end subroutine build_derivative_matrix


end module system_hamiltonian_m
