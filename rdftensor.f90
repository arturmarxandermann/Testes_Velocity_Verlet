module rdftensor_m
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use functions_m
use omp_lib 
contains


subroutine createRwMatrix(omega_matrix, Rw_matrix)
    implicit none
!subrotina que cria a matriz Rw que relaciona os vetores rhodot e rho
!na forma rho_dot = matmul(Rw, rho).
!Rw é composta por 9 blocos de matrizes, que são numerados na seguinte ordem
! Rw = ( ( 1 ) ( 2 )  ( 3 )  )
!      ( ( 4 ) ( 5 )  ( 6 )  )
!      (  (7 ) ( 8 )  ( 9 )  )
!Como o bloco 1, 4 e 7 relaciona elementos rhoreal_aa, rhoreal_ab, rhoimag_ab
!com rhoreal_aa
!seu tamanho é dim_el. Os outros blocos tem tamanho non_diag pois relacionam os
!elementos
!de antes com rhoreal_ab, rhoimag_ab

    ! args
    REAL*8,              INTENT(IN)    :: omega_matrix(dim_el, dim_el) 
    REAL*8, allocatable, INTENT(OUT)   :: Rw_matrix(:,:)
    
    !local 
    INTEGER                        :: non_diag, rhosquared, aux1
    TYPE(tensor_product_matrices)  :: Rw_blocks(9)
    !TRmtc são matrizes utilziadas no produto tensorial para montar Rw
    
    non_diag   = (dim_el * ( dim_el -1 ) ) / 2
    rhosquared = dim_el * dim_el
    aux1       = dim_el + non_diag
    
    
    ALLOCATE( Rw_matrix(rhosquared, rhosquared)        , source = 0.d0 )
    ALLOCATE( Rw_blocks(1)%elements(dim_el, dim_el)    , source = 0.d0 )
    ALLOCATE( Rw_blocks(2)%elements(dim_el, non_diag)  , source = 0.d0 )
    ALLOCATE( Rw_blocks(3)%elements(dim_el, non_diag)  , source = 0.d0 )
    ALLOCATE( Rw_blocks(4)%elements(non_diag, dim_el)  , source = 0.d0 )
    ALLOCATE( Rw_blocks(5)%elements(non_diag, non_diag), source = 0.d0 )
    ALLOCATE( Rw_blocks(6)%elements(non_diag, non_diag), source = 0.d0 )
    ALLOCATE( Rw_blocks(7)%elements(non_diag, dim_el)  , source = 0.d0 )
    ALLOCATE( Rw_blocks(8)%elements(non_diag, non_diag), source = 0.d0 )
    ALLOCATE( Rw_blocks(9)%elements(non_diag, non_diag), source = 0.d0 )
    
    l = 0
    do m = 1, dim_el
      do n = 1, dim_el
        if (n > m) then
          l = l + 1
          Rw_blocks(6)%elements(l, l) = omega_matrix(n, m) 
        endif 
      enddo
    enddo
    
    
    Rw_blocks(8)%elements = - Rw_blocks(6)%elements 
   
    !----- montando a matriz total --------
    Rw_matrix(1:dim_el, 1:dim_el)                   = Rw_blocks(1)%elements
    Rw_matrix(1:dim_el, dim_el+1:aux1)              = Rw_blocks(2)%elements
    Rw_matrix(1:dim_el, aux1+1:rhosquared)          = Rw_blocks(3)%elements
    
    Rw_matrix(dim_el+1:aux1, 1:dim_el)              = Rw_blocks(4)%elements
    Rw_matrix(dim_el+1:aux1, dim_el+1:aux1)         = Rw_blocks(5)%elements
    Rw_matrix(dim_el+1:aux1, aux1+1:rhosquared)     = Rw_blocks(6)%elements
    
    Rw_matrix(aux1+1:rhosquared, 1:dim_el)          = Rw_blocks(7)%elements
    Rw_matrix(aux1+1:rhosquared, dim_el+1:aux1)     = Rw_blocks(8)%elements
    Rw_matrix(aux1+1:rhosquared, aux1+1:rhosquared) = Rw_blocks(9)%elements
    
    
    !call print_mat2(Rw_matrix, rhosquared, rhosquared) 
    end subroutine createRwMatrix


end module rdftensor_m
