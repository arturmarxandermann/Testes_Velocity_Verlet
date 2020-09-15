module system_operators_m
use constants_m
use types_m
use overlap_m
use parameters_m
use functions_m 

contains


subroutine print_mat3(aa, nn, mm)
implicit none
integer, intent(in) :: nn, mm
real*8, dimension(nn, mm), intent(in) :: aa
integer :: i, j
do, i=1,mm
    write(*,'(100g10.3)') ( aa(i,j), j=1,nn )
enddo

end subroutine print_mat3


subroutine create_system_operators(el_or_hl, psize, phi, ovlpm, OP_system_ham)
!subrotina para criar os operadores Q^\alpha para eletrons e buracos
implicit none
integer, INTENT(IN) :: el_or_hl, psize
real*8, dimension(psize, psize), intent(in) :: phi
!real*8, dimension(psize, psize), intent(in) :: Ztrm
!real*8, dimension(psize, psize), intent(in) :: SInverse
real*8, dimension(psize, psize) :: ovlpm 
TYPE(operators), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: OP_system_ham


TYPE(operators), DIMENSION(:), ALLOCATABLE :: OP_system_sites


real*8, dimension(:, :), allocatable :: phi_transpose 
real*8, DIMENSION(:, :), ALLOCATABLE :: rec_matrix, CouplingInterSites


real*8 :: taxa_dephasing
real*8, allocatable, dimension(:, :) :: S, OVLPBlocks

integer :: dims, lim_sup, lim_inferior, nstates
integer :: r1, c1, r2, c2 
integer :: valinf, valmax


if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if


if (el_or_hl == 1 ) then
  dims = dim_el
  nstates = nmstates_el
  taxa_dephasing = v_dephasing_el
  !ALLOCATE(OP_system_sites(nsites))
  !ALLOCATE(OP_system_ham(nsites))
  ALLOCATE(OP_system_sites(half_ndim))
  ALLOCATE(OP_system_ham(half_ndim))
  do i = 1, half_ndim
    ALLOCATE(OP_system_sites(i)%elements(dim_el, dim_el), source = 0.d0)
    ALLOCATE(OP_system_ham(i)%elements(dim_el, dim_el), source = 0.d0)
  enddo
end if

if (el_or_hl == 2 ) then
  dims = dim_hl
  nstates = nmstates_hl
  taxa_dephasing = v_dephasing_hl
  
  !ALLOCATE(OP_system_sites(nsites))
  !ALLOCATE(OP_system_ham(nsites))
  ALLOCATE(OP_system_sites(half_ndim))
  ALLOCATE(OP_system_ham(half_ndim))
  do i = 1, half_ndim
    ALLOCATE(OP_system_sites(i)%elements(dim_hl, dim_hl), source = 0.d0)
    ALLOCATE(OP_system_ham(i)%elements(dim_hl, dim_hl), source = 0.d0)
  enddo
end if


ALLOCATE(rec_matrix(nstates, 1), source = 0.d0)
ALLOCATE(phi_transpose(dims, dims), source = 0.d0)
ALLOCATE(S(nstates, nstates), source = 0.d0 ) 
ALLOCATE(OVLPBlocks(nstates, nstates), source = 0.d0 ) 
ALLOCATE(CouplingInterSites(nstates, nstates), source = 0.d0 ) 


phi_transpose = transpose(phi) 




!======== PARTE NÂO DIAGONAL DOS BLOCOS DIAGONAIS ============ 
do i = 1, nstates
CouplingInterSites = taxa_dephasing 
enddo

forall(k = 1:nstates) CouplingInterSites(k, k) = 0.d0

k = 0
do alpha = 1, nsites
  OP_system_sites(alpha)%elements(k*nstates+1:k*nstates+nstates, k*nstates+1:k*nstates+nstates) = &
                                                                               RESHAPE(CouplingInterSites, (/nstates, nstates/))
  k = k + 1
enddo
 
   
!===========================================================

!==== BLOCOS NÂO DIAGONAIS -> TERMOS DE RELAXAÇÃO ===========
alpha = nsites+1

OP_system_sites(alpha)%elements(:, :) = ovlpm * taxa_dephasing

forall(i = 1 : dims) OP_system_sites(alpha)%elements(i, i) = 0.d0

!========================================================================================================



!PASSO DA BASE DOS SITIOS PARA A BASE DO HAMILTONIANO --------
do alpha = 1, half_ndim
OP_system_ham(alpha)%elements(:, :) = matmul(phi_transpose, matmul(OP_system_sites(alpha)%elements(:, :), phi))
!OP_system_ham(alpha)%elements(:, :) = OP_system_sites(alpha)%elements(:, :)
enddo


!call printaletters2('op_si1', 5, dims, dims, dims, OP_system_ham(1)%elements(:, :) )
!call printaletters2('op_si2', 6, dims, dims, dims, OP_system_ham(2)%elements(:, :) )
!call printaletters2('op_si3', 7, dims, dims, dims, OP_system_ham(3)%elements(:, :) )
!!call printaletters2('op_si4', 8, dims, dims, dims, OP_system_ham(4)%elements(:, :) )


!stop 

!do alpha = 1, half_ndim
!    print*, "OP_system", alpha
!    call print_mat3(OP_system_sites(alpha)%elements(:, :), dims, dims)
!enddo
!stop 


DEALLOCATE(rec_matrix, phi_transpose, OP_system_sites, S, CouplingInterSites, OVLPBlocks)
return

end subroutine create_system_operators









end module system_operators_m
