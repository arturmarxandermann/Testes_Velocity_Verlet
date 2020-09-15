module rdftensor_m
use f95_precision
use lapack95
use blas95
use types_m
use parameters_m
use constants_m
use functions_m
use system_operators_m
use omp_lib 
contains


function testFT(frequencia)
implicit none
real*8, intent(in) :: frequencia
real*8 :: testFT
real*8 :: beta_boltz
real*8 :: cte

beta_boltz = 1.d0/(kboltz*temp)
cte = (HB_ev_ps*beta_boltz)/2.d0

if (frequencia < 0.d0) then
testFT = coupling * (abs(frequencia))*exp(HB_ev_ps*frequencia*beta_boltz)*(cosh((abs(frequencia) * cte)) / sinh( abs(frequencia) * cte) )
endif

if (frequencia > 0.d0) then
testFT = coupling * frequencia*(cosh(frequencia * cte)/ sinh(frequencia * cte))
endif

if (frequencia .EQ. 0.d0) then  !lim x->0  x(cotgh(x/t) = t ... t = 2*coupling*kT
testFT = ( (2.d0 * coupling * kboltz * temp) / HB_ev_ps )
endif


return
end function testFT

pure function FT_transf(frequencia)
implicit none
real*8, intent(in) :: frequencia
real*8 :: FT_transf
real*8 :: beta_boltz
real*8 :: cte

beta_boltz = 1.d0/(kboltz*temp)
cte = (HB_ev_ps*beta_boltz)/2.d0

if (frequencia < 0.d0) then
FT_transf = coupling * abs(frequencia)*exp(HB_ev_ps*frequencia*beta_boltz) * (cosh((abs(frequencia) * cte)) / sinh( abs(frequencia) * cte) )
endif

if (frequencia > 0.d0) then
FT_transf = coupling * frequencia*(cosh(frequencia * cte)/ sinh(frequencia * cte))
endif

if (frequencia .EQ. 0.d0) then  !lim x->0  x(cotgh(x/t) = t ... t = 2*coupling*kT
FT_transf = ( (2.d0 * coupling * kboltz * temp ) / HB_ev_ps )
endif

return
end function FT_transf


subroutine create_FT_func
implicit none
real*8 :: omega_min, omega_max
real*8, ALLOCATABLE, DIMENSION(:) :: omega_grid, funcFT

integer, parameter :: wpoints = 10000

ALLOCATE(omega_grid(wpoints), source = 0.d0 )
allocate(funcFT(wpoints), source = 0.d0 )

omega_min = -7000.d0
omega_max =  7000.d0

do i = 1, wpoints
  omega_grid(i) = omega_min + (float(i - 1)/float(wpoints)) * (omega_max - omega_min)
  funcFT(i) = testFT(omega_grid(i))
enddo

open (unit=5, file='FTtransf', action="write", status="replace")

do i = 1, wpoints
  write(5, '(F9.3, F9.3)') omega_grid(i), funcFT(i)
enddo

CLOSE(5)

DEALLOCATE(omega_grid, funcFT)
return
end subroutine create_FT_func



subroutine create_redfield_tensor(el_or_hl, ham_size, autovetores, omega_ab, ovlpm, rdf_tensor)
implicit none
integer, INTENT(IN) :: el_or_hl, ham_size
real*8, DIMENSION(ham_size, ham_size), INTENT(IN) :: autovetores
real*8, DIMENSION(ham_size, ham_size), INTENT(IN) :: omega_ab
!real*8, dimension(ham_size, ham_size), intent(in) :: Ztrm 
!real*8, dimension(ham_size, ham_size), intent(in) :: IS_ovlpm
real*8, dimension(ham_size, ham_size), intent(in) :: ovlpm
real*8, DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(OUT) :: rdf_tensor

real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: gammamais, gammamenos
real*8, DIMENSION(:, :), ALLOCATABLE :: temporaria, autovetoresTR
type(operators), DIMENSION(:), ALLOCATABLE :: OP_system
!type(operators), DIMENSION(:), ALLOCATABLE :: OP_system_mult
!real*8, DIMENSION(:, :, :, :), ALLOCATABLE :: OP_system_sum
integer :: a, b, c, d, alpha, i, j, k, l
real*8, dimension(:, :, :), allocatable :: tempo, temporary 
real*8 :: OP_system_number, soma, soma1, soma2
real*8 :: t1, t2, t3, t4, t5, t6
integer :: thread_num
real*8, dimension(:, :), allocatable :: FT_transf_mais, FT_transf_menos 

!set OMP_NUM_THREADS = 90

if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if

if (el_or_hl == 1 ) then
  dims = dim_el
  call create_system_operators(1, dims, autovetores, ovlpm, OP_system)
end if

if (el_or_hl == 2 ) then
  dims = dim_hl
  call create_system_operators(2, dims, autovetores, ovlpm, OP_system)
end if

!ALLOCATE(OP_system_mult(nsites))
!do alpha = 1, nsites
!  ALLOCATE(OP_system_mult(alpha)%elements(dims, dims))
!enddo

ALLOCATE(gammamais(ham_size, ham_size, ham_size, ham_size), source = 0.d0)
ALLOCATE(gammamenos(ham_size, ham_size, ham_size, ham_size), source = 0.d0)
ALLOCATE(temporaria(ham_size, ham_size), source = 0.d0)
ALLOCATE(autovetoresTR(ham_size, ham_size), source = 0.d0 )
!ALLOCATE(OP_system_sum(ham_size, ham_size, ham_size, ham_size), source = 0.d0)
!allocate(tempo(nsites, ham_size, ham_size))
allocate(tempo(half_ndim, ham_size, ham_size))
allocate(temporary(nsites, ham_size, ham_size))
ALLOCATE(rdf_tensor(ham_size, ham_size, ham_size, ham_size), source = 0.d0)
ALLOCATE(FT_transf_mais(ham_size, ham_size), source = 0.d0)
ALLOCATE(FT_transf_menos(ham_size, ham_size), source = 0.d0)
autovetoresTR = TRANSPOSE(autovetores)




!DIR& PARALLEL
do alpha = 1, half_ndim
   !tempo(alpha, :, :) = MATMUL(autovetoresTR, MATMUL(OP_system(alpha)%elements,   autovetores))
   tempo(alpha, :, :) = OP_system(alpha)%elements(:, :)
enddo 

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(d, c) SHARED(FT_transf_mais, FT_transf_menos)
do d = 1, dims
  do c = 1, dims
  FT_transf_mais(c, d) = FT_transf(omega_ab(d, c) )
  FT_transf_menos(c, d) = FT_transf(omega_ab(c, d) ) 
  enddo
enddo
!$OMP END PARALLEL DO 

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(d, c, b, a) SHARED(tempo, gammamais, gammamenos)
do d = 1, dims
  do c = 1, dims
    do b = 1, dims
      do a = 1, dims
        OP_system_number = SUM( tempo(:, a, b) * tempo(:, c, d) ) 
        !gammamais(a, b, c, d) = half * FT_transf(omega_ab(a, b)) * OP_system_number
        !gammamenos(a, b, c, d) = half * FT_transf(omega_ab(d, c)) * OP_system_number 
        gammamais(a, b, c, d) = half * FT_transf_mais(c, d) * OP_system_number
        gammamenos(a, b, c, d) = half * FT_transf_menos(a, b) * OP_system_number 
      enddo
    enddo
  enddo
enddo
!$OMP END PARALLEL DO 

!!DIR$ PARALLEL
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(d, c, b, a) SHARED(rdf_tensor, gammamais, gammamenos)
do d = 1, dims
  do c = 1, dims
    do b = 1, dims
      do a = 1, dims
        if (a /= c .AND. d /= b) then
          rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c)
        endif
        if (a /= c .AND.  d == b) then
         ! rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c) &
         ! - soma_elemento(ham_size, a, c, gammamenos)
          rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c) &
          - soma_elemento(ham_size, a, c, gammamais)
       
        endif
        if (a == c .AND. d /= b) then
          !rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c) &
          !- soma_elemento(ham_size, d, b, gammamais)
          rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c) &
          - soma_elemento(ham_size, d, b, gammamenos)
        endif
        if (a == c .AND. d == b) then
          rdf_tensor(a, b, c, d) = gammamais(d, b, a, c) + gammamenos(d, b, a, c) &
          - soma_elemento(ham_size, a, c, gammamais) - soma_elemento(ham_size, d, b, gammamenos)
        endif
      enddo
    enddo
  enddo
enddo
!$OMP END PARALLEL DO 

!print*, "rdf tensor calculado"
!call cpu_time(t6)

!print*, "tempo para criar matrix tempo:", t2 - t1
!print*, "tempo para criar gammais e gammamenos:", t4 - t3
!print*, "tempo para criar rdftensor", t6 - t5

return
deallocate(tempo, gammamais, gammamenos, temporaria, autovetoresTR, OP_system, FT_transf_mais, FT_transf_menos)
end subroutine create_redfield_tensor

subroutine createRwMatrix(rhosize, omega_matrix, rtensor, Rw_matrix)
!subrotina que cria a matriz Rw que relaciona os vetores rhodot e rho
!na forma rho_dot = matmul(Rw, rho).
!Rw é composta por 9 blocos de matrizes, que são numerados na seguinte ordem
! Rw = ( ( 1 ) ( 2 )  ( 3 )  )
!      ( ( 4 ) ( 5 )  ( 6 )  )
!      (  (7 ) ( 8 )  ( 9 )  )
!Como o bloco 1, 4 e 7 relaciona elementos rhoreal_aa, rhoreal_ab, rhoimag_ab
!com rhoreal_aa
!seu tamanho é rhosize. Os outros blocos tem tamanho non_diag pois relacionam os
!elementos
!de antes com rhoreal_ab, rhoimag_ab


implicit none
integer, INTENT(IN) :: rhosize
REAL*8, INTENT(IN), DIMENSION(rhosize, rhosize) :: omega_matrix
REAL*8, INTENT(IN), DIMENSION(rhosize, rhosize, rhosize, rhosize) :: rtensor
REAL*8, INTENT(OUT), ALLOCATABLE, DIMENSION(:, :) :: Rw_matrix

INTEGER :: non_diag, rhosizesquared, aux1, counter, i, j, k, l, n
TYPE(tensor_product_matrices), DIMENSION(9) :: Rw_blocks
!TRmtc são matrizes utilziadas no produto tensorial para montar Rw



non_diag = (rhosize * ( rhosize -1 ) ) / 2
rhosizesquared = rhosize * rhosize
aux1 = rhosize + non_diag




ALLOCATE(Rw_matrix(rhosizesquared, rhosizesquared), source = 0.d0 )


ALLOCATE(Rw_blocks(1)%elements(rhosize, rhosize), source = 0.d0 )
ALLOCATE(Rw_blocks(2)%elements(rhosize, non_diag), source = 0.d0 )
ALLOCATE(Rw_blocks(3)%elements(rhosize, non_diag), source = 0.d0 )

ALLOCATE(Rw_blocks(4)%elements(non_diag, rhosize), source = 0.d0 )
ALLOCATE(Rw_blocks(5)%elements(non_diag, non_diag), source = 0.d0 )
ALLOCATE(Rw_blocks(6)%elements(non_diag, non_diag), source = 0.d0 )

ALLOCATE(Rw_blocks(7)%elements(non_diag, rhosize), source = 0.d0 )
ALLOCATE(Rw_blocks(8)%elements(non_diag, non_diag), source = 0.d0 )
ALLOCATE(Rw_blocks(9)%elements(non_diag, non_diag), source = 0.d0 )


!ALLOCATE(popvec_rec(nsites), source = 0.d0 ) 





!---- bloco 1 ------

!!$OMP PARALLEL DO PRIVATE(j, i) SHARED(rtensor, Rw_blocks) 
do j = 1, rhosize
  do i = 1, rhosize
        Rw_blocks(1)%elements(i, j) = rtensor(i, i, j, j) 
  enddo
enddo
!!$OMP END PARALLEL DO 

!-----------------



!--- bloco 2 ------
!!$OMP PARALLEL DO PRIVATE(i, j, n) SHARED(rtensor, Rw_blocks)
n = 0
do i = 1, rhosize
    do j = 1, rhosize
      if (j > i ) then
        n = n + 1
        do l = 1, rhosize
        Rw_blocks(2)%elements(l, n) = rtensor(l, l, i, j) + rtensor(l, l, j, i)
        enddo
      endif
    enddo
enddo
!!$OMP END PARALLEL DO
!----------------





!--- bloco 3 -------
!bloco 3 é composto de zeros, já definido com o source
!-------------------

!--- bloco 4 --------

!!$OMP PARALLEL DO PRIVATE(i, j) SHARED(rtensor, Rw_blocks)
do n = 1, rhosize
  k = 0
  do i = 1, rhosize
    do j = 1, rhosize
      if (j > i ) then
        k = k + 1
        Rw_blocks(4)%elements(k, n) = rtensor(i, j, n, n)
      endif
    enddo
  enddo
enddo
!!$OMP END PARALLEL DO 

!--------------------


!--- blocos 5, 9, 6, 8, bloco 7 é zero


l = 0
counter = 0


do i = 1, rhosize
  do j = 1, rhosize
      l = 0
    if (j > i) then
       counter = counter + 1

      !!$OMP PARALLEL DO PRIVATE(m, n) ORDERED SCHEDULE(DYNAMIC) SHARED(rtensor, Rw_blocks, omega_matrix)
      do m = 1, rhosize
        do n = 1, rhosize
          if (n > m) then
            !counter = counter + 1
            l = l + 1
            !!$OMP ORDERED
            Rw_blocks(5)%elements(l, counter) = rtensor(m, n, i, j) + rtensor(m, n, i, j)
            Rw_blocks(9)%elements(l, counter) = rtensor(m, n, i, j) - rtensor(m, n, i, j) 
            !!$OMP END ORDERED
              if (l == counter) then
                !l == counter => matrizes diagonais
                 !!$OMP ORDERED
                 Rw_blocks(6)%elements(l, l) = omega_matrix(n, m) 
                 !!$OMP END ORDERED
              endif
          endif
         enddo
       enddo
       !!$OMP END PARALLEL DO 
    endif
  enddo
enddo
!!$OMP END PARALLEL DO


Rw_blocks(8)%elements = - Rw_blocks(6)%elements 


!----- montando a matriz total --------
Rw_matrix(1:rhosize, 1:rhosize) = Rw_blocks(1)%elements
Rw_matrix(1:rhosize, rhosize+1:aux1) = Rw_blocks(2)%elements
Rw_matrix(1:rhosize, aux1+1:rhosizesquared) = Rw_blocks(3)%elements

Rw_matrix(rhosize+1:aux1, 1:rhosize) = Rw_blocks(4)%elements
Rw_matrix(rhosize+1:aux1, rhosize+1:aux1) = Rw_blocks(5)%elements
Rw_matrix(rhosize+1:aux1, aux1+1:rhosizesquared) = Rw_blocks(6)%elements

Rw_matrix(aux1+1:rhosizesquared, 1:rhosize) = Rw_blocks(7)%elements
Rw_matrix(aux1+1:rhosizesquared, rhosize+1:aux1) = Rw_blocks(8)%elements
Rw_matrix(aux1+1:rhosizesquared, aux1+1:rhosizesquared) = Rw_blocks(9)%elements


!call print_mat2(Rw_matrix, rhosizesquared, rhosizesquared)


return
end subroutine createRwMatrix



end module rdftensor_m
