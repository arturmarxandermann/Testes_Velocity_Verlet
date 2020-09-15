module edo_solver_m
use f95_precision
use lapack95
use blas95
use functions_m
use constants_m
use parameters_m
use types_m

contains

subroutine createRwMatrix(rhosize, omega_matrix, rtensor, Rw_matrix)
!subrotina que cria a matriz Rw que relaciona os vetores rhodot e rho
!na forma rho_dot = matmul(Rw, rho).
!Rw é composta por 9 blocos de matrizes, que são numerados na seguinte ordem
! Rw = ( ( 1 ) ( 2 )  ( 3 )  )
!      ( ( 4 ) ( 5 )  ( 6 )  )
!      (  (7 ) ( 8 )  ( 9 )  )
!Como o bloco 1, 4 e 7 relaciona elementos rhoreal_aa, rhoreal_ab, rhoimag_ab com rhoreal_aa
!seu tamanho é rhosize. Os outros blocos tem tamanho non_diag pois relacionam os elementos
!de antes com rhoreal_ab, rhoimag_ab


implicit none
integer, INTENT(IN) :: rhosize
REAL*8, INTENT(IN), DIMENSION(rhosize, rhosize) :: omega_matrix
REAL*8, INTENT(IN), DIMENSION(rhosize, rhosize, rhosize, rhosize) :: rtensor
REAL*8, INTENT(OUT), ALLOCATABLE, DIMENSION(:, :) :: Rw_matrix

INTEGER :: non_diag, rhosizesquared, aux1, counter
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


!---- bloco 1 ------
forall(i=1:rhosize, j=1:rhosize) Rw_blocks(1)%elements(i, j) = rtensor(i, i, j, j)
!-----------------


!--- bloco 2 ------
do n = 1, rhosize
  k = 0
  do i = 1, rhosize
    do j = 1, rhosize
      if (j > i ) then
        k = k + 1
        Rw_blocks(2)%elements(n, k) = rtensor(n, n, i, j) + rtensor(n, n, j, i)
      endif
    enddo
  enddo
enddo
!----------------

!--- bloco 3 -------
!bloco 3 é composto de zeros, já definido com o source
!-------------------

!--- bloco 4 --------
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
!--------------------

!--- blocos 5, 9, 6, 8, bloco 7 é zero
l = 0
counter = 0
do i = 1, rhosize
  do j = 1, rhosize
    if (j > i) then
      counter = 0
      l = l + 1
      do m = 1, rhosize
        do n = 1, rhosize
          if (n > m) then
            counter = counter + 1
            Rw_blocks(5)%elements(l, counter) = rtensor(i, j, m, n) + rtensor(i, j, n, m)
            Rw_blocks(9)%elements(l, counter) = rtensor(i, j, m, n) - rtensor(i, j, n, m)
              if (l == counter) then
                !l == counter => matrizes diagonais
                 Rw_blocks(6)%elements(l, counter) = + omega_matrix(m, n)
                 Rw_blocks(8)%elements(l, counter) = - omega_matrix(m, n)
              endif
          endif
        enddo
      enddo
    endif
  enddo
enddo
!---------------------


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
!----------------------------------------

!call print_mat2(Rw_matrix, rhosizesquared, rhosizesquared)

return
end subroutine createRwMatrix

subroutine CoupledEDOSolver(Msize, t_in, t_out, M_coeffi, vec_in, vec_out)
implicit none
integer, intent(in) :: Msize !Msize = neqn
real*8, INTENT(IN) ::t_in, t_out
real*8, dimension(Msize, Msize), intent(in) :: M_coeffi
real*8, dimension(Msize), intent(in) :: vec_in
real*8, dimension(Msize), intent(out) :: vec_out



integer :: info
complex*16, ALLOCATABLE, DIMENSION(:) :: EGvalues
complex*16, ALLOCATABLE, dimension(:, :) :: M_coeffi_temp, EGvectors, vr
complex*16, ALLOCATABLE, DIMENSION(:, :) :: exp_EGsystem_in, exp_EGsystem_out, vec_in_temp, coeffi_INITime

real*8, ALLOCATABLE, DIMENSION(:) :: REGvalues, IEGvalues
REAL*8, ALLOCATABLE, DIMENSION(:, :) :: REGexp_EGsystem_in, IEGexp_EGsystem_in, REGexp_EGsystem_out, IEGexp_EGsystem_out

real*8, ALLOCATABLE, DIMENSION(:, :) :: REGvectors, IEGvectors

real*8, parameter :: lim_EGvec = 705.d0

ALLOCATE(exp_EGsystem_in(Msize, Msize) , source = (0.d0, 0.d0) )
ALLOCATE(exp_EGsystem_out(Msize, Msize), source = (0.d0, 0.d0) )

ALLOCATE(REGexp_EGsystem_in(Msize, Msize), source = (0.d0) )
ALLOCATE(IEGexp_EGsystem_in(Msize, Msize), source = (0.d0) )

ALLOCATE(REGvectors(Msize, Msize), source = 0.d0 )
ALLOCATE(IEGvectors(Msize, Msize), source = 0.d0 )

ALLOCATE(REGexp_EGsystem_out(Msize, Msize), source = (0.d0) )
ALLOCATE(IEGexp_EGsystem_out(Msize, Msize), source = (0.d0) )


ALLOCATE(coeffi_INITime(Msize, Msize) , source = (0.d0, 0.d0) )
ALLOCATE(vec_in_temp(Msize, 1) , source = (0.d0, 0.d0) )

ALLOCATE(EGvectors(Msize, Msize), source = (0.d0, 0.d0) )
ALLOCATE(EGvalues(Msize), source = (0.d0, 0.d0) )
ALLOCATE(REGvalues(Msize), source = 0.d0 )
ALLOCATE(IEGvalues(Msize), source = 0.d0 )

ALLOCATE(M_coeffi_temp(Msize, Msize), source = (0.d0, 0.d0) )
ALLOCATE(vr(Msize, Msize), source = (0.d0, 0.d0) )

forall (i=1:Msize) vec_in_temp(i, 1) = vec_in(i)
M_coeffi_temp = M_coeffi

info = 0

!------------------------------------------------------------------------------
!PASSO 1: DIAGONALIZAR M_coeffi_temp

call geev(M_coeffi_temp, EGvalues, vr, EGvectors, info  )


!geev encontra os autovetores de uma matriz generica real M_coeffi_temp
!os autovetores sao normalizados e armazenados na matriz EGvectors
!a parte real/imaginaria dos autovalores sao armazenadas em R_EGvalues/I_EGvalues
!vr armazena os autovalores em uma forma "right side" => nao utilizo

if ( info /= 0 ) write(*,*) info, "Error to found eigenvalues"
!
REGvalues = real(EGvalues)
IEGvalues = AIMAG(EGvalues)
!!
REGvectors = real(EGvectors)
IEGvectors = aimag(EGvectors)



do j = 1, Msize
    if ( t_in * REGvalues(j) <= -lim_EGvec) then
      REGvalues(j) = -lim_EGvec/t_in
    endif
    if ( t_in * IEGvalues(j) <= -lim_EGvec) then
      IEGvalues(j) = -lim_EGvec/t_in
    endif

  EGvalues(j) = REGvalues(j) + zi * IEGvalues(j)
enddo




forall(i=1:Msize, j=1:Msize)  exp_EGsystem_in(i, j)  = exp(EGvalues(j) * t_in ) * EGvectors(i, j)


!---------------------------------------------------------------------------
!PASSO 2: calcular as matrizes do sistema de equacoes

!exp_EGsystem_in é a matriz com coeficientes que será utilizada para o cálculo
!dos coeficientes iniciais C_n^i do conjunto de equacoes lineares acopladas


coeffi_INITime = exp_EGsystem_in


call gesv(coeffi_INITime, vec_in_temp)



!subrotina para resolver o sistema de equacoes lineares A*X = B
!B pode ter vários RHS => ver no mkl, aqui deixamos apenas com 1 RHS
!Dessa forma a matriz vec_in_temp tem uma coluna => vec_in_temp(Msize, 1)
!a resposta dos coeficientes C1, C2, .. , CN são reescritos na matriz vec_in_temp
if (info /= 0) then
        print*, "Something went wrong in calculating the System of Linear Equations"
        stop
endif

do j = 1, Msize
    if ( t_out * REGvalues(j) <= -lim_EGvec) then
      REGvalues(j) = -lim_EGvec/t_out
    endif
    if ( t_out * IEGvalues(j) <= -lim_EGvec) then
      IEGvalues(j) = -lim_EGvec/t_out
    endif

  EGvalues(j) = REGvalues(j) + zi * IEGvalues(j)
enddo


forall(i=1:Msize, j=1:Msize)  exp_EGsystem_out(i, j)  = exp(EGvalues(j) * t_out) * EGvectors(i, j)

!exp_EGsystem_out = exp_EGsystem_out * exp(t_in)

!------------------------------------------------------------------------

!PASSO 3: calcula os vetores de saída vec_out

forall(i=1:Msize) vec_out(i) = real(sum( vec_in_temp(:, 1) * exp_EGsystem_out(i, :)  ))




DEALLOCATE(exp_EGsystem_in, exp_EGsystem_out, coeffi_INITime, EGvectors, M_coeffi_temp &
          , vr, vec_in_temp, EGvalues, REGexp_EGsystem_in, REGexp_EGsystem_out, IEGexp_EGsystem_in &
          , REGvectors, IEGvectors, IEGexp_EGsystem_out)
return
end subroutine CoupledEDOSolver


end module edo_solver_m
