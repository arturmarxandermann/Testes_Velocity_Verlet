module plot_auto_function_m
use types_m
use constants_m
use parameters_m
use functions_m

public :: plot_fem_autofunction

type hermite_type
    real*8 :: x(xypoints)
    real*8 :: y(xypoints)
end type hermite_type

type cte_normalization
     real*8 :: states(max_hmt, max_hmt)
end type cte_normalization


type psi_quantum_oscl_2d
    real*8 :: xy(xypoints, xypoints)
end type psi_quantum_oscl_2d

contains


subroutine calculate_fem_autofunction(sites_array)
implicit none
type(quantum_site), INTENT(IN) :: sites_array(nm_rows, nm_columns)

type(hermite_type), DIMENSION(:, :, :), allocatable :: hmt_site
type(hermite_type), DIMENSION(:, :),  allocatable :: hmt_grid
type(hermite_type), DIMENSION(:, :), allocatable :: exp_part
type(cte_normalization), DIMENSION(:, :), allocatable :: norm_cte !constante denormalizacao para o sitio(i, j) no estado sitio(i, j)%state(k)
type(psi_quantum_oscl_2d), DIMENSION(:, :, :, :), allocatable :: psi_func
!autofuncao do oscilador quantido 2d no sitio(i, j) e estado psi_func(nm_rows,
!nm_columns, k)
real*8, DIMENSION(:, :, :, :), ALLOCATABLE :: z_value

REAL*8, PARAMETER :: scale_nanometro = 1.d9
real*8  :: lesser_number
real*8, DIMENSION(:, :), ALLOCATABLE :: coeffxy
real*8, DIMENSION(:, :), ALLOCATABLE :: hmt_x, hmt_y
real*8, DIMENSION(:), ALLOCATABLE :: xgrid, ygrid
real*8, DIMENSION(:, :), ALLOCATABLE :: ovlp_plot
integer :: lx, ly, i, j, k, nm_arquivo, y_cut


ALLOCATE(xgrid(xypoints), source = 0.d0)
ALLOCATE(ygrid(xypoints), source = 0.d0)
ALLOCATE(coeffxy(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(hmt_x(xypoints, max_hmt), source = 0.d0)
ALLOCATE(hmt_y(xypoints, max_hmt), source = 0.d0)
ALLOCATE(hmt_site(nm_rows, nm_columns, max_hmt) )
ALLOCATE(hmt_grid(nm_rows, nm_columns) )
ALLOCATE(exp_part(nm_rows, nm_columns) )
ALLOCATE(norm_cte(nm_rows, nm_columns) )
ALLOCATE(psi_func(nm_rows, nm_columns, max_hmt, max_hmt) )
ALLOCATE(ovlp_plot(xypoints, xypoints), source = 0.d0)


do i = 1, xypoints
  xgrid(i) = (xmin + (float(i - 1)/float(xypoints)) * (xmax - xmin) )
  ygrid(i) = (ymin + (float(i - 1)/float(xypoints)) * (ymax - ymin) )
 
enddo


!DO PARA TESTAR OS POL DE HERMITE -- CERTOS
!call h_polynomial_value(xypoints, max_hmt, xgrid, hmt_x)
!do j = 1, xypoints
!write(45, 13) xgrid(j), hmt_x(j, 1)
!write(46, 13) xgrid(j), hmt_x(j, 2)
!write(47, 13) xgrid(j), hmt_x(j, 3)
!write(48, 13) xgrid(j), hmt_x(j, 4)
!enddo
!stop


do j = 1, nm_columns
  do i = 1, nm_rows
    coeffxy(i, j) = SQRT( (me  * sites_array(i, j)%omega * thz_to_hz ) / HB_ev )
  enddo
enddo



do j = 1, nm_columns
  do i = 1, nm_rows
    do l = 1, xypoints
      hmt_grid(i, j)%x(l) = coeffxy(i, j) * (xgrid(l) - sites_array(i, j)%posicao_x * scale_nanometro)
      hmt_grid(i, j)%y(l) = coeffxy(i, j) * (ygrid(l) - sites_array(i, j)%posicao_y * scale_nanometro)
    enddo
  enddo
enddo



do j = 1, nm_columns
  do i = 1, nm_rows
    call h_polynomial_value(xypoints, max_hmt, hmt_grid(i, j)%x, hmt_x)
    call h_polynomial_value(xypoints, max_hmt, hmt_grid(i, j)%y, hmt_y)
       do k = 1, max_hmt
        do l = 1, xypoints
       
          hmt_site(i, j, k)%x(l) = hmt_x(l, k) !polinomio de hermite no eixo xpara o sitio (i, j) no estado k = hmt_site_nm(i, j, k) no ponto x = l
          hmt_site(i, j, k)%y(l) = hmt_y(l, k)


          enddo
      enddo
  enddo
enddo




do k = 1, max_hmt
  do l = 1, max_hmt
    do j = 1, nm_columns
      do i = 1, nm_rows
      norm_cte(i, j)%states(k, l) = (1.d0/sqrt( (2.d0**( float(k-1) + float(l-1)))  * ((GAMMA(float(l)))**2.d0) * &
   ((GAMMA(float(k)))**2.d0)  )) * SQRT( (me * sites_array(i, j)%omega * thz_to_hz ) / (PI * HB_ev) )
      enddo
    enddo
  enddo
enddo


do j = 1, nm_columns
  do i = 1, nm_rows
    do l = 1, xypoints
      exp_part(i, j)%x(l) = exp(-(me/(2.d0 * HB_ev)) * sites_array(i, j)%omega * &
thz_to_hz  * (xgrid(l) - sites_array(i, j)%posicao_x * scale_nanometro)**2.0 )
      exp_part(i, j)%y(l) = exp(-(me/(2.d0 * HB_ev)) * sites_array(i, j)%omega * &
thz_to_hz *  (ygrid(l) - sites_array(i, j)%posicao_y * scale_nanometro)**2.0 )
    enddo
  enddo
enddo


!CÁLCULO DO PSI_FUNC: FUNCAO DE ONDA PARA O OSCILADOR BIDIMENSIONAL TAL QUE
!PSI_FUNC( PONTO EM X, PONTO EM Y,
!N_AUTOFUNCAO)%XY(X_VALUE, Y_VALUE)
do k = 1, max_hmt
  do l = 1, max_hmt
    do j = 1, nm_columns
      do i = 1, nm_rows
        do lx = 1, xypoints
          do ly = 1, xypoints
          psi_func(i, j, k, l)%xy(lx, ly) = norm_cte(i, j)%states(k, l) * &
hmt_site(i, j, k)%x(lx) * hmt_site(i, j, l)%y(ly) * exp_part(i, j)%x(lx) * exp_part(i, j)%y(ly)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo


call plot_fem_autofunction(xgrid, ygrid, psi_func, z_value)


!ARQUIVO PARA O PLOT EM DUAS DIMENSOES

open(40, status = 'replace')
do lx = 1, xypoints
  do ly = 1, xypoints
    write(40, 13) xgrid(lx), ygrid(ly), z_value(lx, ly, 1, 1)
  enddo
    write(40, 13)
enddo
close(40)
!


!ARQUIVOS PARA O CORTE EM Y = 0 VARIANDO EM X
nm_arquivo = 0
y_cut = xypoints/2
do j = 1, nm_columns
  do i = 1, nm_rows
    open(50 + nm_arquivo, status = 'replace')
      do lx = 1, xypoints
        write(50 + nm_arquivo, 13) xgrid(lx) , psi_func(i, j, 1, 1)%xy(lx, y_cut+ 1)
      enddo
    close(50 + nm_arquivo)
    nm_arquivo = nm_arquivo + 1
  enddo
enddo



!nm_arquivo = 0
!k = 0
!do j = 1, nm_columns
!  do i = 1, nm_rows
!    open(70 + nm_arquivo)
!      do lx = 1, xypoints
!        k = j + 1
!        if (k == nm_columns+1) exit
!          lesser_number = min(psi_func(i, j, 1)%xy(lx, y_cut), psi_func(i, k,
!          1)%xy(lx, y_cut))
!          write(70 + nm_arquivo, 13) xgrid(lx), lesser_number
!      enddo
!    close(70 + nm_arquivo)
!    nm_arquivo  = nm_arquivo + 1
!  enddo
!enddo







13 format(60F14.5)

deallocate(xgrid, ygrid, coeffxy, hmt_x, hmt_y, hmt_site, hmt_grid, exp_part, norm_cte, psi_func, ovlp_plot)
return
end subroutine calculate_fem_autofunction

subroutine plot_fem_autofunction(xgrid, ygrid, psi_func, z_value)
implicit none
real*8, intent(in) :: xgrid(xypoints), ygrid(xypoints)
type(psi_quantum_oscl_2d), intent(in) :: psi_func(nm_rows, nm_columns, max_hmt, max_hmt)
real*8,  ALLOCATABLE, intent(out) :: z_value(:, :, :, :)

real*8, parameter :: scale_nanometro = 1.d9
real*8, ALLOCATABLE :: vectors_phi_values(:), vectors_phi_values_abs(:)
integer :: lx, ly, s, k_value, k, l
real*8, parameter :: lim_distance = 100.d0 ! só vejo qual f(z) é maior para sitios dentro de 2 nm

ALLOCATE(z_value(xypoints, xypoints, max_hmt, max_hmt), source = 0.d0)
ALLOCATE(vectors_phi_values(20), source = 0.d0)
ALLOCATE(vectors_phi_values_abs(20), source = 0.d0)


do s = 1, max_hmt
do l = 1, max_hmt
  k = 0 
  do lx = 1, xypoints
  do ly = 1, xypoints
    k = 0
      do j = 1, nm_columns
      do i = 1, nm_rows
              k = k + 1
              vectors_phi_values(k) = psi_func(i, j, s, l)%xy(lx, ly)
              vectors_phi_values_abs(k) = ABS(vectors_phi_values(k))
        enddo
      enddo
      k_value =  MAXLOC(vectors_phi_values_abs, dim=1)  !o valor de k vai ser o da maior (em modulo) funcao fem
    z_value(lx, ly, s, l) = vectors_phi_values(k_value)
  enddo
enddo
enddo
enddo

deallocate(vectors_phi_values, vectors_phi_values_abs)
13 format(60F14.4)
return
end subroutine plot_fem_autofunction

subroutine h_polynomial_value ( m, n, x, p )

  !*****************************************************************************80
  !
  !! H_POLYNOMIAL_VALUE evaluates H(i,x).
  !
  !  Discussion:
  !
  !    H(i,x) is the physicist's Hermite polynomial of degree I.
  !
  !  Differential equation:
  !
  !    Y'' - 2 X Y' + 2 N Y = 0
  !
  !  First terms:
  !
  !      1
  !      2 X
  !      4 X^2     -  2
  !      8 X^3     - 12 X
  !     16 X^4     - 48 X^2     + 12
  !     32 X^5    - 160 X^3    + 120 X
  !     64 X^6    - 480 X^4    + 720 X^2    - 120
  !    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
  !    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
  !    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
  !   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
  !
  !  Recursion:
  !
  !    H(0,X) = 1,
  !    H(1,X) = 2*X,
  !    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
  !
  !  Norm:
  !
  !    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
  !    = sqrt ( PI ) * 2^N * N!
  !
  !    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    01 October 2002
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Milton Abramowitz, Irene Stegun,
  !    Handbook of Mathematical Functions,
  !    National Bureau of Standards, 1964,
  !    ISBN: 0-486-61272-4,
  !    LC: QA47.A34.
  !
  !    Larry Andrews,
  !    Special Functions of Mathematics for Engineers,
  !    Second Edition,
  !    Oxford University Press, 1998.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) M, the number of evaluation points.
  !
  !    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
  !    Note that polynomials 0 through N will be computed.
  !
  !    Input, real ( kind = 8 ) X(M), the evaluation points.
  !
  !    Output, real ( kind = 8 ) P(M,0:N), the values of the first N+1 Hermite
  !    polynomials at the point X.
  !
implicit none

integer, INTENT(IN) :: m
integer, INTENT(IN) :: n
real*8, dimension(m), INTENT(IN) :: x
real*8, dimension(m, n), INTENT(OUT) :: p

integer :: j


if ( n <= 0 ) then
return
end if

p(1:m,1) = 1.0D+00

if ( n == 0 ) then
return
end if

p(1:m,2) = 2.0D+00 * x(1:m)

do j = 3, n
p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
               - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
end do

return
end subroutine h_polynomial_value



end module plot_auto_function_m
