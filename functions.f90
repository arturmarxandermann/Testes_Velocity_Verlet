module functions_m
use f95_precision
use lapack95
use blas95
use omp_lib
use types_m
use parameters_m
use constants_m !usa as constantes do constants_m => variaveis globais.

implicit none

public

contains

subroutine print_matrices(step, rhoSites_in, phi, hMtx) 
    implicit none
    
    !args
    integer,    intent(in) :: step
    complex*16, intent(in) :: rhoSites_in(d_el,d_el)
    real*8,     intent(in) :: phi(d_el, d_el), hMtx(d_el,d_el)

    !local
	  real*8, allocatable    :: rhomtx(:,:)

    allocate(rhomtx(d_el,d_el), source = 0.d0 ) 
    
    if ( step == 1 .OR. step == nm_divisoes/2 .OR. step == nm_divisoes ) then
      rhomtx = 0.d0
      rhomtx = real(rhoSites_in)
      print*, "PARTE REAL COMECO DO PROGRAMA"
      call print_mat2(rhomtx, d_el, d_el)
      rhomtx = 0.d0
      rhomtx = aimag(rhoSites_in)
      print*, "PARTE IMAGINARIA COMECO DO PROGRAMA"
      call print_mat2(rhomtx, d_el, d_el)
      print*, "autovetores"
      call print_mat2(phi, d_el, d_el)
      print*, "h_mtx"
      call print_mat2(hMtx, d_el, d_el)
    endif


deallocate(rhomtx) 
end subroutine print_matrices       
        
subroutine monta_rho(vector, rhoMtx)
    implicit none
      !---- subrotina para montar uma matriz no tipo RHO ou RHOPONTO dado um vetor linha Y ou YP ---
      !A MATRIZ RHO FICA    RHO = (y(1) y(4) y(5)) + i (0    y(7) y(8) ) ou seja, temos dimensao*dimensao = 9 variavies
      !COMO (EX 3X3):             (y(4) y(2) y(6))     (y(7)  0   y(9) )          para resolver
      !                           (y(5) y(6) y(3))     (y(8) y(9)  0   )
    
    !args
    real*8,     intent(in)  :: vector(d_el*d_el)
    complex*16, intent(out) :: rhoMtx(d_el,d_el)

    !local 
    real*8                  :: rho_real(d_el,d_el), rho_imag(d_el,d_el)
    integer                 :: ndim  !numero de eq. pra resolver na matriz rho_real
    
    
    rhoMtx = 0.d0 + zi * 0.d0
    ndim = ((d_el*(d_el + 1)) / 2) !numero de eq. pra resolver na matriz rho_real
    rho_real = 0.d0 
    rho_imag = 0.d0
    

    k = d_el + 1
    l = ndim + 1
    do j = 1, d_el-1
      do i = j+1, d_el
        rho_real(i, j) = vector(k) !PARTE NÃO DIAGONAL REAL 
        rho_imag(i, j) = vector(l) !PARTE NÃO DIAGONAL IMAG
        l = l + 1
        k = k + 1
      enddo
    enddo
    
    !==============================
    forall(i = 1:d_el) rho_imag(i, i) = 0.d0 !PARTE DIAGONAL IMAG 
    !=======================================
    
    rhoMtx = rho_real + transpose(rho_real) + zi * ( rho_imag - transpose(rho_imag) )
    
    
    do i = 1, d_el
      rhoMtx(i, i) = vector(i) !PARTE DIAGONAL REAL 
    enddo
    !===============================
    
    !======= VERIFICA SE RHO ESTÁ HERMITIANO ======
    do j = 1, d_el-1
      do i = j+1, d_el
          if ( rhoMtx(i, j) /= conjg(rhoMtx(j, i) )) then
            print*, "RHO NÃO ESTÁ HERMITIANO"
            stop 
          endif
       enddo
    enddo 
    !==============================================

end subroutine monta_rho


!---- subrotina para montar um vetor do tipo Y OU YP  dado uma matriz do tipo RHO OU RHO PONTO ---
! O VETOR Y FICA        y(1) = rho_real(1, 1),  y(4) = rho_real(2, 1), y(7) = rho_imag(2, 1)
! ESCRITO NESSA FORMA   y(2) = rho_real(2, 2),  y(5) = rho_real(3, 1), y(8) = rho_imag(3, 1)
! COMO: (EX 3X3):       y(3) = rho_real(3, 3),  y(6) = rho_real(3, 2), y(9) = rho_imag(3, 2)


subroutine monta_y(rhoMtx, vector)
    implicit none

    !args
    complex*16, intent(in)  :: rhoMtx(d_el, d_el)
    real*8,     intent(out) :: vector(d_el*d_el)
    
    !local
    integer :: ndim  !numero de eq. pra resolver na matriz rho_real
    real*8  :: rho_real(d_el,d_el), rho_imag(d_el, d_el)
    
    

    ndim = ((d_el*(d_el + 1)) / 2) !numero de eq. pra resolver na matriz rho_real

    
    forall(k=1:d_el) vector(k) = real(rhoMtx(k, k)) !PARTE DIAGONAL REAL
    
    k = d_el + 1
    l = ndim + 1
    do j = 1, d_el-1
      do i = j+1, d_el
        vector(k) = real(rhoMtx(i, j)) !PARTE NÃO DIAGONAL REAL
        vector(l) = aimag(rhoMtx(i, j)) !PARTE NÃO DIAGONAL IMAG
        k = k + 1
        l = l + 1
      enddo
    enddo
    
end subroutine monta_y

subroutine rhoMtx_to_pop(rhoMtx, pop_matrix)
    implicit none
!subroutina que devolvolve a populacao do  sitio, ou seja
!matriz_pop(3, 1) = real(rho(3, 3)), matriz_pop(nrows, ncolumns) = real(rho(nr*nc, nr*nc))

    !args
    complex*16, intent(in) :: rhoMtx(d_el, d_el)
    real*8, intent(out)    :: pop_matrix(nr, nc)
    
    !local
    REAL*8 :: soma_temp
    

    pop_matrix = 0.d0
    soma_temp = 0.d0


    k = 1
    do j = 1, nc !-1 !vetor coluna que vai pegar os elementos da diagional principal do rho
      do i = 1, nr
        do l = 1, ns_el
          soma_temp = soma_temp + real(rhoMtx(k, k))
          k = k + 1
        enddo
        pop_matrix(i, j) = soma_temp
        soma_temp = 0.d0
      enddo
    enddo
    
end subroutine rhoMtx_to_pop

subroutine rhosite_TO_rhoham(EGvectors, transpose_EGvectors, rhosite, rhoham)
    implicit none


    !args
    REAL*8,       INTENT(IN) :: EGvectors(d_el, d_el), transpose_EGvectors(d_el, d_el)
    COMPLEX*16,   INTENT(IN) :: rhosite(d_el, d_el)
    COMPLEX*16,   INTENT(OUT) :: rhoham(d_el, d_el)
    
    
    rhoham = matmul( matmul( transpose_EGvectors, rhosite ), EGvectors )   
        
end subroutine rhosite_TO_rhoham


subroutine rhoham_TO_rhosite(EGvectors, transpose_EGvectors, rhoham, rhosite)
    implicit none

    ! args
    REAL*8,     INTENT(IN) :: EGvectors(d_el, d_el), transpose_EGvectors(d_el, d_el)
    COMPLEX*16, INTENT(in) :: rhoham(d_el, d_el)
    COMPLEX*16, INTENT(OUT) :: rhosite(d_el, d_el)
 
    
    rhosite = matmul( matmul( EGvectors, rhoham ), transpose_EGvectors )
    
end subroutine rhoham_TO_rhosite




subroutine printa_resultado(nm_arquivo, t, rhoMtx)
    implicit none
    !SUBROUTINE PARA ESCREVER OS RESULTADOS DO RHO SOMANDO OS ESTADOS PARA CADA SITIO
    !OU SEJA, PARA NSTATES = 3, PROB DE ENCONTRARMOS O ELETRON NO SITIO 1 = RHO(1,1)+RHO(2,2)+RHO(3,3)
    ! args     
    integer,    intent(in) :: nm_arquivo
    real*8,     intent(in) :: t
    complex*16, intent(in) :: rhoMtx(d_el, d_el)
    
    ! local
    real*8, allocatable    :: matriz_temporaria(:)
    real*8                 :: temp_soma
    integer                :: verificador, counter
     

    ALLOCATE(matriz_temporaria(nsites), source = 0.d0) !matriz  temporaria guarda as populacoes de cada sitio
     
     
    counter = 1
    temp_soma = 0.d0
    do  i = 1, nsites  !-1 !vetor coluna que vai pegar os elementos da diagional principal do rho
      do j = 1, ns_el
        temp_soma = real(rhoMtx(counter, counter)) + temp_soma
        counter = counter + 1
      enddo
      verificador = int(temp_soma)
      if (verificador >= 5) then
        print*, "As populacoes estao divergindo! Algo esta errado! Verifique fort.{14, 15, 16, 17}!"
        STOP
      endif
      matriz_temporaria(i) = temp_soma
      temp_soma = 0.d0
    enddo
     
     
     
    write ( nm_arquivo,  '(60F12.5)', advance='no' ) t
    write ( nm_arquivo,  '(60F12.5)', advance = 'no') (  (matriz_temporaria(i)),  i = 1, nsites )
    write ( nm_arquivo, '(60F12.5)' ) sum(matriz_temporaria(:)) 
     
     
    DEALLOCATE(matriz_temporaria)
end subroutine printa_resultado


subroutine open_write_files
implicit none

    open(15,  file =  'popSiteBasis',     status = 'replace')  !eletron base local
    open(14,  file =  'popNonSiteBasis',  status = 'replace')  !eletron base desloc.
    open(100, file =  'radius',           status = 'replace')
    open(101, file =  'Vforces',          status = 'replace')
    open(102, file =  'Eforces',          status = 'replace')
    open(103, file =  'energiatotal',     status = 'replace')
    open(84,  file =  'energiacinetica',  status = 'replace')
    open(85,  file =  'energiapotencial', status = 'replace')
    open(86,  file =  'energiazero',      status = 'replace')
    open(87,  file =  'energiael',        status = 'replace')
    open(80,  file =  'energiaham',       status = 'replace')
    open(105, file =  'temperature',      status = 'replace')


return
end subroutine open_write_files

subroutine close_write_files
implicit none

close(14)
close(15)
close(100)
close(101)
close(102) 
close(103)
close(84)
close(85) 
close(86)
close(87)
close(80) 
close(105) 

return
end subroutine close_write_files


subroutine fderiv(d_el, x_vector, y_vector, deriv)
    implicit none
    !funcao para calcular a derivada de um vetor

    ! args
    integer,   intent(in) :: d_el
    real*8,    intent(in) :: x_vector(d_el), y_vector(d_el)
    real*8,    intent(out) :: deriv(d_el)

    ! local 
    integer   :: i, j
    real*8    :: delx(d_el-1)
    
    delx = 0.d0
    deriv = 0.d0
    
    
    do i = 2, d_el-1
      delx(i) =  x_vector(i) - x_vector(i-1)
      deriv(i) = (y_vector(i+1) - y_vector(i-1)) / (2.d0 * delx(i) )
    enddo
    
    deriv(1) = deriv(2)
    deriv(d_el) = deriv(d_el-1)
end subroutine fderiv


subroutine print_mat2(aa, nn, mm)
implicit none

    ! args
    integer,   intent(in) :: nn, mm
    real*8,    intent(in) :: aa(nn, mm)

    !local
    integer :: i, j

    do i=1,mm
      write(*,'(100g12.4)') ( aa(i,j), j=1,nn )
    enddo
end subroutine print_mat2





!funcao para fazer a integral numerica de um vetor
!===================================
function sumtrap(i1,i2,eixox,eixoy)
!===================================
implicit none 

    ! args 
    integer,  intent(in) :: i1 , i2
    real*8,   intent(in) :: eixox(:)
    real*8,   intent(in) :: eixoy(:)
    
    ! local 
    real*8  :: sumtrap

!------------------------------------------------------------------------------
! CALCULA A INTEGRAL DA FUNCAO Y(I) PELO METODO DO TRAPEZIO COM PASSO VARIAVEL
!------------------------------------------------------------------------------

     sumtrap  = sum( (eixox(i1+1:i2)-eixox(i1:i2-1)) * (eixoy(i1+1:i2)+eixoy(i1:i2-1)) ) / 2.0d0

end function sumtrap


subroutine printaletters2(nomearquivo, nmfile, matrizsize, nmlines, nmcol, matriz)
    implicit none

    ! args 
    character(len=6), intent(in) :: nomearquivo
    integer,          intent(in) :: nmfile, matrizsize, nmlines, nmcol
    real*8,           intent(in) :: matriz(matrizsize, matrizsize)
    
    
    open(file = nomearquivo, status = "replace", unit = nmfile)
    do j = 1, nmcol
      do i = 1, nmlines
        write(nmfile, 13, advance = "no") matriz(i, j)
      enddo
      write(nmfile, 13, advance = "yes")
    enddo
    
    close(nmfile)
    
    
    13 format (10F8.3)
end subroutine printaletters2


function lm_polynomial(n, l, x) 
    implicit none

    !args
    integer, intent(in) :: n, l
    real*8, intent(in) :: x
    real*8 :: lm_polynomial

    !local
    real*8 :: lr

    lr = float(l)

    
    SELECT CASE (n)
      CASE(0)
        lm_polynomial = 1.d0
      CASE(1)
        lm_polynomial = -x + lr + 1.d0
      CASE(2)
        lm_polynomial = HALF * (x * x - 2.d0 * x * ( lr + 2.d0 ) + ( (lr + 1.d0) * (lr + 2.d0) ) ) 
      CASE(3)
        lm_polynomial = 0.16666d0 * ( - x**3.0 + (3.d0 * x * x * ( lr + 3.d0 )) - (3.d0 * x * ( lr + 2.d0 ) * ( lr + 3.d0 ))  + &
                                      ( ( lr + 1.d0 ) * ( lr + 2.d0 ) * ( lr + 3.d0 ) ) ) 


    END SELECT

end function lm_polynomial 


subroutine plot_wavefunction(nn, ll, points) 
    implicit none

    !args
    integer, intent(in) :: nn, ll, points

    !local
    real*8 :: rr(points), theta(points)
    real*8 :: factor, norm, lr, fnorm1, fnorm2, abslr
    complex*16, allocatable, dimension(:,:) :: phipolar, phipolarsq

    allocate(phipolar(points, points), source = (0.d0, 0.d0)) 
    allocate(phipolarsq(points, points), source = (0.d0, 0.d0)) 

    
    lr = float(ll)
    abslr = abs(lr)
    factor = (me * site(1, 1)%omega * 1.d12) / hbar
    fnorm1 = factor**( HALF * (abslr + 1.d0))
    fnorm2 = sqrt( Gamma(float(nn)+1.d0) / ( PI * Gamma(float(nn) + abslr + 1.d0)  ) )
    norm = fnorm1 * fnorm2


    print*, factor, "factor"
    print*, fnorm1, "fnorm1"
    print*, fnorm2, "fnorm2"
    print*, norm, "norm" 


    forall(i = 1 : points) rr(i) = (( float(i-1) * 2.d0 ) / float(points-1))
    forall(i = 1 : points) theta(i) = ((float(i-1) * 2.d0 * PI ) / float(points-1))



    rr = rr * 1.d-9 


    do j = 1, points
      do i = 1, points
        phipolar(i, j) = norm * (cos(lr * theta(j)) + zi * sin(lr * theta(j))) * (rr(i)**abslr) * exp(- HALF * factor * rr(i)**2.0) * &
                         lm_polynomial(nn, ll, factor * rr(i)**2.0)
      enddo
    enddo

    phipolarsq = phipolar * conjg(phipolar) 

    open(200, file = 'phipolar', status = 'replace' ) 
    open(201, file = 'radialfunc', status = 'replace' )


    do j = 1, points
      do i = 1, points
      write(200, '(60F12.5)') rr(i) * 1.d9 * cos(theta(j)) , rr(i) * 1.d9 * sin(theta(j)) , 4.d0 * PI * rr(i)**2.0 * real(phipolarsq(i, j))
      enddo
      write(200, '(60F12.5)') 
    enddo


    do i = 1, points
      write(201, '(60F12.5)') -rr(i) * 1.d9, 4.d0 * PI * rr(i)**2.d0 * real(phipolarsq(i, points/2))
    enddo

    do i = 1, points
      write(201, '(60F12.5)') rr(i) * 1.d9, 4.d0 * PI * rr(i)**2.d0 * real(phipolarsq(i, 1))
    enddo

    close(200) 
    close(201)

    deallocate(phipolarsq) 

end subroutine plot_wavefunction  



end module functions_m
