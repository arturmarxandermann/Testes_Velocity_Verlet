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
  integer, intent(in) :: nn, mm
  real*8, dimension(nn, mm), intent(in) :: aa
  integer :: i, j
  do, i=1,mm
      write(*,'(100g12.4)') ( aa(i,j), j=1,nn )
enddo

end subroutine print_mat3


!subrotina que calcula a hamiltoniana para eletron ou buraco.
!se el_or_hl = 1 a ham. output é para eletron
!se el_or_hl = 2 a ham. output é para buraco
subroutine build_hamiltonian(el_or_hl, pop_el, pop_hl, ti, hamiltonian)
implicit none
integer, INTENT(IN) :: el_or_hl
REAL*8, DIMENSION(nm_rows, nm_columns), INTENT(IN) :: pop_el, pop_hl
real*8, intent(in) :: ti 
real*8, INTENT(OUT), allocatable, DIMENSION(:, :) :: hamiltonian

integer :: nmstates, c1, r1, s_c1, s_r1, s1, s2, c2, r2
integer :: dims
real*8 :: mult_fator, V1, V2, w1, w2, dist
real*8, allocatable, DIMENSION(:, :) :: S
real*8, ALLOCATABLE, DIMENSION(:, :) :: coup_temp
real*8, DIMENSION(nm_rows, nm_columns) :: site_type, omegas
real*8, DIMENSION(:, :), ALLOCATABLE :: energy_factor, elhl_coupling

real*8, DIMENSION(6) :: Qnumbers_energy
INTEGER, DIMENSION(6) :: Qnumbers_row
real*8, dimension(:, :), allocatable :: pop_fixada, pop_naofixada, ksi
real*8, dimension(:), allocatable :: ksi_sum
real*8 :: ksi_summ, CteSum 
real*8, dimension(:), allocatable :: elhl_plot

!criando o vetor associado aos valores dos numeros quanticos dos estados do osciladores harmonicos
Qnumbers_row = (/ 00, 01, 10, 02, 11, 20 /)

!criando o vetor com as energias associado a cada numero quantico
Qnumbers_energy = (/ 1.d0, 2.d0, 2.d0, 3.d0, 3.d0, 3.d0 /)

if (el_or_hl /= 1 .AND.  el_or_hl /= 2 ) then
  print*, "Número correspondente ao elétron ou buraco errado."
  print*, "1 => elétron, 2=> buraco"
  STOP
end if

allocate(pop_fixada(nm_rows, nm_columns), source = 0.d0)
allocate(pop_naofixada(nm_rows, nm_columns), source = 0.d0)

if (el_or_hl == 1 ) then
  mult_fator = 1.d0
  nmstates = nmstates_el
  dims = dim_el
  site_type(:, :) = sites_array(:, :)%site_type_el
  omegas(:, :) = sites_array(:, :)%omega
  allocate( hamiltonian(dim_el , dim_el), source = 0.d0 )
  allocate( coup_temp(nmstates, nmstates), source = 0.d0 )
  pop_fixada = pop_el
  pop_naofixada = pop_hl
end if

if (el_or_hl == 2 ) then
  mult_fator = 1.d0
  nmstates = nmstates_hl
  dims = dim_hl
  site_type(:, :) = sites_array(:, :)%site_type_hl
  omegas(:, :) = sites_array(:, :)%omega
  allocate( hamiltonian(dim_hl , dim_hl), source = 0.d0 )
  allocate( coup_temp(nmstates, nmstates), source = 0.d0 )
  pop_fixada = pop_hl
  pop_naofixada = pop_el
end if

ALLOCATE(energy_factor(nm_rows, nm_columns), source = 0.d0)
ALLOCATE(elhl_coupling(nm_rows, nm_columns), source = 0.d0)
allocate(ksi(nm_rows, nm_columns), source = 0.d0)
allocate(ksi_sum(nm_columns), source = 0.d0)
allocate(elhl_plot(nsites), source = 0.d0 )
allocate(S(nmstates, nmstates), source = 0.d0 ) 

do c1 = 1, nm_columns
    do r1 = 1, nm_rows  
        do c2 = 1, nm_columns
            do r2 = 1, nm_rows 
               if (c1 /= c2 .OR. r1 /= r2) then
                    dist = sqrt( (sites_array(r1, c1)%posicao_x - sites_array(r2, c2)%posicao_x)**2.0 + &
                    (sites_array(r1, c1)%posicao_y - sites_array(r2, c2)%posicao_y)**2.0 ) * 1.d9
                    ksi(r2, c2) =  pop_naofixada(r2, c2) / dist 
               else
                    ksi(r2, c2) = pop_naofixada(r2, c2) 
               endif 
            enddo
            ksi_sum(c2) = sum(ksi(:, c2))
        enddo
        ksi_summ = sum(ksi_sum(:))
        elhl_coupling(r1, c1) = elhlCouplingCte * sites_array(r1, c1)%Egap * pop_fixada(r1, c1) * ksi_summ
    enddo   
enddo



k = 1
do j = 1, nm_columns
  do i = 1, nm_rows
    do s1 = 1, nmstates
       hamiltonian(k, k) = site_type(i, j) + HB_ev_ps * omegas(i, j) *  Qnumbers_energy(s1) !- elhl_coupling(i, j) 
       k = k + 1
     enddo
   enddo
enddo






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
                                                                  RESHAPE(coup_temp, (/nmstates, nmstates/) )
          hamiltonian(l*nmstates+1:l*nmstates+nmstates, k*nmstates+1:k*nmstates+nmstates) = &
                                                                  RESHAPE(transpose(coup_temp), (/nmstates, nmstates/) )
        endif
      k = k + 1
      enddo
    enddo
    l = l + 1
  enddo
enddo

DEALLOCATE(energy_factor, elhl_coupling, pop_fixada, pop_naofixada, ksi, ksi_sum, elhl_plot, S)
return

end subroutine build_hamiltonian

subroutine calculate_eigenvectors(ham_size, hamiltoniana, energias, phi, phi_transpose, omega_matrix)
!calcula os autovetores (phi) e autovalores (energias) da hamiltoniana. omega_matriz é a matriz de frequencias angulares
implicit none
integer, intent(in) :: ham_size
real*8, intent(in)  :: hamiltoniana(ham_size, ham_size)
real*8, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: energias
real*8, allocatable, DIMENSION(:, :), intent(out) :: phi, phi_transpose, omega_matrix

real*8, dimension(:), allocatable :: ovlp_eigenvalues
real*8, dimension(:, :), ALLOCATABLE ::  hamiltoniana_diag, ovlp_diag  !matriz hamiltoniana que vai ser diagonalizada e virar os autovetores


character*1 :: jobz, uplo
integer :: info, is, js 
character(len=6) :: file_name
character(len=2) :: file_id


allocate(phi(ham_size, ham_size), source = 0.d0)
allocate(phi_transpose(ham_size, ham_size), source = 0.d0)
allocate(omega_matrix(ham_size, ham_size), source = 0.d0)
allocate(energias(ham_size), source = 0.d0)
allocate(hamiltoniana_diag(ham_size, ham_size), source = 0.d0 )

jobz = "V"
uplo = "L" 
info = 0
hamiltoniana_diag = hamiltoniana




call syevd(hamiltoniana_diag, energias, jobz, uplo, info)

!If info = 0, contains the n orthonormal eigenvectors, stored by columns.
!(The i-th column corresponds to the ith eigenvalue.)
!Como o iézimo autovetor é disposto na iézima coluna, C_n^N = <n|N> é o coeficiente da expansão do autovetor

if ( info /= 0 ) write(*,*) info, "Error to found eigenvalues"


do i = 1, ham_size
print*, "ENERGIA", i, energias(i)
enddo


do n = 1, ham_size
  phi(n,:) = hamiltoniana_diag(n,:)
enddo



phi_transpose = transpose(phi)




!DIR$ PARALLEL
do j = 1, ham_size
  do i = 1, ham_size
    omega_matrix(i, j) = (energias(i) - energias(j)) / HB_ev_ps !omegas em 1/s
  enddo
enddo




13 format(3es14.5E3)


deallocate(hamiltoniana_diag) 
return

end subroutine calculate_eigenvectors


subroutine build_overlap_matrix(nmstates, hamsize, ovlpmatrix)
implicit none

integer, intent(in) :: nmstates, hamsize
real*8, ALLOCATABLE, DIMENSION(:, :), INTENT(out) :: ovlpmatrix

real*8, allocatable, DIMENSION(:, :) :: S
integer :: c1, c2, r1, r2
INTEGER, DIMENSION(6) :: Qnumbers_row
ALLOCATE(ovlpmatrix(hamsize, hamsize), source = 0.d0)


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


do j = 1, hamsize
    ovlpmatrix(j, j) = 1.d0 
enddo 

return
end subroutine build_overlap_matrix


subroutine calculate_derivative_matrix(nmstates, hamsize, dermtx)
implicit none
integer, intent(in) :: nmstates, hamsize
real*8, allocatable, dimension(:, :), intent(out) :: dermtx
real*8, allocatable, dimension(:, :) :: derivativeMtx
integer :: c1, r1, c2, r2

allocate(dermtx(hamsize, hamsize), source = 0.d0 )


l = 0
k = 0
do c1 = 1, nm_columns
  do r1 = 1, nm_rows
      k = 0
    do c2 = 1, nm_columns
      do r2 = 1, nm_rows
        if ( k > l ) then
        call DerivativeOverlap(sites_array, nmstates_el, r2, c2, r1, c1, derivativeMtx)
        dermtx(k*nmstates+1:k*nmstates+nmstates, l*nmstates+1:l*nmstates+nmstates) =  &
                                                                RESHAPE(derivativeMtx, (/nmstates, nmstates/) )

        dermtx(l*nmstates+1:l*nmstates+nmstates, k*nmstates+1:k*nmstates+nmstates) =  &
                                                                RESHAPE(transpose(derivativeMtx), (/nmstates, nmstates/) )
                                
                                        
                                        
        endif 


        k = k + 1
      enddo
    enddo
      l = l + 1
  enddo
enddo







return
end subroutine calculate_derivative_matrix


subroutine Newbuild_overlap_matrix(nmstates, hamsize, ovlpmatrix)
implicit none

integer, intent(in) :: nmstates, hamsize
real*8, ALLOCATABLE, DIMENSION(:, :), INTENT(out) :: ovlpmatrix

real*8, allocatable, DIMENSION(:, :) :: S
integer :: c1, c2, r1, r2
INTEGER, DIMENSION(6) :: Qnumbers_row
ALLOCATE(ovlpmatrix(hamsize, hamsize), source = 0.d0)


Qnumbers_row = (/ 00, 01, 10, 02, 11, 20 /)



l = 0
k = 0
do c1 = 1, nm_columns
  do r1 = 1, nm_rows
      k = 0
    do c2 = 1, nm_columns
      do r2 = 1, nm_rows
        if ( k > l ) then
        call NewOverlap(sites_array, nmstates, r2, c2, r1, c1, Qnumbers_row, S)
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


do j = 1, hamsize
    ovlpmatrix(j, j) = 1.d0
enddo

return
end subroutine Newbuild_overlap_matrix

subroutine printaletters3(nomearquivo, nmfile, matrizsize, nmlines, nmcol, matriz)
implicit none
character*6, INTENT(IN) :: nomearquivo
integer, intent(in) :: nmfile, matrizsize, nmlines, nmcol
real*8, intent(in), dimension(matrizsize, matrizsize) :: matriz


open(file = nomearquivo, status = "replace", unit = nmfile)
do i = 1, nmlines
  do j = 1, nmcol
    write(nmfile, 13, advance = "no") matriz(i, j)
  enddo
  write(nmfile, 13, advance = "yes")
enddo

close(nmfile)


13 format (8F6.2)
return
end subroutine printaletters3





end module system_hamiltonian_m
