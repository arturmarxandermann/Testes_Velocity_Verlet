module gnugraphs_m
use parameters_m
use constants_m
use overlap_m
use types_m


contains
subroutine newcalcula_cor(charge_site, pop_el_matrix, pop_hl_matrix,  colorintensity, hue_color) 
real*8, intent(in) :: charge_site, pop_el_matrix, pop_hl_matrix
real*8, intent(out) :: colorintensity, hue_color
real*8 :: zeropop, lim_exciton


zeropop = 0.01d0
lim_exciton = 0.05d0

colorintensity = 0.0
hue_color = 0.0
    


 !FAZENDO A PARTE DO ELETRON 
 
  if ( charge_site >= zeropop .AND. pop_el_matrix >= zeropop) then !.AND. pop_hl_matrix <= zeropop ) then
    hue_color = 0.3d0

    if ( pop_el_matrix >= 0.66d0 ) then
    colorintensity = 1.0d0
    endif

    if ( pop_el_matrix < 0.66d0 .AND. pop_el_matrix >= 0.33d0) then
    colorintensity = 0.5d0
    endif

    if ( pop_el_matrix < 0.33d0 .AND. pop_el_matrix >= 0.2d0 ) then
    colorintensity = 0.4d0
    endif

    if ( pop_el_matrix < 0.2d0 .AND. pop_el_matrix >= 0.1d0 ) then
    colorintensity = 0.35d0
    endif

    if ( pop_el_matrix < 0.1d0 .AND. pop_el_matrix >= 0.03d0 ) then
    colorintensity = 0.3d0
    endif
    
    if ( pop_el_matrix < 0.03d0 .AND. pop_el_matrix >= 0.02d0 ) then
    colorintensity = 0.25d0
    endif

    if ( pop_el_matrix < 0.02d0 .AND. pop_el_matrix >= zeropop ) then
    colorintensity = 0.2d0
    endif

    if ( pop_el_matrix < 0.01d0 .AND. pop_el_matrix >= 0.005d0 ) then
    colorintensity = 0.1d0
    endif

    if ( pop_el_matrix < 0.005d0 ) then
    colorintensity = 0.d0
    endif


  endif

   !if ( charge_site > 0.007d0 .AND. charge_site <= 0.02d0 ) then
   !     colorintensity = 0.25d0
   !     hue_color = 0.3
   !endif 


   ! if ( charge_site > 0.02d0 .AND. charge_site <= 0.025d0 ) then
   !      colorintensity = 0.5d0 
   !      hue_color = 0.3  
   ! endif


   ! if ( charge_site > 0.025d0 .AND. charge_site <= 0.03d0 ) then
   !      colorintensity = 0.55d0
   !      hue_color = 0.3  
   ! endif
   ! 
   ! if ( charge_site > 0.03d0 .AND. charge_site <= 0.035d0 ) then
   !      colorintensity = 0.6d0
   !      hue_color = 0.3   
   ! endif

   !if ( charge_site > 0.035d0 .AND. charge_site <= 0.04d0 ) then
   !      colorintensity = 0.65d0
   !      hue_color = 0.3   
   !endif

   !if ( charge_site > 0.04d0 .AND. charge_site <= 0.045d0 ) then
   !     colorintensity = 0.7d0
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.045d0 .AND. charge_site <= 0.05d0 ) then
   !     colorintensity = 0.73d0
   !     hue_color = 0.3
   !endif
  
   !if ( charge_site > 0.05d0 .AND. charge_site <= 0.06d0 ) then 
   !     colorintensity = 0.76d0
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.06d0 .AND. charge_site <= 0.07d0 ) then 
   !     colorintensity = 0.80
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.07d0 .AND. charge_site <= 0.08d0 ) then 
   !     colorintensity = 0.80
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.08d0 .AND. charge_site <= 0.09d0 ) then 
   !     colorintensity = 0.80
   !     hue_color = 0.38
   !endif

   !if ( charge_site > 0.09d0 .AND. charge_site <= 0.1d0 ) then 
   !     colorintensity = 0.80
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.1d0 .AND. charge_site <= 0.2d0 ) then 
   !     colorintensity = 0.80
   !     hue_color = 0.3
   !endif
   !
  !if ( charge_site > 0.2d0 .AND. charge_site <= 0.25d0 ) then 
   !     colorintensity = 0.83
   !     hue_color = 0.3
   !endif

  !if ( charge_site > 0.25d0 .AND. charge_site <= 0.3d0 ) then 
   !     colorintensity = 0.86
   !     hue_color = 0.3
   !endif

   !if ( charge_site > 0.3d0 .AND. charge_site <= 0.34d0 ) then 
   !     colorintensity = 0.90
   !     hue_color = 0.3
   !endif
 
   !if ( charge_site > 0.34d0 .AND. charge_site <= 0.36d0 ) then 
   !     colorintensity = 0.94
   !     hue_color = 0.3
   !endif
 
   !if ( charge_site > 0.36d0 .AND. charge_site <= 0.4d0 ) then 
   !     colorintensity = 0.98
   !     hue_color = 0.3
   !endif
  
   !if ( charge_site > 0.4d0 .AND. charge_site <= 1.0d0 ) then 
   !     colorintensity = 1.0
   !     hue_color = 0.3
   !endif
 

!FAZENDO A PARTE DO BURACO


  if ( charge_site <= - zeropop ) then
    hue_color = 0.d0 

    if ( pop_hl_matrix >= 0.66d0 ) then
    colorintensity = 1.0d0
    endif

    if ( pop_hl_matrix < 0.66d0 .AND. pop_hl_matrix >= 0.33d0) then
    colorintensity = 0.5d0
    endif

    if ( pop_hl_matrix < 0.33d0 .AND. pop_hl_matrix >= 0.2d0 ) then
    colorintensity = 0.4d0
    endif

    if ( pop_hl_matrix < 0.2d0 .AND. pop_hl_matrix >= 0.1d0 ) then
    colorintensity = 0.35d0
    endif

    if ( pop_hl_matrix < 0.1d0 .AND. pop_hl_matrix >= 0.03d0 ) then
    colorintensity = 0.3d0
    endif

    if ( pop_hl_matrix < 0.03d0 .AND. pop_hl_matrix >= 0.02d0 ) then
    colorintensity = 0.25d0
    endif
    
    if ( pop_hl_matrix < 0.02d0 .AND. pop_hl_matrix >= zeropop ) then
    colorintensity = 0.2d0
    endif


    if ( pop_hl_matrix < zeropop ) then
    colorintensity = 0.d0
    endif

  !  if ( pop_hl_matrix >= 0.66d0 ) then
  !  colorintensity = 1.d0
  !  endif

  !  if ( pop_hl_matrix < 0.66d0 .AND. pop_hl_matrix >= 0.33d0) then
  !  colorintensity = 0.5d0
  !  endif

  !  if ( pop_hl_matrix < 0.33d0 .AND. pop_hl_matrix >= zeropop ) then
  !  colorintensity = 0.3d0
  !  endif

  !  if ( pop_hl_matrix < zeropop ) then
  !  colorintensity = 0.d0
  !  endif


  endif


  ! if ( charge_site >= -0.05d0 .AND. charge_site < -0.03d0 ) then 
  !     colorintensity = 0.55
  !     hue_color = 0.d0 
  ! endif
  ! 
  ! if ( charge_site >= -0.06d0 .AND. charge_site < -0.04d0 ) then
  !     colorintensity = 0.6
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -0.07d0 .AND. charge_site < -0.06d0 ) then
  !     colorintensity = 0.65
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -0.1d0 .AND. charge_site < -0.07d0 ) then
  !     colorintensity = 0.7
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -0.2d0 .AND. charge_site < -0.1d0 ) then
  !     colorintensity = 0.75
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -0.3d0 .AND. charge_site < -0.2d0 ) then
  !     colorintensity = 0.8
  !     hue_color = 0.d0
  ! endif
  !
  ! if ( charge_site >= -0.4d0 .AND. charge_site < -0.3d0 ) then
  !     colorintensity = 0.85
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -0.5d0 .AND. charge_site < -0.4d0 ) then
  !     colorintensity = 0.9
  !     hue_color = 0.d0
  ! endif
 
  ! if ( charge_site >= -0.7d0 .AND. charge_site < -0.5d0 ) then
  !     colorintensity = 0.93
  !     hue_color = 0.d0
  ! endif
  ! 
  ! if ( charge_site >= -0.8d0 .AND. charge_site < -0.7d0 ) then
  !     colorintensity = 0.96
  !     hue_color = 0.d0
  ! endif

  ! if ( charge_site >= -1.0d0 .AND. charge_site < -0.8d0 ) then
  !     colorintensity = 1.0
  !     hue_color = 0.d0
  ! endif


  !FAZENDO PARTE DO EXCITON 


  if ( abs(charge_site) < lim_exciton .AND. pop_el_matrix > zeropop .AND. pop_hl_matrix > zeropop ) then
  
  hue_color = 0.15d0

  if ( pop_el_matrix >= 0.66d0 ) then
  colorintensity = 1.d0
  endif

  if ( pop_el_matrix >= 0.33d0 .AND. pop_el_matrix < 0.66d0 ) then
  colorintensity = 0.55d0
  endif


  if ( pop_el_matrix >= 0.2d0 .AND. pop_el_matrix < 0.33d0 ) then
  colorintensity = 0.45d0
  endif


  if ( pop_el_matrix >= 0.02d0 .AND. pop_el_matrix < 0.2d0 ) then
  colorintensity = 0.35d0
  endif


  if ( pop_el_matrix >= zeropop .AND. pop_el_matrix < 0.02d0 ) then
  colorintensity = 0.25d0
  endif



  endif








!  if ( pop_el_matrix > 0.02d0 .AND. pop_hl_matrix > 0.02d0 .AND. &    
!  pop_el_matrix <= 0.1d0 .AND. pop_hl_matrix <= 0.1d0 ) then
!     colorintensity = 0.3d0
!     hue_color = 0.15
!  endif
! 
!  if ( pop_el_matrix > 0.1d0 .AND. pop_hl_matrix > 0.1d0 .AND. &    
!  pop_el_matrix <= 0.2d0 .AND. pop_hl_matrix <= 0.2d0 ) then
!     colorintensity = 0.4d0
!     hue_color = 0.15
!  endif
! 
!  if ( pop_el_matrix > 0.2d0 .AND. pop_hl_matrix > 0.2d0 .AND. &    
!  pop_el_matrix <= 0.3d0 .AND. pop_hl_matrix <= 0.3d0 ) then
!     colorintensity = 0.5d0
!     hue_color = 0.15
!  endif
! 
!  if ( pop_el_matrix > 0.3d0 .AND. pop_hl_matrix > 0.3d0 .AND. &
!  pop_el_matrix <= 0.4d0 .AND. pop_hl_matrix <= 0.4d0 ) then
!     colorintensity = 0.6d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.4d0 .AND. pop_hl_matrix > 0.4d0 .AND. &
!  pop_el_matrix <= 0.5d0 .AND. pop_hl_matrix <= 0.5d0 ) then
!     colorintensity = 0.75d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.5d0 .AND. pop_hl_matrix > 0.5d0 .AND. &
!  pop_el_matrix <= 0.6d0 .AND. pop_hl_matrix <= 0.6d0 ) then
!     colorintensity = 0.8d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.6d0 .AND. pop_hl_matrix > 0.6d0 .AND. &
!  pop_el_matrix <= 0.7d0 .AND. pop_hl_matrix <= 0.7d0 ) then
!     colorintensity = 0.85d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.7d0 .AND. pop_hl_matrix > 0.7d0 .AND. &
!  pop_el_matrix <= 0.8d0 .AND. pop_hl_matrix <= 0.8d0 ) then
!     colorintensity = 0.9d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.8d0 .AND. pop_hl_matrix > 0.8d0 .AND. &
!  pop_el_matrix <= 0.9d0 .AND. pop_hl_matrix <= 0.9d0 ) then
!     colorintensity = 0.95d0
!     hue_color = 0.15
!  endif
!
!  if ( pop_el_matrix > 0.9d0 .AND. pop_hl_matrix > 0.9d0 .AND. &
!  pop_el_matrix <= 1.d0 .AND. pop_hl_matrix <= 1.d0 ) then
!     colorintensity = 1.d0
!     hue_color = 0.15
!  endif



return
end subroutine newcalcula_cor



subroutine calcula_cor(pop_electron, pop_buraco, charge_site, hue_color, color_intensity) 
implicit none
real*8, intent(in) :: pop_electron, pop_buraco
real*8, intent(in) :: charge_site
real*8, intent(out) :: hue_color, color_intensity


real*8, parameter :: lim_exciton = 0.02d0 !(limite para definir se é um exciton ou cargas)


if ( abs(charge_site) <= lim_exciton ) then
    hue_color = 0.15d0
    if ( pop_electron >= pop_buraco ) color_intensity = pop_electron
    if ( pop_electron < pop_buraco ) color_intensity = pop_buraco
endif

if ( charge_site > lim_exciton ) then 
    hue_color = 0.3d0
    color_intensity = charge_site
endif

if ( charge_site < - lim_exciton ) then
    hue_color = 0.d0
    color_intensity = abs(charge_site) 
endif

color_intensity = 5.d0 * color_intensity
if (color_intensity > 1.d0 ) color_intensity = 1.d0 

return

end subroutine calcula_cor




subroutine GnuFiles(sites_array, pop_el_matrix, pop_hl_matrix, figure_label, time)
implicit none
type(quantum_site), intent(inout) :: sites_array(nm_rows, nm_columns)
real*8, intent(in) :: pop_el_matrix(nm_rows, nm_columns)
real*8, intent(in) :: pop_hl_matrix(nm_rows, nm_columns)
integer, intent(in) :: figure_label
real*8, intent(in) :: time

real*8, parameter :: valor_min_cor = 0.01d0
real*8, parameter :: scale_nanometro = 1.d-9
real*8, parameter :: scale_freq = 1.d13
integer, parameter :: alpha_max = 1
integer :: p
real*8 :: charge_site


real*8 :: hue_color, colorintensity
real*8 :: time_position_x, time_position_y
real*8, parameter :: clr = 1.d0
!xmin = (sites_array(1, 1)%posicao_x - sites_array(1, 1)%radius)/scale_nanometro
!xmax = (sites_array(1, nm_columns)%posicao_x + sites_array(1, nm_columns)%radius)/scale_nanometro
!ymin = (sites_array(1, 1)%posicao_y - sites_array(1, 1)%radius)/scale_nanometro
!ymax = (sites_array(nm_rows,1)%posicao_y + sites_array(nm_rows, 1)%radius)/scale_nanometro
!XMIN, XMAX, YMIN, YMAX DEFINIDOS EM OVERLAP. DEFINO DIRETO EM OVERLAP PARA PODER USAR NA SUBROTINA PLOT_AUTO_FUNCTION

time_position_x = sites_array(1, 1)%posicao_x/scale_nanometro
time_position_y = ymax - 2.d0


!write(20, '(a)') '#!/usr/bin/gnuplot'

write(20, '(a)') 'clear'
write(20, '(a)') 'reset'
write(20, '(a)') "set terminal pngcairo dashed size 1080, 700 enhanced font 'Verdana, 20'"
write(20, '(a)') 'set label'
write(20, '(a)') 'set key'
write(20, '(a)') "set encoding iso_8859_1"
write(20, '(a)') "set size ratio -1" !deixa tudo na mesma escala
write(20, '(a)') "unset cbrange"
write(20, '(a)') "unset cbtics"
write(20, '(a)') "unset colorbox"


write(20, '(a, F8.3, a, F8.3, a)') 'set xrange [',xmin,':',xmax + 10.0,']'
write(20, '(a, F8.3, a, F8.3, a)') 'set yrange [',ymin,':',ymax,']'
write(20, '(a)') 'set xlabel "x(nm)" '
write(20, '(a)') 'set ylabel "y(nm)" '

!write(20, '(a)') 'set palette model HSV defined (-1.0 0 1 1, -0.1002 0 0.15 1, -0.1001 0.2 0.15 1, -0.001 0.2 1 1,  0.001 0.2 1 1, &
 !                  0.1001 0.2 0.15 1, 0.1002 0.6 0.15 1,  1 0.6 1 1)'
!Paleta de cores esta no modo HSV e a numeracao indica: (ponto_inicial_paleta matiz(hue) saturação brilho,
!                                                          ponto_final_paleta  matiz saturação brilho)


!write(20, '(a)') 'set palette model HSV defined (-1.0 0 1 1, 1 0.3 1 1)'
!write(20, '(a)') 'set palette model HSV defined (-1.0 0 1 1, 1 0.3 1 1)'
write(20, '(a)') 'set style line 1 lw 2 dt 2 lc rgb "black"'
write(20, '(a)') 'set style line 2 lw 2 dt 3 lc rgb "black"'
write(20, '(a)') 'set style line 3 lw 2 dt 4 lc rgb "black"'
write(20, '(a)') 'set style line 4 lw 2 dt 5 lc rgb "black"'


!write(20, '(a)') 'set cbrange[-1:1]' !cbrange indica as propriedades da palheta de cores
!write(20, '(a)') 'set cbtics 0.2'
!write(20, '(a)') 'set cblabel " {/Symbol r}_{el}  - {/Symbol r}_{hl}" '


write(20, '(a, F10.4, a, F8.3, a, F8.3, a)') 'set label 1 sprintf("TIME = %4.4f ps",  ',time,') at ',time_position_x,',',time_position_y,' font "arialbd, 22"'


write(20,fmt='(a)')        "system('mkdir -p animacao')"
write(20,fmt='(a, I5, a)') "fname = sprintf('animacao/animation_%04d.png',",figure_label,")"
write(20, '(a)') 'set output fname'

!ANODOS : DASHTYPE 2
!DOADORES : DASHTYPE 3
!ACEITADORES : DASHTYPE 4
!CATODOS : DASHTYPE 5


k = 1
do c = 1, nm_columns
  do r = 1, nm_rows
    charge_site = pop_el_matrix(r, c) - pop_hl_matrix(r, c)
    !colorintensity = 0.d0
    !hue_color = 0.d0  

    !call calcula_cor( pop_el_matrix(r, c), pop_hl_matrix(r, c), charge_site, hue_color, colorintensity )   
    call newcalcula_cor( charge_site, pop_el_matrix(r, c), pop_hl_matrix(r, c), colorintensity, hue_color )   



if (c >= anode_layer_b .AND. c <= anode_layer_e ) then !====== FAZENDO O CIRCULO PARA O ÂNODO: DASHTYPE 2 ==================
      write(20,fmt='(a, I3, a, F8.3, a, F8.3, a, F7.3, a, F8.3, a, F8.3, a )') &
'set object ',k,' circle center',sites_array(r, c)%posicao_x/scale_nanometro,',',sites_array(r, c)%posicao_y/scale_nanometro,' size &
',sites_array(r, c)%radius/scale_nanometro , 'fc rgb hsv2rgb(',hue_color,',', colorintensity,', 1.0) &
fs solid 3.0 border rgb "black" lw 2.0 dashtype 2'
endif
if (c > anode_layer_e .AND. c <= donnor_layer_e) then !======= FAZENDO O CIRCULO PARA O MATERIAL DOADOR: DASHTYPE 3  =================
      write(20,fmt='(a, I3, a, F8.3, a, F8.3, a, F7.3, a, F8.3, a, F8.3, a )') &
'set object ',k,' circle center',sites_array(r, c)%posicao_x/scale_nanometro,',',sites_array(r, c)%posicao_y/scale_nanometro,' size &
',sites_array(r, c)%radius/scale_nanometro , 'fc rgb hsv2rgb(',hue_color,',', colorintensity,', 1.0) &
fs solid 3.0 border rgb "black" lw 2.0 dashtype 3'
endif
if (c > donnor_layer_e .AND. c <= acceptor_layer_e) then !======= FAZENDO O CIRCULO PARA O MATERIAL ACEITADOR: DASHTYPE 4  =================
      write(20,fmt='(a, I3, a, F8.3, a, F8.3, a, F7.3, a, F8.3, a, F8.3, a )') &
'set object ',k,' circle center',sites_array(r, c)%posicao_x/scale_nanometro,',',sites_array(r, c)%posicao_y/scale_nanometro,' size &
',sites_array(r, c)%radius/scale_nanometro , 'fc rgb hsv2rgb(',hue_color,',', colorintensity,', 1.0) &
fs solid 3.0 border rgb "black" lw 2.0 dashtype 4'
endif
if (c > acceptor_layer_e) then !======= FAZENDO O CIRCULO PARA O CATODO: DASHTYPE 5  =================
      write(20,fmt='(a, I3, a, F8.3, a, F8.3, a, F7.3, a, F8.3, a, F8.3, a )') &
'set object ',k,' circle center',sites_array(r, c)%posicao_x/scale_nanometro,',',sites_array(r, c)%posicao_y/scale_nanometro,' size &
',sites_array(r, c)%radius/scale_nanometro , 'fc rgb hsv2rgb(',hue_color,',', colorintensity,', 1.0) &
fs solid 3.0 border rgb "black" lw 2.0 dashtype 5'
endif

    k = k + 1
  enddo
enddo


write(20, '(a)') 'plot 1/0 w l palette title "" , \'
write(20, '(a)') '     1/0 ls 1 title "ANODE (ITO)"  , \'
write(20, '(a)') '     1/0 ls 2 title "DONNOR (P3HT)" , \'
write(20, '(a)') '     1/0 ls 3 title "ACCEPTOR (PCBM)" , \'
write(20, '(a)') '     1/0 ls 4 title "CATHODE (Ag)" '


end subroutine GnuFiles




end module gnugraphs_m
