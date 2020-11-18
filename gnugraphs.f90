module gnugraphs_m
use parameters_m
use constants_m
use overlap_m
use types_m


contains

subroutine calcula_cor(charge_site, hue_color, color_intensity) 
    implicit none
    !args
    real*8, intent(in) :: charge_site(nr, nc)
    real*8, intent(out) :: hue_color, color_intensity(nr, nc)

    
    hue_color = 0.7d0


    color_intensity = charge_site

end subroutine calcula_cor




subroutine GnuFiles(pop_el, figure_label)
    implicit none
    !args
    real*8, intent(in) :: pop_el(nr, nc)
    integer, intent(in) :: figure_label


    !local
    real*8, parameter :: scale_nanometro = 1.d-9
    real*8, parameter :: scale_freq = 1.d13
    integer :: p
    real*8 :: charge_site


    real*8 :: xmin, xmax, ymin, ymax
    real*8 :: hue_color, color_intensity(nr, nc)
    real*8 :: time_position_x, time_position_y



    xmin = -1.d0
    xmax = 2.d0 * float(nsites) * raioZero / scale_nanometro + 2.d0 * distance_molecule / scale_nanometro
    ymin = -1.d0
    ymax = 2.d0 
    
    time_position_x = site(1, 1)%xPos/scale_nanometro
    time_position_y = ymax - 0.5d0
    
    !write(20, '(a)') '#!/usr/bin/gnuplot'
    
    write(20, '(a)') 'clear'
    write(20, '(a)') 'reset'
    write(20, '(a)') "set terminal pngcairo dashed size 1080, 700 enhanced font 'Verdana, 20'"
    write(20, '(a)') 'set label'
    write(20, '(a)') 'set key'
    write(20, '(a)') "set encoding iso_8859_1"
    write(20, '(a)') "set size ratio -1" !deixa tudo na mesma escala
    write(20, '(a)') "set cbrange[0:1]"
    write(20, '(a)') "set cbtics"
    write(20, '(a)') "set colorbox"
    write(20, '(a)') 'set palette model HSV defined (0 0.7 0 1, 1 0.7 1 1)'
    
    write(20, '(a, F8.3, a, F8.3, a)') 'set xrange [',xmin,':',xmax,']'
    write(20, '(a, F8.3, a, F8.3, a)') 'set yrange [',ymin,':',ymax,']'
    write(20, '(a)') 'set xlabel "x(nm)" '
    write(20, '(a)') 'set ylabel "y(nm)" '

    write(20, '(a)') 'set cblabel "{/Symbol r}_{el}" '
    
    
    write(20, '(a)') 'set style line 1 lw 2 dt 2 lc rgb "black"'
    write(20, '(a)') 'set style line 2 lw 2 dt 3 lc rgb "black"'
    write(20, '(a)') 'set style line 3 lw 2 dt 4 lc rgb "black"'
    write(20, '(a)') 'set style line 4 lw 2 dt 5 lc rgb "black"'
    
    
    
    write(20, '(a, F10.4, a, F8.3, a, F8.3, a)') 'set label 1 sprintf("TIME = %4.4f ps",  ',time,') at ',time_position_x,',',time_position_y,' font "arialbd, 22"'
    
    
    write(20,fmt='(a)')        "system('mkdir -p animacao')"
    write(20,fmt='(a, I5, a)') "fname = sprintf('animacao/animation_%04d.png',",figure_label,")"
    write(20, '(a)') 'set output fname'
    
    
        call calcula_cor( pop_el, hue_color, color_intensity )   
    
    k = 1
    do c = 1, nc
        do r = 1, nr
    
            write(20,fmt='(a, I3, a, F8.3, a, F8.3, a, F7.3, a, F8.3, a, F8.3, a )') &
        'set object ',k,' circle center',site(r, c)%xPos/scale_nanometro,',',site(r, c)%yPos/scale_nanometro,' size &
        ',site(r, c)%radius/scale_nanometro , 'fc rgb hsv2rgb(',hue_color,',', color_intensity(r, c),', 1.0) &
        fs solid 3.0 border rgb "black" lw 2.0 dashtype 1'
    
        k = k + 1
        enddo
    enddo
    
    
    write(20, '(a)') 'plot 1/0 w l palette title "" , \'
    write(20, '(a)') '     1/0 ls 1 title ""  '


end subroutine GnuFiles




end module gnugraphs_m
