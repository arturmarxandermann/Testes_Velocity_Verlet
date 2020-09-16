module overlap_m

use types_m
use constants_m
use parameters_m

contains


subroutine define_sitios(sites_array)
!subrotina que define todos os sitios da celula OPV
implicit none
type(quantum_site),  DIMENSION(:, :), allocatable, intent(out) :: sites_array


real*8, parameter :: scale_nanometro = 1.d-9

!----------------------------------------

allocate(sites_array(nm_rows, nm_columns))



do j = 1, nm_columns
  do i = 1, nm_rows
      sites_array(i, j)%t = siteCoupling
  enddo
enddo

  

do j = 1, nm_columns
  do i = 1, nm_rows
      sites_array(i, j)%homo_energy = 0.d0 
      sites_array(i, j)%lumo_energy = 0.d0 
  enddo
enddo

!SITE RADIUS EM NANOMETROS
do j = 1, nm_columns
  do i = 1, nm_rows
      sites_array(i, j)%mass = siteMass  
      sites_array(i, j)%radiuszero = raioZero
      sites_array(i, j)%radius = sites_array(i, j)%radiuszero
      sites_array(i, j)%radial_vel = 0.d0 
  enddo
enddo

do j = 1, nm_columns
  do i = 1, nm_rows
     sites_array(i, j)%omega = ( 2.d0 * hbar / ( me  * (sites_array(i, j)%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
     sites_array(i, j)%omegazero = sites_array(i, j)%omega * thz_to_hz 
     !omegazero é igual ao omega no instante inicial
  enddo
enddo


do j = 1, nm_columns
  do i = 1, nm_rows
    sites_array(i, j)%site_type_el = (sites_array(i, j)%lumo_energy) 
  enddo
enddo




do j = 1, nm_columns
  do  i = 1, nm_rows
    !PRIMEIRO SITIO LOCALIZADO EM 0,0
    if (i == 1 .AND. j == 1) then
      sites_array(i, j)%posicao_x = 0.d0
      sites_array(i, j)%posicao_y = 0.d0

    else
      sites_array(i, j)%posicao_x = sites_array(i , j-1)%posicao_x + sites_array(i, j-1)%radius + sites_array(i, j)%radius &
        + distance_molecule
      sites_array(i, j)%posicao_y = 0.d0

    endif
  enddo
enddo



do j = 1, nm_columns
  do i = 1, nm_rows
print*, "tamanho molecula", i, j, "em nm é", 2.d0 * sites_array(i, j)%radius * 1.d9
  enddo
enddo

do j = 1, nm_columns
  do i = 1, nm_rows
print*, "omega da molecula", i, j, "em Thz",  sites_array(i, j)%omega
  enddo
enddo

do j = 1, nm_columns
  do i = 1, nm_rows
print*, "V type da molecula el", i, j, "em ev é",  sites_array(i, j)%site_type_el
  enddo
enddo


do j = 1, nm_columns
  do i = 1, nm_rows
print*, "HBAR * omega da molecula:", i, j, "em ev é",  HB_ev_ps*sites_array(i, j)%omega
  enddo
enddo


print*, "distancia entre os sitios utilizada:", distance_molecule


return

end subroutine define_sitios

subroutine DerivativeOverlap(sites_array, nstates, l1, c1, l2, c2, DerOvlp)        
implicit none
!SUBROTINA PARA CALCULAR OS ELEMENTOS <sitio(l1, c1, x-x0, y-y0) | d/dR sitio(l2, c2, x, y) >
type(quantum_site), intent(in) :: sites_array(nm_rows, nm_columns)
integer, intent(in) :: nstates, c1, l1, c2, l2 
real*8, allocatable, intent(out) :: DerOvlp(:, :)

real*8 :: x1, y1, w1, x2, y2, w2, x12, y12, d, d2, d3, d4, sinth, costh, sinth2, costh2, sinth3, costh3
real*8 :: sqmhbar, mhbar, sqhbarm, hbarm
real*8 :: wproduct, wproduct2, wproduct3, wproduct4, wproduct5 
real*8 :: wden1, wden2, wden3, wden4, wden5, wden6, wden7, wden8, wden9, wden10
real*8 :: omega1MinusOmega2, omega1PlusOmega2
real*8 :: omegaMinusPlus, omegaMinusPlus2
real*8 :: w2Sw1half, w2Sw1half2, w1Sw2half, w1Sw2half2
integer :: Q_numbers
real*8 :: EXPd2, omegaDerivative
real*8 :: me2, hbar2, sin2th, cos2th


integer, DIMENSION(6) :: Qnumbers_row !fornece o nm quantico dos estados

allocate(DerOvlp(nstates, nstates), source = R_zero )

Qnumbers_row = (/ 00, 01, 10, 02, 11, 20 /)

x1 = 0.d0 !sites_array(l1, c1)%posicao_x !0.d0 
y1 = 0.d0 !sites_array(l1, c1)%posicao_y !0.d0 

w1 = sites_array(l1, c1)%omega*thz_to_hz

x2 = sites_array(l2, c2)%posicao_x - sites_array(l1, c1)%posicao_x
y2 = sites_array(l2, c2)%posicao_y - sites_array(l1, c1)%posicao_y

w2 = sites_array(l2, c2)%omega*thz_to_hz


x12 = x1 - x2
y12 = y1 - y2



d   = dsqrt(x12*x12 + y12*y12)


d2 = d**2.0
d3 = d**3.0
d4 = d**4.0

if (d > tol) then
    sinth = y12 / d
    costh = x12 / d
else
    sinth = R_ZERO
    costh = R_ONE
endif



omegaDerivative = - ((4.d0 * hbar) / (me * (sites_array(l2, c2)%radius)**3.0)  )

sinth2 = sinth**2.0d0
sinth3 = sinth**3.0d0

costh2 = costh**2.0d0
costh3 = costh**3.0d0

sin2th = 2.d0 * sinth * costh

cos2th = costh2 - sinth2

sqmhbar = sqrt( (me/hbar) )
mhbar = sqmhbar**2.0d0


sqhbarm = sqrt( (hbar/me) )
hbarm = sqhbarm**2.0d0

wproduct = sqrt(w1*w2)
wproduct2 = wproduct**2.0d0
wproduct3 = wproduct**3.0d0
wproduct4 = wproduct**4.0d0
wproduct5 = wproduct**5.0d0

omega1MinusOmega2 = w1 - w2
omega1PlusOmega2 = w1 + w2
omegaMinusPlus = (w1 - w2) * (w1 + w2)
omegaMinusPlus2 = (w1 - w2) * ( (w1 + w2)**2.0 )

w1Sw2half2 = (w1/w2)
w1Sw2half =  sqrt(w1Sw2half2)

w2Sw1half2 = (w2/w1)
w2Sw1half = sqrt(w2Sw1half2)

wden1 = 1.d0 / ( sqrt(omega1PlusOmega2) )
wden2 = wden1**2.d0
wden3 = wden1**3.d0
wden4 = wden1**4.d0
wden5 = wden1**5.d0
wden6 = wden1**6.d0
wden7 = wden1**7.d0
wden8 = wden1**8.d0
wden9 = wden1**9.d0
wden10 = wden1**10.d0

hbar2 = hbar**2.d0
me2 = me**2.d0



EXPd2 = EXP( - HALF * d2 * mhbar * wden2 * wproduct2)



do j = 1 , nstates
    do i = 1 , nstates
  Q_numbers = 100*Qnumbers_row(i) + 1*Qnumbers_row(j)
  select case (Q_numbers)

case(0000)

DerOvlp(i, j) = -((EXPd2*w1Sw2half*(d2*me*w1**2.0*w2 + hbar*(-w1**2.0 + w2**2.0))*wden6)/hbar)

case(0001)

DerOvlp(i, j) = ( SQRT2*d*EXPd2*me*sinth*w1**1.5*(d2*Sqrt(hbar*me)*w1**2.0*w2+2.d0*hbar*sqhbarm*(-w1**2.0 + w2**2.0))*wden8)/hbar**2.0

case(0010)
DerOvlp(i, j) = SQRT2*multFactorDer*costh*d*EXPd2*Sqrt(me)*(w1/hbar)**1.5d0*(-(d2*me*w1**2.d0*w2) + 2.d0*hbar*(w1**2.0 - w2**2.0))*wden8


case(0100) 

DerOvlp(i, j) = -((SQRT2*d*EXPd2*me*sinth*w1*Sqrt(w2)*(-(hbar*omega1PlusOmega2*sqhbarm*(3.d0*w1 - w2)) + & 
                d2*Sqrt(hbar*me)*w1**2.0*w2)*wden8)/hbar**2.0)        



case(0101)

DerOvlp(i, j) = (2.d0*EXPd2*w1*wden10*(2.d0*hbar2*omegaMinusPlus2 + d4*me2*sinth2*w1**3.0*w2**2.0 + &
      d2*hbar*me*omega1PlusOmega2*(-3.d0*w1 + cos2th*(2.d0*w1 - w2) + w2)*wproduct2))/hbar**2.0


case(0110)


DerOvlp(i, j) = (d2*EXPd2*me*sin2th*w1**2.0*w2*(-2.d0*hbar*omega1PlusOmega2*(2.d0*w1 - w2) + d2*me*w1**2.0*w2)*wden10)/hbar**2.0


case(1000)
DerOvlp(i, j) = -((multFactorDer*SQRT2*EXPd2*costh*d*w1*Sqrt(me*w2)*(-(hbar*omega1PlusOmega2*(3.d0*w1 - w2)) + &
                d2*me*w1**2.d0*w2)*wden8)/(hbar**1.5d0)) 

case(1001) 

DerOvlp(i, j) = (d2*EXPd2*me*sin2th*w1**2.0*w2*(-2.d0*hbar*omega1PlusOmega2*(2.d0*w1 - w2) + d2*me*w1**2.0*w2)*wden10)/hbar**2.0

case(1010)

DerOvlp(i, j) = (2.d0*w1* EXPd2*wden10*(costh2*d4*me2*w1**3.0*w2**2.0 + hbar*omega1PlusOmega2*(2.d0*hbar*(w1**2.0 - w2**2.0) + & 
        d2*me*(-3.d0*w1 + w2)*wproduct2 + cos2th*d2*me*(-2.d0*w1 + w2)*wproduct2)))/(hbar**2) 

    end select


    enddo
enddo

DerOvlp = DerOvlp * omegaDerivative 

return
end subroutine DerivativeOverlap









subroutine NewOverlap(sites_array, nstates, l1, c1, l2, c2, Qnumbers_row, S)
implicit none
type(quantum_site), intent(in) :: sites_array(nm_rows, nm_columns)
integer, intent(in) :: nstates, c1, l1, c2, l2
integer, DIMENSION(6), INTENT(IN) :: Qnumbers_row !Qnumbers_row faz o papel do antigo basis -> fornece o nm quantico dos estados
real*8, allocatable, intent(out) :: S(:, :)

real*8 :: x1, y1, w1, x2, y2, w2, x12, y12, d, d2, d3, sinth, costh, sinth2, costh2, sinth3, costh3
real*8 :: sqmhbar, mhbar, mhbar2, sqhbarm, hbarm
real*8 :: wproduct, wproduct2, wproduct3, wproduct4, wproduct5 
real*8 :: wden1, wden2, wden3, wden4, wden5, wden6, wden7, wden8, wden9, wden10
integer :: Q_numbers
real*8 :: EXPd2


allocate(S(nstates, nstates), source = R_zero )


x1 = 0.d0 !sites_array(l1, c1)%posicao_x
y1 = 0.d0 !sites_array(l1, c1)%posicao_y
w1 = sites_array(l1, c1)%omega*thz_to_hz 
x2 = sites_array(l2, c2)%posicao_x - sites_array(l1, c1)%posicao_x
y2 = sites_array(l2, c2)%posicao_y - sites_array(l1, c1)%posicao_y
w2 = sites_array(l2, c2)%omega*thz_to_hz 
x12 = x1 - x2
y12 = y1 - y2



d   = dsqrt(x12*x12 + y12*y12)


d2 = d**2.0
d3 = d**3.0



if (d > tol) then
    sinth = y12 / d
    costh = x12 / d
else
    sinth = R_ZERO
    costh = R_ONE
endif


sinth2 = sinth**2.0
sinth3 = sinth**3.0

costh2 = costh**2.0
costh3 = costh**3.0

sqmhbar = sqrt( (me/hbar) )
mhbar = sqmhbar**2.0
mhbar2 = sqmhbar**4.0


sqhbarm = sqrt( (hbar/me) )
hbarm = sqhbarm**2.0

wproduct = sqrt(w1*w2)
wproduct2 = wproduct**2.0
wproduct3 = wproduct**3.0
wproduct4 = wproduct**4.0
wproduct5 = wproduct**5.0

wden1 = 1.d0 / ( sqrt(w1 + w2) )
wden2 = wden1**2.0
wden3 = wden1**3.0
wden4 = wden1**4.0
wden5 = wden1**5.0
wden6 = wden1**6.0
wden7 = wden1**7.0
wden8 = wden1**8.0
wden9 = wden1**9.0
wden10 = wden1**10.0


EXPd2 = EXP( - HALF * d2 * mhbar * wden2 * wproduct2)



do j = 1 , nstates
    do i = 1 , nstates
  Q_numbers = 100*Qnumbers_row(i) + 1*Qnumbers_row(j)
  select case (Q_numbers)

  
case(0000)
S(i, j) = 2.d0*EXPd2*wden2*wproduct

case(0001)
S(i, j) = (2.d0*SQRT2*d*EXPd2*sinth*Sqrt(mhbar*w2)*(-(hbar*Sqrt(w1 + w2)) + me*sqhbarm*w2*Sqrt(hbarm*wden2))*wden3*wproduct)/hbar

case(0010)

S(i, j) = -2.d0*SQRT2*costh*d*EXPd2*sqmhbar*(w1**1.5d0)*w2*wden4*multFactorOvlp


case(0100) 
S(i, j) = ((2.d0*SQRT2)*d*EXPd2*sinth*Sqrt(w1)*Sqrt(mhbar*w1)*(w2**1.5)*Sqrt(hbarm*wden2)*wden3)/Sqrt(hbarm)

case(0101)
S(i, j) = (4.d0*EXPd2*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*( (d2*me*Sqrt(hbar*me)*sinth2*sqhbarm*(w2**2.0)) + ((hbar**2.0)*(w1 + w2)) - &
          (d2*hbar*me*sinth2*w2*(w1 + w2)))*Sqrt(hbarm*wden2)*wden5*wproduct)/((hbarm**1.5)*(me**2.0))


case(0110)
S(i, j) = -4.d0*costh*d2*EXPd2*sinth*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*wden6*wproduct3

case(1000)

S(i, j) = 2.d0*SQRT2*costh*d*EXPd2*sqmhbar*w1*(w2**1.5d0)*wden4*multFactorOvlp 


case(1001) 
S(i, j) = (-4.d0*costh*d2*EXPd2*sinth*Sqrt(w1)*Sqrt(mhbar*w1)*(w2**1.5)*Sqrt(mhbar*w2)*wden7* &
        ((w1**2.0)*Sqrt(hbarm*wden2) + w2*(-(sqhbarm*Sqrt(w1 + w2)) + w2*Sqrt(hbarm*wden2)) + 2.d0*Sqrt(hbarm*wden2)*wproduct2))/Sqrt(hbarm)

case(1010)
S(i, j) = (4.d0*EXPd2*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*Sqrt(hbarm*wden2)*wden5*wproduct*(hbar*sqhbarm*(w1 + w2) - &
           costh2*d2*Sqrt(hbar*me)*wproduct2))/hbar 


    end select


    enddo
enddo

return
end subroutine NewOverlap





subroutine Overlap(sites_array, nstates, l1, c1, l2, c2, Qnumbers_row, S)
!calcula o overlap S entre o sitio1(x, y) e o sitio2(x, y)
!para todos os estados harmonicos
implicit none
type(quantum_site), intent(in) :: sites_array(nm_rows, nm_columns)
integer, intent(in) :: nstates, c1, l1, c2, l2
integer, DIMENSION(6), INTENT(IN) :: Qnumbers_row !Qnumbers_row faz o papel do antigo basis -> fornece o nm quantico dos estados
real*8, allocatable, intent(out) :: S(:, :)

!variaveis locais

integer, allocatable  :: basis(:, :)
integer ::  i , j , Q_numbers
real*8 :: x1 , y1 , x2 , y2 , x12 , y12 , d, d2, d3, d4 , d5 , d6 , d7 , d8 , renso , renso2 , renso3 , renso4 , renso5 , renso6
real*8 :: w1 , w2 , mdhbar , sq_mdhbar, sq_w1 , sq_w2 ,  w1_2 , w2_2, w1_3 , w2_3, w1_4 , w2_4 , w1_5 , w2_5 , w1_6 , w2_6 , w1_7 , w2_7, w1_8 , w2_8
real*8 :: nu0 , nu2 , nu4, nubar , nubar2 , nubar3 , nubar4 , nubar5, mdhbar2, mdhbar3 , mdhbar4
real*8 :: sin_th , sin_th2 , sin_th3 , sin_th4 , sin_th5 , sin_th6 , sin_th7 , sin_th8  , sin_2th , sin_2th2, sin_2th3 , sin_2th4
real*8 :: cos_th , cos_th2 , cos_th3 , cos_th4 , cos_th5 , cos_th6, cos_th7 , cos_th8 , EXPd2 , renso7 , renso8


allocate(S(nstates, nstates), source = R_zero )
allocate(basis(nstates, 2), source = 0)


x1 = 0.d0 !sites_array(l1, c1)%posicao_x - sites_array(l2, c2)%posicao_x
y1 = 0.d0 !sites_array(l1, c1)%posicao_y - sites_array(l2, c2)%posicao_y
w1 = sites_array(l1, c1)%omega*thz_to_hz
x2 = sites_array(l2, c2)%posicao_x - sites_array(l1, c1)%posicao_x
y2 = sites_array(l2, c2)%posicao_y - sites_array(l1, c1)%posicao_y
w2 = sites_array(l2, c2)%omega*thz_to_hz
x12 = x1- x2
y12 = y1- y2

nu0 = (w1*w2)/(w1+w2)
nubar = sqrt((w1*w2))/(w1+w2)


mdhbar = me/hbar

mdhbar2 = mdhbar*mdhbar
mdhbar3 = mdhbar*mdhbar*mdhbar
mdhbar4 = mdhbar*mdhbar*mdhbar*mdhbar


sq_mdhbar = dsqrt(me/hbar)
!print*,"sq_mdhbar", sq_mdhbar

d   = dsqrt(x12*x12 + y12*y12)



if (d > tol) then
    sin_th = y12 / d
    cos_th = x12 / d
else
    sin_th = R_ZERO
    cos_th = R_ONE
endif


w1_2 = w1*w1
w2_2 = w2*w2

w1_3 = w1*w1*w1
w2_3 = w2*w2*w2

w1_4 = w1*w1*w1*w1
w2_4 = w2*w2*w2*w2

w1_5 = w1*w1*w1*w1*w1
w2_5 = w2*w2*w2*w2*w2

w1_6 = w1*w1*w1*w1*w1*w1
w2_6 = w2*w2*w2*w2*w2*w2

w1_7 = w1*w1*w1*w1*w1*w1*w1
w2_7 = w2*w2*w2*w2*w2*w2*w2

w1_8 = w1*w1*w1*w1*w1*w1*w1*w1
w2_8 = w2*w2*w2*w2*w2*w2*w2*w2


sq_w1 = dsqrt(w1)
sq_w2 = dsqrt(w2)

nu2 = nu0 * nu0
nu4 = nu0 * nu0 * nu0 * nu0

nubar2 = nubar*nubar
nubar3 = nubar*nubar*nubar
nubar4 = nubar*nubar*nubar*nubar
nubar5 = nubar*nubar*nubar*nubar*nubar

renso = R_ONE/(w1+w2)
renso2 = renso*renso
renso3 = renso*renso*renso
renso4 = renso*renso*renso*renso
renso5 = renso*renso*renso*renso*renso
renso6 = renso*renso*renso*renso*renso*renso
renso7 = renso*renso*renso*renso*renso*renso*renso
renso8 = renso*renso*renso*renso*renso*renso*renso*renso

d2 = d*d
d3 = d2*d
d4 = d2*d2
d5 = d4*d
d6 = d5*d
d7 = d6*d
d8 = d7*d

sin_th2 = sin_th*sin_th
sin_th3 = sin_th2*sin_th
sin_th4 = sin_th2*sin_th2
sin_th5 = sin_th4*sin_th
sin_th6 = sin_th5*sin_th
sin_th7 = sin_th6*sin_th
sin_th8 = sin_th7*sin_th

sin_2th = R_TWO*cos_th*sin_th
sin_2th2 = sin_2th*sin_2th
sin_2th3 = sin_2th2*sin_2th
sin_2th4 = sin_2th3*sin_2th

cos_th2 = cos_th*cos_th
cos_th3 = cos_th2*cos_th
cos_th4 = cos_th2*cos_th2
cos_th5 = cos_th4*cos_th
cos_th6 = cos_th5*cos_th
cos_th7 = cos_th6*cos_th
cos_th8 = cos_th7*cos_th

EXPd2 = EXP(-(mdhbar*nu0*d*d)/2.d0)


!Q_NUMBERS DESSA FORMA FORNECE UMA MATRIZ DE OVERLAP NA FORMA
!S(i, j) = (0000 0001 0010)
!          (0100 0101 0110)
!          (1000 1001 1010)

do j = 1 , nstates
    do i = 1 , nstates

! if (j>=i) then

Q_numbers = 100*Qnumbers_row(i) + 1*Qnumbers_row(j)

select case (Q_numbers)

case( 0000 )
S(i,j) = R_TWO*nubar*EXPd2

case( 0001 )
S(i,j) = -R_TWO*SQRT2*sq_mdhbar*sq_w1*nubar2*EXPd2*d*sin_th

case( 0002 )
S(i,j) = SQRT2*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 0010 )
S(i,j) = -R_TWO*SQRT2*sq_mdhbar*sq_w1*nubar2*EXPd2*d*cos_th

case( 0011)
S(i,j) = R_TWO*mdhbar*w1*d2*sin_2th*nubar3*EXPd2

case( 0020 )
S(i,j) = SQRT2*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 0100 )
S(i,j) = R_TWO*SQRT2*sq_mdhbar*sq_w2*nubar2*EXPd2*d*sin_th

case( 0101 )
S(i,j) = R_FOUR*nubar2*(R_ONE-mdhbar*nu0*d2*sin_th2)*EXPd2

case( 0102)
S(i,j) = R_TWO*d*sin_th*nubar2*sq_mdhbar*sq_w2*(w2_2-R_FOUR*w1*w2+w1_2*(R_TWO*mdhbar*w2*d2*sin_th2-R_FIVE))*EXPd2*renso2

case( 0110 )
S(i,j) = -R_TWO*nubar2*mdhbar*nu0*d2*sin_2th*EXPd2

case( 0111 )
S(i,j) = -R_FOUR*SQRT2*d*cos_th*nubar2*w1*sq_mdhbar*sq_w2*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 0120 )
S(i,j) = R_TWO*d*sin_th*nubar2*sq_mdhbar*sq_w2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 0200 )
S(i,j) = SQRT2*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso2

case( 0201 )
S(i,j) = R_TWO*sq_mdhbar*sq_w1*d*sin_th*nubar2*(R_FOUR*w1*w2-w1_2+w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso2

case( 0202 )
S(i,j) =nubar*(R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_FOUR-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 0210 )
S(i,j) = R_TWO*sq_mdhbar*sq_w1*d*cos_th*nubar2*(w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)-w1_2)*EXPd2*renso2

case( 0211 )
S(i,j) = SQRT2*d2*sin_2th*mdhbar*w1*nubar3*EXPd2*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*renso2

case( 0220 )
S(i,j) = EXPd2*renso4*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))

case( 1000 )
S(i,j) = R_TWO*SQRT2*sq_mdhbar*sq_w2*nubar2*EXPd2*d*cos_th

case( 1001 )
S(i,j) = -R_TWO*nubar2*mdhbar*nu0*d2*sin_2th*EXPd2

case( 1002 )
S(i,j) = R_TWO*sq_mdhbar*sq_w2*d*cos_th*nubar2*(w2_2 - w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 1010 )
S(i,j) = R_FOUR*nubar2*(R_ONE-mdhbar*nu0*d2*cos_th2)*EXPd2

case( 1011 )
S(i,j) = -R_FOUR*SQRT2*sq_mdhbar*sq_w2*w1*nubar2*d*sin_th*(w2+w1*(R_ONE - mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 1020 )
S(i,j ) = R_TWO*sq_mdhbar*sq_w2*nubar2*d*cos_th*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 1100 )
S(i,j) = R_TWO*mdhbar*w2*d2*sin_2th*nubar3*EXPd2

case( 1101 )
S(i,j) = R_FOUR*SQRT2*d*cos_th*nubar2*w2*sq_mdhbar*sq_w1*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 1102 )
S(i,j) = SQRT2*mdhbar*w2*d2*sin_2th*nubar3*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 1110 )
S(i,j) = R_FOUR*SQRT2*sq_mdhbar*sq_w1*w2*d*sin_th*nubar2*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 1111 )
S(i,j) = R_EIGHT*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 1120 )
S(i,j) = SQRT2*mdhbar*w2*nubar3*d2*sin_2th*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 2000 )
S(i,j) = SQRT2*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso2

case( 2001 )
S(i,j) = -R_TWO*sq_mdhbar*sq_w1*d*sin_th*nubar2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso2

case( 2002 )
S(i,j) = nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 2010 )
S(i,j) = -R_TWO*d*cos_th*nubar2*sq_mdhbar*sq_w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 2011 )
S(i,j) = SQRT2*mdhbar*w1*nubar3*d2*sin_2th*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 2020 )
S(i,j) = nubar*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 0003 )
S(i,j) =R_TWO*SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w1*(-R_TREE*w2_2+w1_2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 0012 )
S(i,j) = -R_TWO*sq_mdhbar*sq_w1*nubar2*d*cos_th*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case( 0021 )
S(i,j) = -R_TWO*sq_mdhbar*sq_w1*nubar2*d*sin_th*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 0030 )
S(i,j) = -R_TWO*SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w1*(+R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case( 0103 )
S(i,j) = R_TWO*SQRT2*SQRT1d3*nubar2*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso3

case( 0112 )
S(i,j) = - SQRT2*mdhbar*nubar4*d2*sin_2th*(-R_FOUR*w1*w2 +w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso

case( 0121 )
S(i,j) = R_TWO*SQRT2*nubar2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 0130 )
S(i,j) = -SQRT2*SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso

case( 0203 )
S(i,j) = -SQRT2*SQRT1d3*d*sin_th*sq_mdhbar*sq_w1*nubar2*(-R_TREE*w1_4+R_TWO*w1_3*w2*(R_EIGHTEEN+mdhbar*w1*d2*sin_th2)+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w1*d2*sin_th2)+R_TREE*w2_4*(-R_FIVE+R_TWO*mdhbar*w1*d2*sin_th2)+R_TWO*w1*w2_3*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)*(R_SIX-mdhbar*w1*d2*sin_th2))*EXPd2*renso4

case( 0212 )
S(i,j) = SQRT2*sq_mdhbar*sq_w1*nubar2*d*cos_th*(w2_4+w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_FOUR-mdhbar*w2*d2*sin_th2)-R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 0221 )
S(i,j) = SQRT2*sq_mdhbar*sq_w1*nubar2*d*sin_th*(-w1_2+R_FOUR*w1*w2+w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 0230 )
S(i,j) = -SQRT2*SQRT1d3*d*cos_th*sq_mdhbar*sq_w1*nubar2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 0300 )
S(i,j) = R_TWO*SQRT1d3*d*sin_th*sq_mdhbar*sq_w2*nubar2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso2

case( 0301 )
S(i,j) =R_TWO*SQRT2*SQRT1d3*nubar2*(R_TREE*w1_3+R_TREE*w1_2*w2*(R_ONE-mdhbar*w1*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_3*(-R_TREE+mdhbar*w1*d2*sin_th2*(R_NINE-R_TWO*mdhbar*w1*d2*sin_th2)))*EXPd2*renso3

case( 0302 )
S(i,j) = SQRT2*SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w2*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_SIX-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 0303 )
S(i,j) = -R_TWO*ONE_THIRD*nubar2*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*sin_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN - R_THIRTY_THREE*mdhbar*w2*d2*sin_th2+R_SEVEN*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*sin_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar2*mdhbar*w2_3*d6*sin_th6))*EXPd2*renso5

case( 0310 )
S(i,j) = -SQRT2*SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso

case( 0311 )
S(i,j) = -R_FOUR*SQRT1d3*d*cos_th*nubar2*w1*sq_mdhbar*sq_w2*(R_TREE*w1_3+R_TREE*w1_2*w2*(R_ONE-mdhbar*w1*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_3*(-R_TREE+mdhbar*w1*d2*sin_th2*(R_NINE-R_TWO*mdhbar*w1*d2*sin_th2)))*EXPd2*renso4

case( 0312 )
S(i,j) = SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w2_4+R_TREE*w1_4*(R_FIVE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_SIX-mdhbar*w2*d2*sin_th2)-R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_SIX*w1_2*w2_2*(-R_ELEVEN+R_FIVE*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 0320 )
S(i,j) = SQRT2*SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 0321 )
S(i,j) = R_TWO*SQRT1d3*nubar2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case( 0330 )
S(i,j) = -ONE_THIRD*d2*sin_2th*mdhbar*nubar4*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case( 1003 )
S(i,j) = SQRT2*SQRT1d3*d2*sin_2th*mdhbar*nubar4*(-R_TREE*w2_2+w1_2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso

case( 1012 )
S(i,j) = R_TWO*SQRT2*nubar2*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 1021 )
S(i,j) = -SQRT2*mdhbar*nubar4*d2*sin_2th*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso

case( 1030 )
S(i,j) = R_TWO*SQRT2*SQRT1d3*nubar2*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso3

case( 1103 )
S(i,j) = R_FOUR*SQRT1d3*sq_mdhbar*sq_w1*w2*nubar2*d*cos_th*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(-R_TREE+mdhbar*w2*d2*sin_th2*(R_NINE-R_TWO*mdhbar*w2*d2*sin_th2)))*EXPd2*renso4

case( 1112 )
S(i,j) = R_FOUR*sq_mdhbar*sq_w2*nubar3*d*sin_th*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 1121 )
S(i,j) = R_FOUR*sq_mdhbar*sq_w2*nubar3*d*cos_th*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 1130 )
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar2*w2*sq_mdhbar*sq_w1*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case( 1200 )
S(i,j) = R_TWO*sq_mdhbar*sq_w2*nubar2*d*cos_th*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso2

case( 1201 )
S(i,j) = SQRT2*mdhbar*nubar4*d2*sin_2th*(-w1_2+R_FOUR*w1*w2+w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso

case( 1202 )
S(i,j) = SQRT2*sq_mdhbar*sq_w2*nubar2*d*cos_th*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_FOUR-mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 1203 )
S(i,j) = SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w1_4-R_TWO*w1_3*w2*(R_EIGHTEEN+mdhbar*w1*d2*sin_th2)+R_SIX*w1_2*w2_2*(-R_ELEVEN+R_FIVE*mdhbar*w1*d2*sin_th2)-R_TWO*w1*w2_3*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)*(R_SIX-mdhbar*w1*d2*sin_th2)+R_TREE*w2_4*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso3

case( 1210 )
S(i,j) = R_TWO*SQRT2*nubar2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case( 1211 )
S(i,j) = R_FOUR*sq_mdhbar*sq_w1*nubar3*d*sin_th*(-w1_2+R_FOUR*w1*w2+w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case( 1212 )
S(i,j) = R_TWO*nubar2*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_FOUR-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2))*EXPd2*renso5

case( 1220 )
S(i,j) = SQRT2*sq_mdhbar*sq_w2*nubar2*d*cos_th*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 1221 )
S(i,j) =mdhbar*nubar4*d2*sin_2th*(-w1_2+R_FOUR*w1*w2+w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case( 1230 )
S(i,j) =R_TWO*SQRT1d3*nubar2*(-w1_2+w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5


case( 2003 )
S(i,j) = -SQRT2*SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w1*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(R_TREE*w2_2+w1_2*(-R_TREE + R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 2012 )
S(i,j) = -SQRT2*d*cos_th*nubar2*sq_mdhbar*sq_w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 2021 )
S(i,j) = SQRT2*sq_mdhbar*sq_w1*d*sin_th*nubar2*(w2_4+w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)-R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 2030 )
S(i,j) = -SQRT2*SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w1*(-R_TREE*w1_4+R_TWO*w1_3*w2*(R_EIGHTEEN+mdhbar*w1*d2*cos_th2)+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w1*d2*cos_th2)+R_TWO*w1*w2_3*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)*(R_SIX-mdhbar*w1*d2*cos_th2)+R_TREE*w2_4*(-R_FIVE+R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso4

case( 2100 )
S(i,j) = R_TWO*sq_mdhbar*sq_w2*nubar2*d*sin_th*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso2

case( 2101 )
S(i,j) = R_TWO*SQRT2*nubar2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 2102 )
S(i,j) = SQRT2*sq_mdhbar*sq_w2*nubar2*d*sin_th*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 2103 )
S(i,j) = -R_TWO*SQRT1d3*nubar2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(R_TREE+mdhbar*w2*d2*sin_th2*(-R_NINE+R_TWO*mdhbar*w2*d2*sin_th2)))*EXPd2*renso5

case( 2110 )
S(i,j) = -SQRT2*mdhbar*nubar4*d2*sin_2th*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*EXPd2*renso

case( 2111 )
S(i,j) = -R_FOUR*sq_mdhbar*sq_w1*nubar3*d*cos_th*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 2112 )
S(i,j) = -mdhbar*nubar4*d2*sin_2th*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(-R_FOUR*w1*w2 +w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 2120 )
S(i,j) = SQRT2*sq_mdhbar*sq_w2*nubar2*d*sin_th*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 2121 )
S(i,j) = -R_TWO*nubar2*(w2_4+w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)-R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso5

case( 2130 )
S(i,j) = SQRT1d3*mdhbar*nubar4*d2*sin_2th*(R_TREE*w1_4-R_TWO*w1_3*w2*(R_EIGHTEEN+mdhbar*w1*d2*cos_th2)+R_SIX*w1_2*w2_2*(-R_ELEVEN+R_FIVE*mdhbar*w1*d2*cos_th2)-R_TWO*w1*w2_3*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)*(R_SIX-mdhbar*w1*d2*cos_th2)+R_TREE*w2_4*(R_FIVE-R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso3

case( 3000 )
S(i,j) = R_TWO*SQRT1d3*d*cos_th*sq_mdhbar*sq_w2*nubar2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso2

case( 3001 )
S(i,j) = -SQRT2*SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso

case( 3002 )
S(i,j) = SQRT2*SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 3003 )
S(i,j) = -ONE_THIRD*d2*sin_2th*mdhbar*nubar4*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case( 3010 )
S(i,j) = R_TWO*SQRT2*SQRT1d3*nubar2*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso3

case( 3011 )
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar2*w1*sq_mdhbar*sq_w2*(R_TREE*w2_3+R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)-R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(-R_TREE-R_SIX*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case( 3012 )
S(i,j) = R_TWO*SQRT1d3*nubar2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case( 3020 )
S(i,j) = SQRT2*SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_SIX-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case( 3021 )
S(i,j) = SQRT1d3*d2*sin_2th*mdhbar*nubar4*(R_TREE*w2_4+R_TREE*w1_4*(R_FIVE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_SIX-mdhbar*w2*d2*cos_th2)-R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_SIX*w1_2*w2_2*(-R_ELEVEN+R_FIVE*mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case( 3030 )
S(i,j) = -R_TWO*ONE_THIRD*nubar2*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*cos_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN - R_THIRTY_THREE*mdhbar*w2*d2*cos_th2+R_SEVEN*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*cos_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar2*mdhbar*w2_3*d6*cos_th6))*EXPd2*renso5

case( 0004 )
S(i,j) = SQRT1d6*nubar*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*sin_th2*(-R_TREE+mdhbar*w2*d2*sin_th2)))*EXPd2*renso4

case( 0013 )
S(i,j) = SQRT2*SQRT1d3*nubar3*mdhbar*w1*d2*sin_2th*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso2

case(0022)
S(i,j) = nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case( 0031 )
S(i,j) = SQRT2*SQRT1d3*nubar3*mdhbar*w1*d2*sin_2th*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso2

case(0040)
S(i,j) = SQRT1d6*nubar*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*cos_th2*(-R_TREE+mdhbar*w2*d2*cos_th2)))*EXPd2*renso4


case(0104)
S(i,j) = SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w2*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(0113)
S(i,j) = R_FOUR*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w1*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso3

case(0122)
S(i,j) = SQRT2*d*sin_th*nubar2*sq_mdhbar*sq_w2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(0131)
S(i,j) = -R_FOUR*SQRT1d3*nubar3*d*cos_th*sq_mdhbar*sq_w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case(0140)
S(i,j) = SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(0204)
S(i,j) = HALF*SQRT1d3*nubar*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*sin_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*sin_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*sin_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0213)
S(i,j) = SQRT1d3*d2*sin_2th*mdhbar*w1*nubar3*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*sin_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(0222)
S(i,j) = SQRT1d2*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(0231)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(0240)
S(i,j) = HALF*SQRT1d3*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(0304)
S(i,j) = ONE_THIRD*SQRT1d2*d*sin_th*nubar2*sq_mdhbar*sq_w2*(-R_NINE*w2_6+R_FIFTEEN*w1_2*w2_4*(R_THIRTEEN-R_TWELVE*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_5*(R_THIRTY_SIX+mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_3*w2_3*(-R_SIXTY+R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_6*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)-R_TREE*w1_4*w2_2*(R_345-R_280*mdhbar*w2*d2*sin_th2+R_THIRTY_SIX*mdhbar2*w2_2*d4*sin_th4) +R_TWO*w1_5*w2*(-R_132+R_195*mdhbar*w2*d2*sin_th2-R_SIXTY*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0313)
S(i,j) = R_TWO*ONE_THIRD*SQRT2*d*cos_th*nubar3*sq_mdhbar*sq_w1*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*sin_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*sin_th2+R_SEVEN*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*sin_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso5

case(0322)
S(i,j) = SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(0331)
S(i,j) = R_TWO*SQRT2*ONE_THIRD*d*cos_th*nubar3*sq_mdhbar*sq_w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*(R_TREE*w2_3+R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)-R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(-R_TREE-R_SIX*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case(0340)
S(i,j) = ONE_THIRD*SQRT1d2*d*sin_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6


case(0400)
S(i,j) = SQRT1d6*nubar*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*sin_th2*(-R_TREE+mdhbar*w1*d2*sin_th2)))*EXPd2*renso4

case(0401)
S(i,j) = SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w1*(-R_TREE*w1_4+R_TWENTY_FOUR*w1_3*w2+R_SIX*w2_2*w1_2*(R_FIVE-R_TWO*mdhbar*w1*d2*sin_th2)+R_EIGHT*w1*w2_3*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(-R_TWENTY_SEVEN-R_FOUR*mdhbar*w1*d2*sin_th2*(-R_SEVEN+mdhbar*w1*d2*sin_th2)))*EXPd2*renso4

case(0402)
S(i,j) = HALF*SQRT1d3*nubar*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*sin_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*sin_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0403)
S(i,j) = -ONE_THIRD*SQRT1d2*d*sin_th*nubar2*sq_mdhbar*sq_w1*(R_EIGHTY_ONE*w2_6-R_TWELVE*w1*w2_5*(R_TWENTY_TWO+R_SEVEN*mdhbar*w2*d2*sin_th2)+w1_6*(-R_NINE+R_SIX*mdhbar*w2*d2*sin_th2)-R_120*w1_3*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TWELVE*w1_5*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_2*w2_4*(-R_345+R_130*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2_2*(R_195+R_180*mdhbar*w2*d2*sin_th2-R_108*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0404)
S(i,j) = ONE_TWELFTH*nubar*( R_NINE*w2_8 -R_THIRTY_SIX*w1*w2_7*(R_EIGHT+mdhbar*w2*d2*sin_th2) +R_TWELVE*w1_3*w2_5*(R_152-R_TWENTY_THREE*mdhbar*w2*d2*sin_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*sin_th4) -R_TWELVE*w1_2*w2_6*(R_NINETEEN-R_NINETY_NINE*mdhbar*w2*d2*sin_th2-mdhbar2*w2_2*d4*sin_th4) +R_TREE*w1_8*(R_TREE-R_TWELVE*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4) +R_FOUR*w1_5*w2_3*(R_456-R_1755*mdhbar*w2*d2*sin_th2+R_840*mdhbar2*w2_2*d4*sin_th4-R_SEVENTY_SIX*mdhbar3*w2_3*d6*sin_th6)-R_TWELVE*w1_7*w2*(R_TWENTY_FOUR-R_NINETY_NINE*mdhbar*w2*d2*sin_th2+R_FOURTY_FOUR*mdhbar2*w2_2*d4*sin_th4-R_FOUR*mdhbar3*w2_3*d6*sin_th6) +R_SIX*w1_4*w2_4*(R_585-R_1170*mdhbar*w2*d2*sin_th2+R_190*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6)+R_FOUR*w1_6*w2_2*(-R_FIFTY_SEVEN-R_SIXTY_NINE*mdhbar*w2*d2*sin_th2+R_285*mdhbar2*w2_2*d4*sin_th4-R_SEVENTY_SIX*mdhbar3*w2_3*d6*sin_th6+R_FOUR*mdhbar4*w2_4*d8*sin_th8) )*EXPd2*renso8

case(0410)
S(i,j) = -SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*sin_th2*(-R_TREE+mdhbar*w1*d2*sin_th2)))*EXPd2*renso4

case(0411)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*sin_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(0412)
S(i,j) = -SQRT1d6*d*cos_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*sin_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*sin_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0413)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w1*(R_EIGHTY_ONE*w2_6-R_TWELVE*w1*w2_5*(R_TWENTY_TWO+R_SEVEN*mdhbar*w2*d2*sin_th2)+w1_6*(-R_NINE+R_SIX*mdhbar*w2*d2*sin_th2)-R_120*w1_3*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TWELVE*w1_5*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_2*w2_4*(-R_345+R_130*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2_2*(R_195+R_180*mdhbar*w2*d2*sin_th2-R_108*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(0420)
S(i,j) = HALF*SQRT1d3*nubar*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*sin_th2*(-R_TREE+mdhbar*w1*d2*sin_th2)))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso6

case(0421)
S(i,j) = -SQRT1d6*d*sin_th*nubar2*sq_mdhbar*sq_w1*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*sin_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(0422)
S(i,j) = HALF*SQRT1d6*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2))*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*sin_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*sin_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4+R_EIGHT*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso8

case(0430)
S(i,j) = -ONE_THIRD*SQRT1d2*d*cos_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*sin_th2*(-R_TREE+mdhbar*w1*d2*sin_th2)))*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso6

case(0431)
S(i,j) = ONE_SIXTH*d2*sin_2th*EXPd2*nubar3*mdhbar*w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2))*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*sin_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*renso6

case(0440)
S(i,j) = ONE_TWELFTH*nubar*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*sin_th2*(-R_TREE+mdhbar*w1*d2*sin_th2)))*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso8

case(1004)
S(i,j) = SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*sin_th2*(-R_TREE+mdhbar*w2*d2*sin_th2)))*EXPd2*renso4

case(1013)
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w1*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_2+w1_2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case(1022)
S(i,j) = SQRT2*d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(1031)
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w1*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso3

case(1040)
S(i,j) = SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*cos_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(1104)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(1113)
S(i,j) = R_FOUR*SQRT2*SQRT1d3*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(1122)
S(i,j) = d2*sin_2th*nubar3*mdhbar*w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(1131)
S(i,j) = R_FOUR*SQRT2*SQRT1d3*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(1140)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*cos_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(1204)
S(i,j) = SQRT1d6*d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*sin_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*sin_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*sin_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(1213)
S(i,j) = -R_TWO*SQRT2*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w1*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*sin_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case(1222)
S(i,j) = d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)*(R_FOUR-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2))*EXPd2*renso6

case(1231)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(1240)
S(i,j) = SQRT1d6*d*cos_th*nubar2*sq_mdhbar*sq_w2*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*cos_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(1300)
S(i,j) = SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso2

case(1301)
S(i,j) = SQRT1d3*R_FOUR*d*cos_th*nubar3*sq_mdhbar*sq_w2*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso3

case(1302)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*d2*w2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(1303)
S(i,j) = -R_TWO*ONE_THIRD*SQRT2*d*cos_th*nubar3*sq_mdhbar*sq_w2*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*sin_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*sin_th2+R_SEVEN*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*sin_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso5

case(1304)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w2*(-R_NINE*w2_6+R_FIFTEEN*w1_2*w2_4*(R_THIRTEEN-R_TWELVE*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_5*(R_THIRTY_SIX+mdhbar*w2*d2*sin_th2)+R_TWELVE*w1_3*w2_3*(-R_SIXTY+R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_6*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)-R_TREE*w1_4*w2_2*(R_345-R_280*mdhbar*w2*d2*sin_th2+R_THIRTY_SIX*mdhbar2*w2_2*d4*sin_th4) +R_TWO*w1_5*w2*(-R_132+R_195*mdhbar*w2*d2*sin_th2-R_SIXTY*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(1310)
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*EXPd2*renso3

case(1311)
S(i,j) = R_FOUR*SQRT2*SQRT1d3*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso4

case(1312)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w2*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case(1313)
S(i,j) = -R_FOUR*ONE_THIRD*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*cos_th2))*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*sin_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*sin_th2+mdhbar2*w2_2*d4*sin_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*sin_th2+R_SEVEN*mdhbar2*w2_2*d4*sin_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*sin_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso6

case(1320)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case(1321)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case(1322)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*sin_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(1330)
S(i,j) = R_TWO*SQRT2*ONE_THIRD*d*sin_th*nubar3*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(1331)
S(i,j) = R_FOUR*ONE_THIRD*nubar3*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*sin_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(1340)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*sin_th2))*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*cos_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(2004)
S(i,j) = HALF*SQRT1d3*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*sin_th2*(-R_TREE+mdhbar*w2*d2*sin_th2)))*EXPd2*renso6

case(2013)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(2022)
S(i,j) = SQRT1d2*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(2031)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(-R_TREE*w1_4+R_TWO*w1_3*w2*(R_EIGHTEEN+mdhbar*w1*d2*cos_th2)+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w1*d2*cos_th2)+R_TWO*w1*w2_3*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)*(R_SIX -mdhbar*w1*d2*cos_th2)+R_TREE*w2_4*(-R_FIVE+R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso4

case(2040)
S(i,j) = HALF*SQRT1d3*nubar*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*cos_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*cos_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*cos_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(2104)
S(i,j) = SQRT1d6*d*sin_th*nubar2*sq_mdhbar*sq_w2*(w1_2-w2_2+R_TWO*mdhbar*w1*w2_2*d2*cos_th2)*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(2113)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(R_TREE-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso5

case(2122)
S(i,j) = d*sin_th*nubar2*sq_mdhbar*sq_w2*(-R_FOUR*w1*w2 +w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*EXPd2*renso6

case(2131)
S(i,j) = -R_TWO*SQRT2*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w1*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*cos_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(2140)
S(i,j) = SQRT1d6*d*sin_th*sq_mdhbar*sq_w2*nubar2*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*cos_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*cos_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*cos_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(2200)
S(i,j) = nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*EXPd2*renso4

case(2201)
S(i,j) = -SQRT2*d*sin_th*nubar2*sq_mdhbar*sq_w1*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2))*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(2202)
S(i,j) = SQRT1d2*nubar*(w1_2-w2_2+R_TWO*mdhbar*w1*w2_2*d2*cos_th2)*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(2203)
S(i,j) = -SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w1*(w1_2-w2_2+R_TWO*mdhbar*w1*w2_2*d2*cos_th2)*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*sin_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(2204)
S(i,j) = HALF*SQRT1d6*nubar*(w1_2-w2_2+R_TWO*mdhbar*w1*w2_2*d2*cos_th2)*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*sin_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*sin_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*sin_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*sin_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*sin_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*sin_th4+R_FOUR*mdhbar3*w2_3*d6*sin_th6))*EXPd2*renso8

case(2210)
S(i,j) = -SQRT2*d*cos_th*nubar2*sq_mdhbar*sq_w1*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*EXPd2*renso4

case(2211)
S(i,j) = d2*sin_2th*nubar3*mdhbar*w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(2212)
S(i,j)= -d*cos_th*nubar2*sq_mdhbar*sq_w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE - mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR +mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(2213)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*cos_th2))*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*sin_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*sin_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(2220)
S(i,j)= SQRT1d2*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE - mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR +mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(2221)
S(i,j) = d*sin_th*nubar2*sq_mdhbar*sq_w1*(w2_4+w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)-R_TWO*w1_3*w2*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)*(R_FOUR-mdhbar*w2*d2*cos_th2)-R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2))*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*EXPd2*renso6

case(2222)
S(i,j) = HALF*nubar*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*(-w2_4-w1_4*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_EIGHTEEN*w1_2*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+R_TWO*w1*w2_3*(R_FOUR+mdhbar*w2*d2*sin_th2)+R_TWO*w1_3*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*sin_th2+R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso8

case(2230)
S(i,j) = -SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w1*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*cos_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(2231)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(w1_2-R_FIVE*w2_2+R_TWO*w1*w2*(-R_TWO+mdhbar*w2*d2*sin_th2))*(-R_FIFTEEN*w2_4+R_TWO*w1_2*w2_2*(R_THIRTY_THREE-R_THIRTEEN*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_3*(R_TWO+mdhbar*w2*d2*cos_th2)+w1_4*(-R_TREE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(2240)
S(i,j) = HALF*SQRT1d6*nubar*(w1_2-w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*sin_th2))*(-R_TREE*w2_6+R_TREE*w1_2*w2_4*(R_THIRTY_FIVE-R_THIRTY_SIX*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_5*(R_EIGHT+mdhbar*w2*d2*cos_th2)+R_TWELVE*mdhbar*w1_3*w2_4*d2*cos_th2*(-R_ONE+R_TWO*mdhbar*w2*d2*cos_th2)+w1_4*w2_2*(-R_HUNDRED_FIVE+R_312*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4)+w1_6*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+R_TWO*w1_5*w2*(-R_TWENTY_FOUR+R_NINETY_NINE*mdhbar*w2*d2*cos_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso8

case(3004)
S(i,j) = ONE_THIRD*SQRT1d2*d*cos_th*nubar2*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(R_TREE*w2_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*sin_th2*(-R_TREE+mdhbar*w2*d2*sin_th2)))*EXPd2*renso6

case(3013)
S(i,j) = R_TWO*SQRT2*ONE_THIRD*d*sin_th*nubar3*sq_mdhbar*sq_w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2))*(R_TREE*w2_3+R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)-R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(-R_TREE-R_SIX*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(3022)
S(i,j) = SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w2*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(3031)
S(i,j) = R_TWO*ONE_THIRD*SQRT2*d*sin_th*nubar3*sq_mdhbar*sq_w1*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*cos_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*cos_th2+R_SEVEN*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*cos_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso5

case(3040)
S(i,j) = ONE_THIRD*SQRT1d2*d*cos_th*nubar2*sq_mdhbar*sq_w2*(-R_NINE*w2_6+R_FIFTEEN*w1_2*w2_4*(R_THIRTEEN-R_TWELVE*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_5*(R_THIRTY_SIX+mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_3*w2_3*(-R_SIXTY+R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_6*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)-R_TREE*w1_4*w2_2*(R_345-R_280*mdhbar*w2*d2*cos_th2+R_THIRTY_SIX*mdhbar2*w2_2*d4*cos_th4) +R_TWO*w1_5*w2*(-R_132+R_195*mdhbar*w2*d2*cos_th2-R_SIXTY*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(3100)
S(i,j) = SQRT2*SQRT1d3*d2*sin_2th*mdhbar*w2*nubar3*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*EXPd2*renso2

case(3101)
S(i,j) = R_FOUR*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*EXPd2*renso3

case(3102)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso4

case(3103)
S(i,j) = -R_TWO*SQRT2*ONE_THIRD*nubar3*sq_mdhbar*sq_w2*d*cos_th*(R_TREE*w1_2+w2_2*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2))*(-R_TREE*w2_3+R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)-R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(R_TREE+mdhbar*w2*d2*sin_th2*(-R_NINE+R_TWO*mdhbar*w2*d2*sin_th2)))*EXPd2*renso5

case(3104)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w2*(R_TREE*w1_2-R_TREE*w2_2+R_TWO*mdhbar*w1*w2_2*d2*cos_th2)*(-R_TWENTY_FOUR*w1*w2_3+R_TREE*w2_4+R_EIGHT*w1_3*w2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2)-R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*sin_th2+R_FOUR*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(3110)
S(i,j) = R_FOUR*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w2*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso3

case(3111)
S(i,j) = R_FOUR*SQRT2*SQRT1d3*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(3112)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*sin_th*nubar3*sq_mdhbar*sq_w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(3113)
S(i,j) = R_FOUR*ONE_THIRD*nubar3*(-R_TREE*w2_3-R_TREE*w1*w2_2*(R_ONE-R_TREE*mdhbar*w2*d2*cos_th2)+R_TREE*w1_3*(R_ONE-mdhbar*w2*d2*cos_th2)+w1_2*w2*(R_TREE+R_SIX*mdhbar*w2*d2*cos_th2-R_TWO*mdhbar2*w2_2*d4*cos_th4))*(R_TREE*w2_3-R_TREE*w1_2*w2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+R_TREE*w1*w2_2*(R_ONE-mdhbar*w2*d2*sin_th2)+w1_3*(-R_TREE+R_NINE*mdhbar*w2*d2*sin_th2-R_TWO*mdhbar2*w2_2*d4*sin_th4))*EXPd2*renso6

case(3120)
S(i,j) = SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(3121)
S(i,j) = R_TWO*SQRT2*SQRT1d3*d*cos_th*nubar3*sq_mdhbar*sq_w2*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso5

case(3122)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w2*(-R_FOUR*w1*w2+w2_2+w1_2*(-R_FIVE+R_TWO*mdhbar*w2*d2*sin_th2))*(-R_TREE*w2_4+R_SIX*w1_2*w2_2*(R_ELEVEN-R_FIVE*mdhbar*w2*d2*cos_th2)+R_TWO*w1*w2_3*(R_EIGHTEEN+mdhbar*w2*d2*cos_th2)+R_TREE*w1_4*(-R_FIVE+R_TWO*mdhbar*w2*d2*cos_th2)+R_TWO*w1_3*w2*(R_SIX-R_THIRTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(3130)
S(i,j) = -R_TWO*ONE_THIRD*SQRT2*d*sin_th*nubar3*sq_mdhbar*sq_w2*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*cos_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*cos_th2+R_SEVEN*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*cos_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso5

case(3131)
S(i,j) =-R_FOUR*ONE_THIRD*nubar3*(w2+w1*(R_ONE-mdhbar*w2*d2*sin_th2))*(R_NINE*w2_5-R_TREE*w1*w2_4*(R_FIVE+R_NINE*mdhbar*w2*d2*cos_th2)+R_SIX*w1_2*w2_3*(-R_FIFTEEN+R_TWELVE*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_5*(R_TREE-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)-R_SIX*w1_3*w2_2*(R_FIFTEEN-R_THIRTY_THREE*mdhbar*w2*d2*cos_th2+R_SEVEN*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2*(-R_FIFTEEN+R_SEVENTY_TWO*mdhbar*w2*d2*cos_th2-R_FOURTY_TWO*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(3140)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w2*(-R_NINE*w2_6+R_FIFTEEN*w1_2*w2_4*(R_THIRTEEN-R_TWELVE*mdhbar*w2*d2*cos_th2)+R_SIX*w1*w2_5*(R_THIRTY_SIX+mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_3*w2_3*(-R_SIXTY+R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_6*(R_TWENTY_SEVEN-R_TWENTY_EIGHT*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)-R_TREE*w1_4*w2_2*(R_345-R_280*mdhbar*w2*d2*cos_th2+R_THIRTY_SIX*mdhbar2*w2_2*d4*cos_th4) +R_TWO*w1_5*w2*(-R_132+R_195*mdhbar*w2*d2*cos_th2-R_SIXTY*mdhbar2*w2_2*d4*cos_th4+R_FOUR*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(4000)
S(i,j) = SQRT1d6*nubar*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*cos_th2*(-R_TREE+mdhbar*w1*d2*cos_th2)))*EXPd2*renso4

case(4001)
S(i,j) = -SQRT1d3*d*sin_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*cos_th2*(-R_TREE+mdhbar*w1*d2*cos_th2)))*EXPd2*renso4

case(4002)
S(i,j) = HALF*SQRT1d3*nubar*(R_TREE*w1_4-R_SIX*w2_2*w1_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*cos_th2*(-R_TREE+mdhbar*w1*d2*cos_th2)))*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso6

case(4003)
S(i,j) = ONE_THIRD*SQRT1d2*d*sin_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w1_4-R_SIX*w1_2*w2_2*(R_ONE-R_TWO*mdhbar*w1*d2*cos_th2)+w2_4*(R_TREE+R_FOUR*mdhbar*w1*d2*cos_th2*(-R_TREE+mdhbar*w1*d2*cos_th2)))*(-R_TREE*w2_2+w1_2*(R_TREE-R_TWO*mdhbar*w2*d2*sin_th2))*EXPd2*renso6

case(4004)
S(i,j) = ONE_TWELFTH*nubar*(R_TREE*w1_4+R_TWELVE*mdhbar*w1_3*w2_2*d2*cos_th2+R_TREE*w2_4-R_TWELVE*mdhbar*w1*w2_4*d2*cos_th2+w1_2*(-R_SIX*w2_2+R_FOUR*mdhbar2*w2_4*d4*cos_th4))*(R_TREE*w2_4-R_SIX*w2_2*w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2)+w1_4*(R_TREE+R_FOUR*mdhbar*w2*d2*sin_th2*(-R_TREE+mdhbar*w2*d2*sin_th2)))*EXPd2*renso8

case(4010)
S(i,j) = SQRT1d3*d*cos_th*nubar2*sq_mdhbar*sq_w1*(-R_TREE*w1_4+R_TWENTY_FOUR*w1_3*w2 +R_SIX*w1_2*w2_2*(R_FIVE-R_TWO*mdhbar*w1*d2*cos_th2)+R_EIGHT*w1*w2_3*(-R_TREE+R_TWO*mdhbar*w1*d2*cos_th2)+w2_4*(-R_TWENTY_SEVEN-R_FOUR*mdhbar*w1*d2*cos_th2*(-R_SEVEN+mdhbar*w1*d2*cos_th2)))*EXPd2*renso4

case(4011)
S(i,j) = HALF*SQRT2*SQRT1d3*d2*sin_2th*nubar3*mdhbar*w1*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*cos_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso4

case(4012)
S(i,j) = -SQRT1d6*d*cos_th*nubar2*sq_mdhbar*sq_w1*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*cos_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(4013)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w1*(R_TREE*w2_2+w1_2*(-R_TREE+R_TWO*mdhbar*w2*d2*sin_th2))*(R_TREE*w1_4+R_TWENTY_SEVEN*w2_4+R_FOUR*w1*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_3*w2*(-R_TWO+mdhbar*w2*d2*cos_th2)+R_TWO*w1_2*w2_2*(-R_FIFTEEN-R_EIGHT*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4))*EXPd2*renso6

case(4020)
S(i,j) = HALF*SQRT1d3*nubar*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*cos_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*cos_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(4021)
S(i,j) =-SQRT1d6*d*sin_th*nubar2*sq_mdhbar*sq_w1*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*cos_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*cos_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(4022)
S(i,j) = HALF*SQRT1d6*nubar*(w2_2-w1_2*(R_ONE-R_TWO*mdhbar*w2*d2*sin_th2))*(R_TREE*w2_6+R_EIGHT*mdhbar*w1_3*w2_4*d2*cos_th2*(R_THIRTY_NINE-R_ELEVEN*mdhbar*w2*d2*cos_th2)-R_TREE*w1_6*(R_ONE-R_TWO*mdhbar*w2*d2*cos_th2)-R_TWELVE*w1*w2_5*(R_FOUR+mdhbar*w2*d2*cos_th2)+R_TWELVE*w1_5*w2*(R_FOUR-R_NINE*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+w1_2*w2_4*(-R_HUNDRED_FIVE+R_198*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2_2*(R_HUNDRED_FIVE-R_TWELVE*mdhbar*w2*d2*cos_th2-R_SIXTY_EIGHT*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso8

case(4030)
S(i,j) = -ONE_THIRD*SQRT1d2*d*cos_th*nubar2*sq_mdhbar*sq_w1*(R_EIGHTY_ONE*w2_6-R_TWELVE*w1*w2_5*(R_TWENTY_TWO+R_SEVEN*mdhbar*w2*d2*cos_th2)+w1_6*(-R_NINE+R_SIX*mdhbar*w2*d2*cos_th2)-R_120*w1_3*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TWELVE*w1_5*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_2*w2_4*(-R_345+R_130*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2_2*(R_195+R_180*mdhbar*w2*d2*cos_th2-R_108*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(4031)
S(i,j) = ONE_SIXTH*d2*sin_2th*nubar3*mdhbar*w1*(R_EIGHTY_ONE*w2_6-R_TWELVE*w1*w2_5*(R_TWENTY_TWO+R_SEVEN*mdhbar*w2*d2*cos_th2)+w1_6*(-R_NINE+R_SIX*mdhbar*w2*d2*cos_th2)-R_120*w1_3*w2_3*(R_SIX-R_SEVEN*mdhbar*w2*d2*cos_th2+mdhbar2*w2_2*d4*cos_th4)+R_TWELVE*w1_5*w2*(R_EIGHTEEN-R_FIFTEEN*mdhbar*w2*d2*cos_th2+R_TWO*mdhbar2*w2_2*d4*cos_th4)+R_TREE*w1_2*w2_4*(-R_345+R_130*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4)+w1_4*w2_2*(R_195+R_180*mdhbar*w2*d2*cos_th2-R_108*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6))*EXPd2*renso6

case(4040)
S(i,j) = ONE_TWELFTH*nubar*( R_NINE*w2_8 -R_THIRTY_SIX*w1*w2_7*(R_EIGHT+mdhbar*w2*d2*cos_th2) +R_TWELVE*w1_3*w2_5*(R_152-R_TWENTY_THREE*mdhbar*w2*d2*cos_th2-R_FOURTY_FOUR*mdhbar2*w2_2*d4*cos_th4) -R_TWELVE*w1_2*w2_6*(R_NINETEEN-R_NINETY_NINE*mdhbar*w2*d2*cos_th2-mdhbar2*w2_2*d4*cos_th4) +R_TREE*w1_8*(R_TREE-R_TWELVE*mdhbar*w2*d2*cos_th2+R_FOUR*mdhbar2*w2_2*d4*cos_th4) +R_FOUR*w1_5*w2_3*(R_456-R_1755*mdhbar*w2*d2*cos_th2+R_840*mdhbar2*w2_2*d4*cos_th4-R_SEVENTY_SIX*mdhbar3*w2_3*d6*cos_th6)-R_TWELVE*w1_7*w2*(R_TWENTY_FOUR-R_NINETY_NINE*mdhbar*w2*d2*cos_th2+R_FOURTY_FOUR*mdhbar2*w2_2*d4*cos_th4-R_FOUR*mdhbar3*w2_3*d6*cos_th6) +R_SIX*w1_4*w2_4*(R_585-R_1170*mdhbar*w2*d2*cos_th2+R_190*mdhbar2*w2_2*d4*cos_th4+R_EIGHT*mdhbar3*w2_3*d6*cos_th6)+R_FOUR*w1_6*w2_2*(-R_FIFTY_SEVEN-R_SIXTY_NINE*mdhbar*w2*d2*cos_th2+R_285*mdhbar2*w2_2*d4*cos_th4-R_SEVENTY_SIX*mdhbar3*w2_3*d6*cos_th6+R_FOUR*mdhbar4*w2_4*d8*cos_th8) )*EXPd2*renso8



        end select

    enddo
   

    enddo

deallocate(basis)

end subroutine Overlap

end module overlap_m
