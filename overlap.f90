module overlap_m

use types_m
use constants_m
use parameters_m

contains


subroutine define_sitios(site, site_point)
!subrotina que define todos os sitios da celula OPV
    implicit none

    ! args
    type(quantum_site), target, allocatable, intent(out) :: site(:,:)
    type(obj_pointer), allocatable, intent(out) :: site_point(:)
    
    ! local
    real*8, parameter :: scale_nanometro = 1.d-9
    
    
    allocate(site(nr, nc))
    allocate(site_point(nsites)) !parametrizando os sitios com um único índice
    
    
    site%t       = siteCoupling
    site%mass    = siteMass  
    site%radius0 = raioZero
    site%radius  = site%radius0
    site%vel     = 0.d0 
    site%omega   = ( 2.d0 * hbar / ( me * (site%radius)**2.0 ) ) * hz_to_thz   !largura a meia altura
    site%omega0  = site%omega * thz_to_hz 
    site%V0      = 0.d0

    
    do j = 1, nc
      do i = 1, nr
        !PRIMEIRO SITIO LOCALIZADO EM 0,0
        if (i == 1 .AND. j == 1) then
          site(i, j)%xPos = 0.d0
          site(i, j)%yPos = 0.d0
    
        else
          site(i, j)%xPos = site(i , j-1)%xPos + site(i, j-1)%radius + site(i, j)%radius &
            + distance_molecule
          site(i, j)%yPos = 0.d0
    
        endif
      enddo
    enddo
    
    do j = 1, nc
      do i = 1, nr
        print*, "tamanho molecula", i, j, "em nm é", 2.d0 * site(i, j)%radius * 1.d9
        enddo
    enddo
    
    do j = 1, nc
      do i = 1, nr
        print*, "omega da molecula", i, j, "em Thz",  site(i, j)%omega
        enddo
    enddo
    
    do j = 1, nc
      do i = 1, nr
        print*, "V type da molecula el", i, j, "em ev é",  site(i, j)%V0
      enddo
    enddo
    
    
    do j = 1, nc
      do i = 1, nr
        print*, "HBAR * omega da molecula:", i, j, "em ev é",  HB_ev_ps*site(i, j)%omega
      enddo
    enddo
    
    
    print*, "distancia entre os sitios utilizada:", distance_molecule
    
    k = 1
    do j = 1, nc
      do i = 1, nr
        site_point(k)%np => site(i, j) !site_pointer
        k = k + 1
      enddo
    enddo



end subroutine define_sitios

subroutine Basis_Builder_Blocks
    implicit none

    !local
    integer :: l, k

    k = 0
    l = 0     

    do j = 1, nsites
      do i = 1, nsites        
          basis(i, j)%rmin = k*ns_el+1
          basis(i, j)%rmax = k*ns_el+ns_el
          basis(i, j)%cmin = l*ns_el+1
          basis(i, j)%cmax = l*ns_el+ns_el

          k = k + 1
      enddo
     k = 0 
     l = l + 1 
    enddo

end subroutine Basis_Builder_Blocks  


subroutine Basis_Builder_hMtx
    implicit none


    !local
    real*8, allocatable :: SMtx(:,:) 
    integer             :: i, j, k


    do j = 1, nsites-1
     do i = j + 1, nsites
      call Overlap(i, j, SMtx)
      basis(i, j)%hMtx(:,:) = SMtx(:,:) * site_point(i)%np%t 
      deallocate(SMtx) 
      enddo
    enddo
  
    do i = 1, nsites
      do k = 1, ns_el 
        basis(i, i)%hMtx(k, k) = site_point(i)%np%V0 + HB_ev_ps * site_point(i)%np%omega * Qn_erg(k)
      enddo
    enddo

end subroutine Basis_Builder_hMtx


subroutine Basis_Builder_DerMtx
    implicit none


    !local
    real*8, allocatable :: DMtx(:,:) 
    integer             :: i, j, k


    do j = 1, nsites-1
     do i = j + 1, nsites
      call DerivativeOverlap(i, j, DMtx)
      basis(i, j)%DerMtx(:,:) = DMtx
      enddo
    enddo


end subroutine Basis_Builder_DerMtx



subroutine Overlap(s1, s2, SMtx)
    implicit none

    ! args
    integer,             intent(in)  :: s1, s2
    real*8, allocatable, intent(out) :: SMtx(:,:)
    
    ! local
    real*8  :: x1, y1, w1, x2, y2, w2, x12, y12
    real*8  :: w1_2, w2_2, w1_05, w2_05, w1_15, w2_15
    real*8  :: d, d2, d3
    real*8  :: sinth, costh, sinth2, costh2, sinth3, costh3
    real*8  :: sqmhbar, mhbar, mhbar2, sqhbarm, hbarm, hbarm_15
    real*8  :: wproduct, wproduct2, wproduct3
    real*8  :: wden1, wden2, wden3, wden4, wden5, wden6, wden7
    real*8  :: Omega1PlusOmega2, SQOmega1PlusOmega2
    real*8  :: me2, hbar2
    integer :: Q_numbers
    real*8  :: EXPd2
    
    
    allocate(SMtx(ns_el, ns_el), source = 0.d0 )
    
    
    x1 = 0.d0 !site(l1, c1)%xPos
    y1 = 0.d0 !site(l1, c1)%yPos
    w1 = site_point(s1)%np%omega*thz_to_hz 

    x2 = site_point(s2)%np%xPos - site_point(s1)%np%xPos
    y2 = site_point(s2)%np%yPos - site_point(s1)%np%yPos
    w2 = site_point(s2)%np%omega*thz_to_hz 
    x12 = x1 - x2
    y12 = y1 - y2

    w1_05 = Sqrt(w1)
    w1_2  = w1 * w1
    w1_15 = w1_05 * w1

    w2_05 = Sqrt(w2)  
    w2_2  = w2 * w2 
    w2_15 = w2_05 * w2
        
    d   = Sqrt(x12*x12 + y12*y12)
    d2  = d * d
    d3  = d2 * d
    
    
    sinth = y12 / d
    costh = x12 / d
    
    sinth2 = sinth * sinth
    sinth3 = sinth2 * sinth
    
    costh2 = costh * costh
    costh3 = costh2 * costh
    
    hbar2 = hbar * hbar
    me2 = me * me

    sqmhbar = sqrt( (me/hbar) )
    mhbar   = sqmhbar * sqmhbar
    mhbar2  = mhbar * mhbar
    
    sqhbarm = sqrt( (hbar/me) )
    hbarm   = sqhbarm * sqhbarm
    hbarm_15 = sqhbarm * hbarm
    
    wproduct  = sqrt(w1*w2) 
    wproduct2 = wproduct  * wproduct 
    wproduct3 = wproduct2 * wproduct 
    

    Omega1PlusOmega2 = w1 + w2
    SQOmega1PlusOmega2 = Sqrt(Omega1PlusOmega2) 


    wden1 = 1.d0 /  SQOmega1PlusOmega2
    wden2 = wden1 * wden1 
    wden3 = wden2 * wden1 
    wden4 = wden3 * wden1
    wden5 = wden4 * wden1 
    wden6 = wden5 * wden1 
    wden7 = wden6 * wden1 
    
    
    EXPd2 = EXP( - HALF * d2 * mhbar * wden2 * wproduct2)
    
    do j = 1 , ns_el
      do i = 1 , ns_el
        Q_numbers = 100*Qn(i) + 1*Qn(j)
        select case (Q_numbers)
          case(0000)
          SMtx(i, j) = 2.d0*EXPd2*wden2*wproduct
    
          case(0001)
          SMtx(i, j) = (2.d0*SQRT2*d*EXPd2*sinth*Sqrt(mhbar*w2)*(-(hbar*SQOmega1PlusOmega2) + &
                   me*sqhbarm*w2*Sqrt(hbarm*wden2))*wden3*wproduct)/hbar
          
          case(0010)
          SMtx(i, j) = -2.d0*SQRT2*costh*d*EXPd2*sqmhbar*w1_15*w2*wden4*multFactorOvlp
          
          
          case(0100) 
          SMtx(i, j) = ((2.d0*SQRT2)*d*EXPd2*sinth*w1_05*Sqrt(mhbar*w1)*w2_15*Sqrt(hbarm*wden2)*wden3)/sqhbarm
          
          case(0101)
          SMtx(i, j) = (4.d0*EXPd2*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*( (d2*me*Sqrt(hbar*me)*sinth2*sqhbarm*w2_2) + ( hbar2 * &
                 Omega1PlusOmega2) - (d2*hbar*me*sinth2*w2*Omega1PlusOmega2))*Sqrt(hbarm*wden2)*wden5*wproduct)/( hbarm_15 * me2 )
          
          
          case(0110)
          SMtx(i, j) = -4.d0*costh*d2*EXPd2*sinth*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*wden6*wproduct3
          
          case(1000)
          SMtx(i, j) = 2.d0*SQRT2*costh*d*EXPd2*sqmhbar*w1*w2_15*wden4*multFactorOvlp 
          
          
          case(1001) 
          SMtx(i, j) = (-4.d0*costh*d2*EXPd2*sinth*w1_05*Sqrt(mhbar*w1)*w2_15*Sqrt(mhbar*w2)*wden7* &
                    (w1_2*Sqrt(hbarm*wden2) + w2*(-sqhbarm*SQOmega1PlusOmega2 + w2*Sqrt(hbarm*wden2)) + & 
                    2.d0*Sqrt(hbarm*wden2)*wproduct2))/sqhbarm
          
          case(1010)
          SMtx(i, j) = (4.d0*EXPd2*Sqrt(mhbar*w1)*Sqrt(mhbar*w2)*Sqrt(hbarm*wden2)*wden5*wproduct*(hbar*sqhbarm*Omega1PlusOmega2 - &
                     costh2*d2*Sqrt(hbar*me)*wproduct2))/hbar 
          
          
        end select
    
      enddo
    enddo

end subroutine Overlap

subroutine DerivativeOverlap(s1, s2, DerOvlp)        
    implicit none
    !SUBROTINA PARA CALCULAR OS ELEMENTOS <sitio(l1, c1, x-x0, y-y0) | d/dR sitio(l2, c2, x, y) >

    ! args
    integer,             intent(in)  :: s1, s2
    real*8, allocatable, intent(out) :: DerOvlp(:, :)
    
    ! local
    real*8    :: x1, y1, w1, x2, y2, w2, x12, y12
    real*8    :: w1_2, w1_3, w1_15, w2_2, w2_3 
    real*8    :: d, d2, d3, d4
    real*8    :: sinth, costh, sinth2, costh2, sinth3, costh3
    real*8    :: sqmhbar, mhbar, sqhbarm, hbarm
    real*8    :: wproduct, wproduct2, wproduct3, wproduct4, wproduct5 
    real*8    :: wden1, wden2, wden6, wden8, wden10
    real*8    :: omega1MinusOmega2, omega1PlusOmega2
    real*8    :: omegaMinusPlus, omegaMinusPlus2
    real*8    :: w2Sw1half, w2Sw1half2, w1Sw2half, w1Sw2half2
    integer   :: Q_numbers
    real*8    :: EXPd2, DerTerm
    real*8    :: me2, hbar2, hbar15, sin2th, cos2th
   
    allocate(DerOvlp(ns_el, ns_el), source = R_zero )
    
    
    x1 = 0.d0 
    y1 = 0.d0 
    w1 = site_point(s1)%np%omega * thz_to_hz
    
    x2 = site_point(s2)%np%xPos - site_point(s1)%np%xPos
    y2 = site_point(s2)%np%yPos - site_point(s1)%np%yPos
    w2 = site_point(s2)%np%omega * thz_to_hz

    w1_2 = w1 * w1
    w1_15 = Sqrt(w1) * w1
    w1_3 = w1_2 * w1

    w2_2 = w2 * w2
    w2_3 = w2_2 * w2
    
    x12 = x1 - x2
    y12 = y1 - y2
    
    d = sqrt(x12 * x12 + y12 * y12)
    d2 = d*d
    d3 = d2*d
    d4 = d2*d2
    
    sinth = y12 / d
    costh = x12 / d
    
    
    
    DerTerm = - ((4.d0 * hbar) / (me * (site_point(s2)%np%radius)**3.0)  )
    
    sinth2 = sinth * sinth
    sinth3 = sinth2 * sinth 
    
    costh2 = costh * costh
    costh3 = costh2 * costh
    
    sin2th = 2.d0 * sinth * costh
    cos2th = costh2 - sinth2
    
    sqmhbar = sqrt( (me/hbar) )
    mhbar = sqmhbar * sqmhbar
    
    sqhbarm = sqrt( (hbar/me) )
    hbarm = sqhbarm * sqhbarm
    
    wproduct = sqrt(w1*w2)
    wproduct2 = wproduct  * wproduct
    wproduct3 = wproduct2 * wproduct
    wproduct4 = wproduct2 * wproduct2
    wproduct5 = wproduct4 * wproduct
    
    omega1MinusOmega2 = w1 - w2
    omega1PlusOmega2  = w1 + w2
    omegaMinusPlus    = omega1MinusOmega2 * omega1PlusOmega2
    omegaMinusPlus2   = omega1MinusOmega2 * ( omega1PlusOmega2 * omega1PlusOmega2 )
    
    w1Sw2half2 = (w1/w2)
    w1Sw2half  =  sqrt(w1Sw2half2)
    
    w2Sw1half2 = (w2/w1)
    w2Sw1half  = sqrt(w2Sw1half2)
    
    wden1 = 1.d0 / ( Sqrt(omega1PlusOmega2) )
    wden2 = wden1 * wden1
    wden6 = wden2 * wden2 * wden2
    wden8 = wden6 * wden2
    wden10 = wden8 * wden2
    
    hbar15 = Sqrt(hbar) * hbar
    hbar2 = hbar * hbar
    me2 = me * me
    
    EXPd2 = EXP( - HALF * d2 * mhbar * wden2 * wproduct2)
    
    do j = 1 , ns_el
      do i = 1 , ns_el
        Q_numbers = 100*Qn(i) + 1*Qn(j)
        select case (Q_numbers)
    
        case(0000)
          DerOvlp(i, j) = -((EXPd2*w1Sw2half*(d2*me*w1_2*w2 + hbar*(-w1_2 + w2_2))*wden6)/hbar)
          
        case(0001)
          DerOvlp(i, j) = ( SQRT2*d*EXPd2*me*sinth*w1_15*(d2*Sqrt(hbar*me)*w1_2*w2+2.d0*hbar*sqhbarm*& 
                            (-w1_2 + w2_2))*wden8)/hbar2
          
        case(0010)
          DerOvlp(i, j) = SQRT2*multFactorDer*costh*d*EXPd2*Sqrt(me)*(w1/hbar)**1.5d0*(-(d2*me*w1_2*w2) + &
                            2.d0*hbar*(w1_2 - w2_2))*wden8
          
        case(0100) 
          DerOvlp(i, j) = -((SQRT2*d*EXPd2*me*sinth*w1*Sqrt(w2)*(-(hbar*omega1PlusOmega2*sqhbarm*(3.d0*w1 - w2)) + & 
                             d2*Sqrt(hbar*me)*w1_2*w2)*wden8)/hbar2 )        
          
        case(0101)
          DerOvlp(i, j) = (2.d0*EXPd2*w1*wden10*(2.d0*hbar2*omegaMinusPlus2 + d4*me2*sinth2*w1_3*w2_2 + &
                            d2*hbar*me*omega1PlusOmega2*(-3.d0*w1 + cos2th*(2.d0*w1 - w2) + w2)*wproduct2))/hbar2
          
          
        case(0110)
          DerOvlp(i, j) = (d2*EXPd2*me*sin2th*w1_2*w2*(-2.d0*hbar*omega1PlusOmega2*(2.d0*w1 - w2) + &
                            d2*me*w1_2*w2)*wden10)/hbar2
          
          
        case(1000)
          DerOvlp(i, j) = -((multFactorDer*SQRT2*EXPd2*costh*d*w1*Sqrt(me*w2)*(-(hbar*omega1PlusOmega2*(3.d0*w1 - w2)) + &
                             d2*me*w1_2*w2)*wden8)/(hbar15)) 
          
        case(1001) 
          DerOvlp(i, j) = (d2*EXPd2*me*sin2th*w1_2*w2*(-2.d0*hbar*omega1PlusOmega2*(2.d0*w1 - w2) + &
                            d2*me*w1_2*w2)*wden10)/hbar2
          
        case(1010)
          DerOvlp(i, j) = (2.d0*w1*EXPd2*wden10*(costh2*d4*me2*w1_3*w2_2 + hbar*omega1PlusOmega2*(2.d0*hbar*(w1_2 - w2_2) + & 
                            d2*me*(-3.d0*w1 + w2)*wproduct2 + cos2th*d2*me*(-2.d0*w1 + w2)*wproduct2)))/hbar2 
    
        end select
    
    
      enddo
    enddo
    
    DerOvlp = DerOvlp * DerTerm     
end subroutine DerivativeOverlap






  
end module overlap_m
