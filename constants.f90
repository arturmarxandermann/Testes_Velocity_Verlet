module constants_m
use types_m

!	Declarar Unidades

!========= OBS ===========
!ESCREVER CONSTANTE.0D10 OU CONSTANTED10 É A MESMA COISA SE ESTAMOS QUERENDO FAZER UM NUMERO COM EXPOENTE 10

!============   ================

    real*8, parameter :: vel_luz = 3.d10 !esta em cm/s
    real*8, parameter :: PI = 3.14159265359d0
    
    !============   ==================
    
    real*8, parameter  :: half = 0.5d0, zero = 0.0d0, one = 1.0d0, two = 2.0d0, five = 5.0d0, four = 4.0d0, three = 3.0d0, eight = 8.0d0 , six = 6.0d0
    
    real*8, parameter :: EM0 = 9.1091d-28  !massa do eletron gramas
    real*8, parameter :: PM0 = 1.6726d-24 !massa do proton em gramas
    real*8, parameter :: EC0 = 4.803204d-10 !carga do eletron em statcoulomb
    real*8, parameter :: QE = 4.80326E-10  !carga do eletron em esu
    real*8, parameter :: ns = 1E-9
    real*8, parameter :: FS =  1d-15                           !fentossegundo
    real*8, parameter :: as =  1d-18                           !atossegundo
    real*8, parameter :: zs =  1d-21       !zeptossegundo
    
    
    real*8, parameter :: HB_ev = 6.58211915d-16 ! hcortado em eV.s
    real*8, parameter :: kboltz = 8.6173324d-5 !constante de boltzmann em ev/K
    real*8, parameter :: HB_ev_ps = 6.58211915d-4  !hcortado em eV.ps
    real*8, parameter :: hz_to_thz = 1.d-12 !constante para passar as frequencias de Hz para Thz
    real*8, parameter :: thz_to_hz = 1.d12 !constante para passar as frequencias de Thz para Hz
    real*8, parameter :: cteConvert = 5.333333d0
    real*8, parameter :: joule_to_ev = 6.241509343d+18
    real*8, parameter :: ev_to_joule = 1.602176565d-19
    
    real*8, parameter :: mass   = 0.511003d6
    complex*16, parameter :: zi = (0.0d0,1.d0)
    real*8, parameter :: kboltz_cm = 0.695034 !cte de boltzman em cm⁻1/K
    real*8, parameter :: cm_to_ps = 5.3050d0 !transforma o tempo de cm para ps
    
    integer :: grid_size, i, m, j, l, k, n, r, c, alpha, beta, gama, delta, lambda, mu, kappa, sn
    
    REAL*8 , PARAMETER :: R_ZERO = 0.d0
    
    REAL*8 , PARAMETER :: R_ONE = 1.d0
    REAL*8 , PARAMETER :: R_TWO = 2.d0
    REAL*8 , PARAMETER :: R_TREE = 3.d0
    REAL*8 , PARAMETER :: R_FOUR = 4.d0
    REAL*8 , PARAMETER :: R_FIVE = 5.d0
    REAL*8 , PARAMETER :: R_SIX = 6.d0
    REAL*8 , PARAMETER :: R_SEVEN = 7.d0
    REAL*8 , PARAMETER :: R_EIGHT = 8.d0
    REAL*8 , PARAMETER :: R_NINE = 9.d0
    REAL*8 , PARAMETER :: R_ELEVEN = 11.d0
    REAL*8 , PARAMETER :: R_TWELVE = 12.d0
    REAL*8 , PARAMETER :: R_THIRTEEN = 13.d0
    REAL*8 , PARAMETER :: R_FIFTEEN = 15.d0
    REAL*8 , PARAMETER :: R_SIXTEEN = 16.d0
    REAL*8 , PARAMETER :: R_EIGHTEEN = 18.d0
    REAL*8 , PARAMETER :: R_NINETEEN = 19.d0
    REAL*8 , PARAMETER :: R_TWENTY_TWO = 22.d0
    REAL*8 , PARAMETER :: R_TWENTY_THREE = 23.d0
    REAL*8 , PARAMETER :: R_TWENTY_FOUR = 24.d0
    REAL*8 , PARAMETER :: R_TWENTY_SEVEN = 27.d0
    REAL*8 , PARAMETER :: R_TWENTY_EIGHT = 28.d0
    REAL*8 , PARAMETER :: R_THIRTY_THREE = 33.d0
    REAL*8 , PARAMETER :: R_THIRTY_FIVE = 35.d0
    REAL*8 , PARAMETER :: R_THIRTY_SIX = 36.d0
    REAL*8 , PARAMETER :: R_THIRTY_NINE = 39.d0
    REAL*8 , PARAMETER :: R_FOURTY_TWO = 42.d0
    REAL*8 , PARAMETER :: R_FOURTY_FOUR = 44.d0
    REAL*8 , PARAMETER :: R_FOURTY_EIGHT = 48.d0
    REAL*8 , PARAMETER :: R_FIFTY_SEVEN = 57.d0
    REAL*8 , PARAMETER :: R_FIFTY_NINE = 59.d0
    REAL*8 , PARAMETER :: R_SIXTY = 60.d0
    REAL*8 , PARAMETER :: R_SIXTY_FOUR = 64.d0
    REAL*8 , PARAMETER :: R_SIXTY_EIGHT = 68.d0
    REAL*8 , PARAMETER :: R_SIXTY_NINE = 69.d0
    REAL*8 , PARAMETER :: R_SEVENTY_TWO = 72.d0
    REAL*8 , PARAMETER :: R_SEVENTY_SIX = 76.d0
    REAL*8 , PARAMETER :: R_EIGHTY_ONE = 81.d0
    REAL*8 , PARAMETER :: R_NINETY_SIX = 96.d0
    REAL*8 , PARAMETER :: R_NINETY_NINE = 99.d0
    REAL*8 , PARAMETER :: R_108 = 108.d0
    REAL*8 , PARAMETER :: R_HUNDRED_FIVE = 105.d0
    REAL*8 , PARAMETER :: R_120 = 120.d0
    REAL*8 , PARAMETER :: R_132 = 132.d0
    REAL*8 , PARAMETER :: R_152 = 152.d0
    REAL*8 , PARAMETER :: R_180 = 180.d0
    REAL*8 , PARAMETER :: R_130 = 130.d0
    REAL*8 , PARAMETER :: R_190 = 190.d0
    REAL*8 , PARAMETER :: R_195 = 195.d0
    REAL*8 , PARAMETER :: R_198 = 198.d0
    REAL*8 , PARAMETER :: R_280 = 280.d0
    REAL*8 , PARAMETER :: R_285 = 285.d0
    REAL*8 , PARAMETER :: R_312 = 312.d0
    REAL*8 , PARAMETER :: R_345 = 345.d0
    REAL*8 , PARAMETER :: R_456 = 456.d0
    REAL*8 , PARAMETER :: R_585 = 585.d0
    REAL*8 , PARAMETER :: R_840 = 840.d0
    REAL*8 , PARAMETER :: R_1170 = 1170.d0
    REAL*8 , PARAMETER :: R_1755 = 1755.d0
    
    
    REAL*8 , PARAMETER :: k_WH = 0.5d0
    REAL*8 , PARAMETER :: k_WH_hl = 0.5d0
    REAL*8 , PARAMETER :: electronMass = 9.1093897d-31 !Kg
    REAL*8 , PARAMETER :: hbar_ev = 6.5821192815d-16 !eV.s
    REAL*8 , PARAMETER :: hbar = 1.05457172647d-34 !J.s
    !REAL*8 , PARAMETER :: upc1 = 3.90383d-10
    !REAL*8 , PARAMETER :: omg1 = 1.51927d15
    REAL*8 , PARAMETER :: hbar_fac = 1.519267544
    
    
    REAL*8 , PARAMETER :: INVERSEFOURPI = 1.d0 / (4.d0 * PI) 
    !REAL*8 , PARAMETER :: PI = 3.141592653589793238462643383279502884197d0
    REAL*8 , PARAMETER :: SQRT2 = 1.4142135623730950488d0
    
    REAL*8 , PARAMETER :: ONE_EIGHTH = 0.12500000000000000000d0
    REAL*8 , PARAMETER :: ONE_THIRD = 0.33333333333333333333d0
    REAL*8 , PARAMETER :: ONE_FOURTH = 0.25000000000000000000d0
    REAL*8 , PARAMETER :: ONE_SIXTH = 0.16666666666666666667d0
    REAL*8 , PARAMETER :: ONE_TWELFTH = 0.083333333333333333333d0
    REAL*8 , PARAMETER :: ONE_TWENTY_FOURTH = 0.041666666666666666667d0
    REAL*8 , PARAMETER :: ONE_THIRTIETH_SIXTH = 0.027777777777777777778d0
    REAL*8 , PARAMETER :: ONE_FORTY_EIGHTH = 0.020833333333333333333d0
    REAL*8 , PARAMETER :: ONE_NINETY_SIXTH = 0.010416666666666666667d0
    REAL*8 , PARAMETER :: ONE_TREE_HUNDRED_EIGHTY_FOUTH = 0.0026041666666666666667d0
    
    REAL*8 , PARAMETER :: tol = 1.d-18 ! 1.d-15
    REAL*8 , PARAMETER :: SQRT3 = 1.7320508075688772935d0
    REAL*8 , PARAMETER :: SQRT1d2 = 0.70710678118654752440d0
    REAL*8 , PARAMETER :: SQRT1d6 = 0.40824829046386301637d0
    REAL*8 , PARAMETER :: SQRT1d3 = 0.57735026918962576451d0
    REAL*8 , PARAMETER :: SQRT3d2 = 1.2247448713915890491d0
    REAL*8 , PARAMETER :: SQRT2dPI = 0.79788456080286535588d0
    REAL*8 , PARAMETER :: SQRT1d3PI = 0.32573500793527994772d0
    REAL*8 , PARAMETER :: SQRT8dPI = 1.5957691216057307118d0
    REAL*8 , PARAMETER :: SQRT1dPI = 0.56418958354775628695d0
    REAL*8 , PARAMETER :: SQRT32dPI = 3.1915382432114614235d0
    REAL*8 , PARAMETER :: SQRT1d2PI = 0.39894228040143267794d0
    REAL*8 , PARAMETER :: SQRT1d6PI = 0.23032943298089031951d0
    
    INTEGER , PARAMETER :: I_ZERO = 0
    INTEGER , PARAMETER :: I_ONE = 1
    INTEGER , PARAMETER :: I_TWO = 2
    
    INTEGER , PARAMETER :: I_TEN = 10
    INTEGER , PARAMETER :: I_HUNDRED = 100
    INTEGER , PARAMETER :: I_THOUSAND = 1000
    
    COMPLEX*16 , PARAMETER :: C_ONE = ( 1.d0,0.d0 )    ! complex parameter
    COMPLEX*16 , PARAMETER :: C_ZERO = ( 0.d0,0.d0 )    ! complex parameter


end module constants_m
