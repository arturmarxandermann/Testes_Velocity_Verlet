module types_m

type quantum_site
    real*8 :: t
    real*8 :: homo_energy
    real*8 :: lumo_energy
    real*8 :: radius
    real*8 :: radiuszero
    real*8 :: radial_vel 
    real*8 :: mass
    real*8 :: site_type_el
    real*8 :: omega
    real*8 :: omegazero
    real*8 :: posicao_x
    real*8 :: posicao_y
end type quantum_site

type operators
real*8, ALLOCATABLE, dimension(:, :) :: elements
end type operators


type tensor_product_matrices
  real*8, ALLOCATABLE, DIMENSION(:, :) :: elements

end type






end module types_m
