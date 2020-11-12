module types_m

type quantum_site
    real*8 :: t
    real*8 :: radius
    real*8 :: radius0
    real*8 :: vel 
    real*8 :: mass
    real*8 :: V0
    real*8 :: omega
    real*8 :: omega0 
    real*8 :: xPos
    real*8 :: yPos
end type quantum_site

type obj_pointer
  type(quantum_site), pointer :: np 
end type obj_pointer

type BasisBuild
    integer :: rmin
    integer :: rmax
    integer :: cmin
    integer :: cmax 
    real*8, allocatable  :: hMtx(:,:)
    real*8, allocatable  :: DerMtx(:,:) 
    complex*16, allocatable :: TMtx(:,:) 
end type BasisBuild

type operators
    real*8, ALLOCATABLE, dimension(:,:) :: elements
end type operators


type tensor_product_matrices
    real*8, ALLOCATABLE, DIMENSION(:,:) :: elements
end type

end module types_m
