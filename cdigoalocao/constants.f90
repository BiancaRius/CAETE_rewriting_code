module types
    implicit none
 
    ! FOR THE GNU FORTRAN COMPILER
    integer,parameter,public :: l_1 = 2  ! standart Logical type
    integer,parameter,public :: i_2 = 2  ! 16 bits integer
    integer,parameter,public :: i_4 = 4  ! 32 bits integer
    integer,parameter,public :: r_4 = 4  ! 32 bits float
    integer,parameter,public :: r_8 = 8  ! 64 bits float
 
end module types

module constants
    use types
    implicit none

    !==============================!
    !========= CONSTANTS ==========!
    !==============================!
    real(r_8), parameter :: klatosa = 6000.0
    real(r_8), parameter :: ltor = 0.77302587552347657
    real(r_8), parameter :: k_allom1 = 100.0 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8), parameter :: k_allom2 = 40.0
    real(r_8), parameter :: k_allom3 = 0.50
    real(r_8), parameter :: tol = 0.0000001
    real(r_8), parameter :: pi = 3.1415926536
    real(r_8), parameter :: reinickerp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8), parameter :: nseg = 20
    real(r_8), parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
    real(r_8), parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

    !==============================!

end module constants

module traits
    use types
    implicit none

    real(r_8), parameter :: sla = 0.023
    real(r_8), parameter :: dwood = 0.74

end module traits