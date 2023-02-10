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
    real(r_8), parameter, public :: klatosa =  6000.0
    real(r_8), parameter, public :: ltor = 0.77302587552347657
    real(r_8), parameter, public :: k_allom1 = 100.0 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8), parameter, public :: k_allom2 = 36!40.0
    real(r_8), parameter, public :: k_allom3 = 0.22!0.50
    real(r_8), parameter, public :: tol = 0.0000001
    real(r_8), parameter, public :: pi = 3.1415926536
    real(r_8), parameter, public :: reinickerp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8), parameter, public :: nseg = 20
    real(r_8), parameter, public :: xacc =  0.1     !x-axis precision threshold for the allocation solution
    real(r_8), parameter, public :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

    !==============================!

end module constants

module traits
    use types
    implicit none

    real(r_8), parameter, public :: sla = 15.36 !0.023
    real(r_8), parameter, public :: dwood = 200.!0.74*1.D3!200. !0.74*1.D6

end module traits