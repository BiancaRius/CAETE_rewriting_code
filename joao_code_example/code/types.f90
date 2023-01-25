module types
   implicit none
   ! FOR THE GNU FORTRAN COMPILER
   integer,parameter,public :: l_1 = 2  ! standard Logical type
   integer,parameter,public :: i_2 = 2  ! 16 bits integer
   integer,parameter,public :: i_4 = 4  ! 32 bits integer
   integer,parameter,public :: r_4 = 4  ! 32 bits float
   integer,parameter,public :: r_8 = 8  ! 64 bits float

end module types

module allometric_params
    use types
    implicit none

    real(r_8),parameter,public :: grid_area = 1000.0D0 !m2
    real(r_8),parameter,public :: pi   =  3.14159265D0
    real(r_8),parameter,public :: xacc =  0.1D0     !x-axis precision threshold for the allocation solution
    real(r_8),parameter,public :: yacc =  1.0D-10  !y-axis precision threshold for the allocation solution

    integer(i_4),parameter,public :: nseg = 20 ! segment number to bisection loop
    integer(i_4),parameter,public :: time = 1000
    integer(i_4),parameter,public :: ntl=365

    real(r_8),parameter,public :: allom1 = 100.D0
    real(r_8),parameter,public :: allom2 = 40.0D0
    real(r_8),parameter,public :: allom3 = 0.5D0
    real(r_8),parameter,public :: latosa = 8000.0D0
    real(r_8),parameter,public :: reinickerp = 1.6D0
    real(r_8),parameter,public :: ltor = 0.77302587552347657D0 !leaf:root from Philip

end module allometric_params