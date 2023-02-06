module allocation

    use constants
    use types
    use traits

implicit none

private

public :: alloc,& !(s) calculates carbon pools output from NPP and carbon from preivous step
          height_calc,& !(f)calculates height
          leaf_req_calc,& !(f)leaf mass requeriment to satisfy allometry
          leaf_inc_min_calc,& !(f) minimum leaf increment to satisfy allocation equations
          root_inc_min_calc

contains

    subroutine alloc(leaf_in, root_in, sap_in, heart_in, bminc_in, dens_in,&
        leaf_out, root_out, sap_out, heart_out)

        !VARIABLE INPUTS
        !carbon inputs (kgC/m2)
        real(r_8), intent(in) :: leaf_in
        real(r_8), intent(in) :: root_in
        real(r_8), intent(in) :: sap_in
        real(r_8), intent(in) :: heart_in


        real(r_8), intent(in) :: dens_in !ind/m2 initial density of individuals

        real(r_8), intent(in) :: bminc_in !carbon (NPP) available to be allocated
                                          !basically NPPt - NPPt-1. NPP accumulated in the year/month/day
                                          !gc/ind/time_step
       
        !VARIABLES OUTPUTS 
        !carbon inputs (kgC/m2)
        real(r_8), intent(out) :: leaf_out
        real(r_8), intent(out) :: root_out
        real(r_8), intent(out) :: sap_out
        real(r_8), intent(out) :: heart_out

        !INTERNAL VARIABLES
        real(r_8) :: wood_in_ind !(gC) total wood - sum of heart and sap

        !variable for bisec loop
        integer(i_4) :: x, i
        real(r_8) :: tmp
        real(r_8) :: tmp2
        real(r_8) :: sens

        !carbon (gC) in compartments considering the density (ind/m2)
        real(r_8) :: leaf_in_ind
        real(r_8) :: root_in_ind
        real(r_8) :: sap_in_ind
        real(r_8) :: heart_in_ind

        !carbon increment (gC) in compartments 
        real(r_8) :: leaf_inc
        real(r_8) :: root_inc
        real(r_8) :: sap_inc
        real(r_8) :: heart_inc

        !Functions to the logic
        !dwood !!####****!! ATTENTION: dwood is already transformed to gC/m3 at constants.f90
        real(r_8) :: height
        real(r_8) :: leaf_req
        real(r_8) :: leaf_inc_min
        real(r_8) :: root_inc_min


        !initializing variables
        leaf_in_ind = 0.0D0
        root_in_ind = 0.0D0
        sap_in_ind = 0.0D0
        heart_in_ind = 0.0D0
        wood_in_ind = 0.0D0
        leaf_inc_min = 0.0D0
        root_inc_min = 0.0D0

        leaf_out = 0.0D0
        root_out = 0.0D0
        sap_out  = 0.0D0
        heart_out = 0.0D0

        height = 0.0D0
        leaf_req = 0.0D0

        !carbon (gC) in compartments considering the density (ind/m2)
        leaf_in_ind = (leaf_in/dens_in)!*1.D3
        root_in_ind = (root_in/dens_in)!*1.D3
        sap_in_ind = (sap_in/dens_in)!*1.D3 
        heart_in_ind = (heart_in/dens_in)!*1.D3
        wood_in_ind = sap_in_ind + heart_in_ind

        ! call functions to logic
        height = height_calc(wood_in_ind)
        print*, 'height',height

        leaf_req = leaf_req_calc(sap_in_ind, height)

        print*, 'leaf_red',leaf_req

        leaf_inc_min = leaf_inc_min_calc(leaf_req, leaf_in_ind)
        print*,'leaf inc min', leaf_inc_min,'leaf_in_ind', leaf_in_ind

        root_inc_min = root_inc_min_calc(leaf_req, root_in_ind)
        print*, 'root_inc_min', root_inc_min
        !lminc_min = leaf_min(lm1,leaf_ind)

        ! rminc_min = root_min(lm1,root_ind)
        
        ! lm1 = leaf_mass(sap_ind,height,sla,woodens)
        ! lminc_min = leaf_min(lm1,leaf_ind)
        ! rminc_min = root_min(lm1,root_ind)

        !define bisection method
            !if f(x1) * f(x2) gt 0 then return -2
        !função f - root_bisec
        
        !________________
        !sensitivity test
        x = 0
        tmp  = 0.2
        tmp2 = 1
        sens = 1.e-3

        do i = 1, 200 
            
        !    if (x.eq.200) exit 
        
           if ((abs(tmp - tmp2)).le.sens) then
            print*, 'sensitivity attained'
            exit 
           endif 
           
           !IFS normal/abnormal allocation
           !sensitivity of carbon (yes)
           ! put scape infinity loop (print)
           tmp = tmp2 
           
           x = x + 1
           
           tmp2 = (tmp2 + 4)/2

        enddo

        print*, x
        !_________________



        !Increment C compartments - OUTPUT FINAL (kgC/m²)
        !!############PROVISÓRIO

        leaf_inc = 3.
        root_inc = 3.
        sap_inc = 3.
        heart_inc = 3.

        leaf_out = ((leaf_in_ind + leaf_inc)*dens_in)/1.D3
        root_out = ((root_in_ind + root_inc)*dens_in)/1.D3
        sap_out  = ((sap_in_ind + sap_inc)*dens_in)/1.D3
        heart_out = ((heart_in_ind + heart_inc)*dens_in)/1.D3


    end subroutine alloc

    function height_calc (wood_in_ind) result (height)
        
        real(r_8), intent(in) :: wood_in_ind !gC/ind - total wood (sap + heart) carbon stock
        
        !Trait
        !dwood - wood density
        
        real(r_8) :: height !m - output

        !variable internal
        real(r_8) :: diameter 
        
        

        !Calculo diameter (necessary to height)
        diameter = ((4*wood_in_ind)/(dwood)*pi*k_allom2)**(1/(2+k_allom3))

        !Height 
        height = k_allom2*(diameter**k_allom3)


    end function

    function leaf_req_calc (sap_in_ind, height)  result (leaf_req)
    
        real(r_8), intent(in) :: sap_in_ind !gC - sapwood input
        real(r_8), intent(in) :: height !m
       
        real(r_8) :: leaf_req !gC - output- leaf mass requeriment to satisfy allometry 

        !Trait
        !dwood - wood density (gc/m3) - already transformed in constants.f90
        !sla - specific leaf area

        leaf_req = (klatosa * sap_in_ind / ((dwood) * height * (sla)))

    end function leaf_req_calc

    function leaf_inc_min_calc (leaf_req, leaf_in_ind) result (leaf_inc_min)
        use types
        use constants    

        real(r_8), intent(in) :: leaf_req !gC leaf mass requeriment to satisfy allometry
        real(r_8), intent(in) :: leaf_in_ind !gC leaf input
        
        real(r_8) :: leaf_inc_min !gC -output- minimum leaf increment to satisfy allocation equations

        leaf_inc_min = leaf_req - leaf_in_ind

    end function leaf_inc_min_calc

    function root_inc_min_calc (leaf_req, root_in_ind) result (root_inc_min) !ROOT MASS MINIMO
        use types
        use constants

        !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
        !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

        real(r_8), intent(in) :: root_in_ind !gC root input
        real(r_8), intent(in) :: leaf_req  !gC leaf mass requeriment to satisfy allometry
        
        real(r_8) :: root_inc_min !gC -output- minimum root increment to satisfy allocation equations

        ! real(r_8) :: rm_ind
        ! rm_ind = (rm*1.D3)/nind

        root_inc_min = (leaf_req / ltor - root_in_ind)

    end function root_inc_min_calc

end module allocation

! function f(leaf_prev_alloc, sap_prev_alloc, heart_prev_alloc, root_prev_alloc,&
    ! bminc, sla, dwood, xis) result (root_soluction)
! 
    !carbon previous allocation process
    ! in kgC/m2 but use in gC/ind [transformation below]
    ! real(r_8), intent(in) :: leaf_prev_alloc 
    ! real(r_8), intent(in) :: sap_prev_alloc
    ! real(r_8), intent(in) :: heart_prev_alloc
    ! real(r_8), intent(in) :: root_prev_alloc
! 
    ! real(r_8), intent(in) :: bminc !carbon (NPP) available to be allocated
                                !    basically NPPt - NPPt-1. NPP accumulated in the year/month/day
                                !    gc/ind/time_step
    ! Traits
    ! real(r_8), intent(in) :: sla !COLOCAR UNIDADE
    ! real(r_8), intent(in) :: dwood !COLOCAR UNIDADE
    ! 
    ! real(r_8), intent(in) :: xis !leafmass allocation amount as input !!DE ONDE VEM ESSE VALOR?
! 
! 
! end function f
!initialize sap, leaf, heart, root

!use constants
!verify constants values

!definir função distribution !equações tau1, tau 2....

!definir função bisection method
    !loop tolerancia(use do while)

!função find deltaL, deltaR, deltaS

!definir função f(x) - achar equivalente
