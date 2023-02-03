module allocation

    use constants
    use types
    use traits

implicit none

private

public :: alloc, sensitivity!, f

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
        !variable for bisec loop
        integer(i_4) :: x, i
        real(r_8) :: tmp
        real(r_8) :: tmp2
        real(r_8) :: sens, sens_res

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

        ! real(r_8) :: height
        ! real(r_8) :: lm1
        ! real(r_8) :: lminc_min
        ! real(r_8) :: rminc_min

        !initializing variables
        leaf_in_ind = 0.0D0
        root_in_ind = 0.0D0
        sap_in_ind = 0.0D0
        heart_in_ind = 0.0D0

        leaf_out = 0.0D0
        root_out = 0.0D0
        sap_out  = 0.0D0
        heart_out = 0.0D0


        !!!!!CHECK IF DWOOD IS IN THE CORRECT UNIT


        !carbon (gC) in compartments considering the density (ind/m2)
        leaf_in_ind = (leaf_in/dens_in)*1.D3
        root_in_ind = (root_in/dens_in)*1.D3
        sap_in_ind = (sap_in/dens_in)*1.D3 
        heart_in_ind = (heart_in/dens_in)*1.D3

        ! call functions to logic 
        ! height = alt_height(wood_ind, woodens)
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
            
           if (x.eq.200) exit 
        
           if ((abs(tmp - tmp2)).le.sens) then
            print*, 'sensitivity attained'
            exit 
           endif 
           
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

    function sensitivity(tmp, tmp2) result(tmp_sens)
        
        real(r_8), intent(in) :: tmp
        real(r_8), intent(in) :: tmp2

        real(r_8) :: tmp_sens

        tmp_sens = tmp+tmp2
    
    end function sensitivity
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

end module allocation