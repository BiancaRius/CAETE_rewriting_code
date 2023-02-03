module allocation_reform

    implicit none

    private

    public ::               &
    leaf_mass              ,& ! (f), leaf mass requeriment to satisfy allometry 
    leaf_min               ,& ! (f), leaf mass minimo
    root_min               ,& ! (f), root mass minimo
    root_bisec             ,& ! (f), raizes da bisection
    alt_height             ,& ! (f), height (in m.)
    alloc                     ! (s), allocation module


contains

    function leaf_mass (sm, height, sla, dwood)  result (lm1)
        use types
        use constants
    
        real(r_8), intent(in) :: sm !in kgC/m2 but use in gC [transformation below]
        real(r_8), intent(in) :: height
        real(r_8), intent(in) :: sla, dwood !traits
        real(r_8) :: lm1 !individual

        lm1 = (klatosa*sm/((dwood*1000000)*height*sla))

    end function leaf_mass

    function leaf_min (lm1, lm) result (lminc_min)
        use types
        use constants    

        real(r_8), intent(in) :: lm !in kgC/m2 but use in gC [transformation below]
        real(r_8), intent(in) :: lm1
        real(r_8) :: lminc_min !individual 

        ! real(r_8) :: lm_ind
        ! lm_ind = (lm*1.D3)/nind

        lminc_min = (lm1 - lm)

    end function leaf_min

    function root_min (lm1, rm) result (rminc_min) !ROOT MASS MINIMO
        use types
        use constants

        !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
        !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

        real(r_8), intent(in) :: rm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: lm1
        real(r_8) :: rminc_min !individual 

        ! real(r_8) :: rm_ind
        ! rm_ind = (rm*1.D3)/nind

        rminc_min = (lm1 / ltor - rm)

    end function root_min

    function root_bisec (lm, sm, hm, rm, bminc, sla, dwood, xis) result (root) 
        use types
        use constants
        
        !RAIZES DA BISECTION = PONTOS NA RETA QUE A CURVA DA BISECTION CORTA 
        
        real(r_8), intent(in) :: lm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: sm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: hm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: rm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: xis !leafmass allocation amount as input
        real(r_8), intent(in) :: bminc
        real(r_8), intent(in) :: sla, dwood !traits
        
        real(r_8) :: pi4 !internal
        real(r_8) :: a1  !internal
        real(r_8) :: a2  !internal
        real(r_8) :: a3  !internal
        real(r_8) :: root !output


        ! Calculos:
        pi4 = (pi/4)
        a1 = (2./k_allom3)
        a2 = (1. + k_allom1)
        a3 = (k_allom2**a1)
        
        
        root = a3 * ((sm + bminc - xis - ((lm + xis) / ltor) + rm + hm) / dwood) / pi4 - &
                    &((sm + bminc - xis - ((lm + xis) / ltor) + rm) / ((lm + xis)&
                    & * sla * dwood / klatosa))**a2
        
    
    end function root_bisec 

    function alt_height (cw, dwood) result (height)
        use types
        use constants

        real(r_8), intent(in) :: cw !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: dwood !trait
        real(r_8) :: diameter !variable internal
        
        real(r_8) :: height !output

        !Calculo diameter (necessary to height)
        diameter = ((4*cw)/dwood*pi*k_allom2)**(1/(2+k_allom3))

        !Height 
        height = k_allom2*(diameter**k_allom3)

    end function alt_height


    !SUBROTINA PARA ALOCAÇÃO DE CARBONO (baseado no LPJmL - Fire)
    subroutine alloc (lm, rm, wm, bminc, dwood, sla, lm2, rm2, cw2)

        use types
        use constants 

        !INPUTS_VARIABLES - These variables comes from module_costants.f90
        real(r_8), intent(in) :: lm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: rm !in kgC/m2 but use in gC/ind [transformation below]
        real(r_8), intent(in) :: wm !wood mass (sap and heart juntos)
        real(r_8), intent(in) :: bminc
        real(r_8), intent(in) :: dwood, sla !traits

        ! real(r_8), intent(in) :: sm !in kgC/m2 but use in gC/ind [transformation below]
        ! real(r_8), intent(in) :: hm !in kgC/m2 but use in gC/ind [transformation below]
        ! real(r_8), intent(in) :: nind 
        ! real(r_8), intent(in) :: height


        !OUTPUT_VARIABLES 
        real(r_8), intent(out) :: lm2 !leaf carbon pool (in KgC/m2) after allocation
        real(r_8), intent(out) :: rm2 !root carbon pool
        real(r_8), intent(out) :: cw2 !wood carbon pool
        real(r_8) :: hm2
        real(r_8) :: sm2

        !INTERNAL_VARIABLES
        real(r_8) :: lm_inc !leaf increment
        real(r_8) :: rm_inc !root increment
        real(r_8) :: sm_inc, hm_inc
        real(r_8) :: bminc_ind !biomass increment individual
        real(r_8) :: leaf_ind, root_ind, sap_ind, heart_ind, wood_ind !carbon pools individual (gC/ind)
        real(r_8) :: height
        real(r_8) :: woodens !dwood transformado
        real(r_8) :: lm1 !FUNCTION
        real(r_8) :: lminc_min !FUNCTION
        real(r_8) :: rminc_min !FUNCTION
        real(r_8) :: nind = 3

        !VARIABLES TO BISECTION LOGIC
        logical   :: normal
        integer(i_4) :: i 
        real(r_8) :: x1             !working vars in bisection
        real(r_8) :: x2
        real(r_8) :: rtbis
        real(r_8) :: dx
        real(r_8) :: xmid
        real(r_8) :: fmid 
        real(r_8) :: fx1
        real(r_8) :: sign

        ! ========================!
        ! ========= START ========!
        ! ========================!

        ! 1) Transform bminc (biomass increment) to individual & gC
        bminc_ind = (bminc/nind)*1.D3

        ! 2) Transform carbon pools in gC/ind. 
        leaf_ind = (lm/nind)*1.D3
        root_ind = (rm/nind)*1.D3
        sap_ind = ((wm*0.05)/nind)*1.D3
        heart_ind = ((wm*0.95)/nind)*1.D3
        wood_ind = (wm/nind)*1.D3

        ! 3) Transform dwood 
        woodens = (dwood*1000000)

        ! call functions to logic 
        height = alt_height(wood_ind, woodens)
        lm1 = leaf_mass(sap_ind,height,sla,woodens)
        lminc_min = leaf_min(lm1,leaf_ind)
        rminc_min = root_min(lm1,root_ind)

        if (rminc_min .gt. 0.0D0 .and. lminc_min .gt. 0.0D0 .and. &
            & (rminc_min + lminc_min) .lt. bminc_ind) then 

            !Normal allocation (positive increment to all living C compartments)
            normal = .true. 

            print*, 'NORMAL ALLOCATION'

            !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
            !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
            !to its allometric minimum, because it should be able to be fulfilled, i.e.:

            !Start to find root procedure (relate to bisection method)

            !intervalo onde está minha raiz (o resultado da minha bisection)

            x1 = lminc_min
            x2 = (bminc_ind - (leaf_ind / ltor - root_ind)) / (1. + 1. / ltor)
            
            dx = x2 - x1

            if (dx < 0.01) then !0.01 é a precisão da bisection. 

                !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
                !assume that the leafmass increment is equal to the midpoint between the values and skip 
                !the root finding procedure

                lm_inc = x1 + 0.5 * dx

            else
                !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al., 1986, p 346)
                !There should be exactly one solution (no proof presented, but Steve has managed one).
                    
                dx = dx/nseg

                fx1 = root_bisec(leaf_ind,sap_ind,heart_ind,&
                &root_ind,bminc_ind,sla,woodens,x1)

                !Find approximate location of leftmost root on the interval (x1,x2).
                !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

                fmid = fx1
                xmid = x1

                i = 1

                do
                    xmid = xmid + dx

                    fmid = root_bisec(leaf_ind,sap_ind,heart_ind,&
                    &root_ind,bminc_ind,sla,woodens,xmid)

                    if (fmid * fx1 .le. 0. .or. xmid .ge. x2) exit  !sign has changed or we are over the upper bound

                    i = i + 1
                enddo

                !the interval that brackets zero in f(x) becomes the new bounds for the root search

                x1 = xmid - dx
                x2 = xmid

                !Apply bisection method to find root on the new interval (x1,x2)
                fx1 = root_bisec(leaf_ind,sap_ind,heart_ind,&
                &root_ind,bminc_ind,sla,woodens,x1)

                if (fx1 .gt. 0.) then
                    sign = -1.
                else
                    sign =  1.
                end if

                rtbis = x1
                dx = x2 - x1

                !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
                !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment

                i = 1

                do 
                    dx   = 0.5 * dx
                    xmid = rtbis + dx

                    !calculate fmid = f(xmid) [eqn (22)]

                    fmid = root_bisec(leaf_ind,sap_ind,heart_ind,&
                    &root_ind,bminc_ind,sla,woodens,xmid)

                    if (fmid * sign .le. 0.) then 
                        rtbis = xmid
                    endif

                    if (dx .lt. xacc .or. abs(fmid) .le. yacc) exit

                    i = i + 1
                end do

                !Now rtbis contains numerical solution for lminc_ind given eqn (22)

                    lm_inc = rtbis
            
            end if  !x2-x1 block

            !Calculate increments in other compartments using allometry relationships
            
            rm_inc = (lm + lm_inc) / ltor - root_ind       !eqn (9)
            
            sm_inc = bminc_ind - lm_inc - rm_inc !eqn (1)

        else 

            !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
            
            normal = .false.

            print*, 'ANNORMAL ALLOCATION'

            !Attempt to distribute this year's production among leaves and roots only

            lm_inc = (bminc_ind-leaf_ind/ltor+root_ind)/(1.+1./ltor)  !eqn (33)

            if (lm_inc .gt. 0.) then

                !Positive allocation to leafmass

                rm_inc = bminc_ind - lm_inc  !eqn (31)
                
                !Add killed roots (if any) to below-ground litter

                if (rm_inc .lt. 0.) then

                    lm_inc = bminc_ind
                    rm_inc = (leaf_ind + lm_inc) / ltor - root_ind

                    !litter_bg(pls) = litter_bg(pls) + abs(rminc_ind(pls)) * nind(pls)

                end if
                
                i = 1

            else

                !Negative allocation to leaf mass

                rm_inc = bminc_ind
                lm_inc = (root_ind + rm_inc) * ltor - leaf_ind  !from eqn (9)

                !Add killed leaves to litter

                !litter_ag_fast(pls) = litter_ag_fast(pls) + abs(lminc_ind(pls)) * nind(pls)
                
                i = 2

            endif

            !Calculate sminc_ind (must be negative)
      
            sm_inc = (leaf_ind + lm_inc) * sla /&
            & klatosa * woodens * height - sap_ind  !eqn (35)

            !Convert killed sapwood to heartwood

            hm_inc = heart_ind + abs(sm_inc)

        endif

        !Increment C compartments - OUTPUT FINAL (kgC/m²)

        lm2 = ((lm2 + lm_inc)*nind)/1.D3
        rm2 = ((rm2 + rm_inc)*nind)/1.D3
        sm2 = ((sm2 + sm_inc)*nind)/1.D3
        hm2 = ((hm2 + hm_inc)*nind)/1.D3
        cw2 = (sm2 + hm2)


    end subroutine alloc
    
end module allocation_reform