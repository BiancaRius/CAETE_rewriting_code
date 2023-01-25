 module allocation
    use types
    use allometric_params

implicit none
 
private

public :: allocate

contains

 ! program allocation !module to test allocation logic module of LPJmL-Fire
    subroutine allocate(lm, hm, sm, rm, dwood, sla, nind, bminc_ind, height,&
        & cl_inc, cw_inc, ch_inc, cs_inc, cr_inc, ctotal_inc,&
        & lm_2, ch_2, cs_2, rm_2, cw_2) 


        !VARIABLES [INPUT] - vem do módulo de estabelecimento/FPC
        real(r_8), intent(in) :: lm, rm, hm, sm, height !valor atualizado do carbono do ano passado.
        real(r_8), intent(in) :: sla, nind, dwood !variable traits
        real(r_8), intent(in) :: bminc_ind !valor aleatório usado tbm no estabelecimento efbmin FPC posterior deve ser a NPP acumulada do ano
        

        !VARIABLES [OUTPUT] - o incremento de cada compartimento de carbono neste ano
        real(r_8), intent(out) :: cl_inc !leaf increment (gC)
        real(r_8), intent(out) :: cw_inc !wood increment (gC)
        real(r_8), intent(out) :: cr_inc !root increment (gC)
        real(r_8), intent(out) :: cs_inc !sapwood increment (gC)
        real(r_8), intent(out) :: ch_inc !heartwood increment (gC)
        real(r_8), intent(out) :: ctotal_inc !total carbon increment (gC)
        real(r_8), intent(out) :: lm_2, ch_2, cs_2, rm_2, cw_2
        ! integer(i_4), parameter :: stdout = output_unit

        real(r_8) :: pi4
        real(r_8) :: a1
        real(r_8) :: a2
        real(r_8) :: a3
        real(r_8) :: hm_new = 0.0D0

        
        real(r_8) :: dwood_gm3
        real(r_8)  :: litter_ag_fast
        real(r_8)  :: litter_bg



        !Local Variables
        real(r_8)  :: lminc_ind      !individual leafmass increment this year
        real(r_8)  :: rminc_ind      !individual fineroot mass increment this year
        real(r_8)  :: lminc_ind_min  !min leafmass increment to maintain current sapwood
        real(r_8)  :: rminc_ind_min  !min rootmass increment to support new leafmass
        real(r_8)  :: sminc_ind      !individual sapmass increment this year


        real(r_8)  :: x1             !working vars in bisection
        real(r_8)  :: x2
        real(r_8)  :: rtbis
        real(r_8)  :: dx
        real(r_8)  :: xmid
        real(r_8)  :: root1, root2, root3
        real(r_8)  :: sign
        logical  :: normal

        real(r_8)  :: fx1
        real(r_8)  :: fmid


        real(r_8)  :: leaf_mass_requirement   !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)

        integer(i_4) :: i

        !initializing variables
        cl_inc = 0.
        cw_inc = 0. !wood increment (gC)
        cr_inc = 0. !root increment (gC)
        cs_inc = 0.!sapwood increment (gC)
        ch_inc = 0.!heartwood increment (gC)
        ctotal_inc = 0.!total carbon increment (gC)
        lm_2 = 0.
        ch_2 = 0.
        cs_2 = 0.
        rm_2 = 0.
        cw_2 = 0.

        dwood_gm3 = dwood * 1000000 !*1e6 transforms dwood to gC/m3

            ! ====== TREE ALLOCATION ======
        !allometric leaf mass requirement
        leaf_mass_requirement  = (latosa * sm) / ( dwood * height * sla) 
            
        lminc_ind_min  = leaf_mass_requirement - lm  !eqn (27)

        !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
        !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

        rminc_ind_min  = leaf_mass_requirement  / ltor - rm       !eqn (30)


        if (rminc_ind_min  .gt. 0. .and. lminc_ind_min  .gt. 0. .and. &
            &(rminc_ind_min  + lminc_ind_min ) .le. bminc_ind ) then

            !Normal allocation (positive increment to all living C compartments)
            normal = .true.

            !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
            !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
            !to its allometric minimum, because it should be able to be fulfilled, i.e.:

            !Start to find root procedure (relate to bisection method)

            x1  = lminc_ind_min 
            x2  = (bminc_ind  - (lm  / ltor - rm )) / (1. + 1. / ltor)
            
            dx  = x2  - x1


            if (dx  < 0.01) then

                !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
                !assume that the leafmass increment is equal to the midpoint between the values and skip 
                !the root finding procedure

                lminc_ind  = x1  + 0.5 * dx 

            else
                !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
                !There should be exactly one solution (no proof presented, but Steve has managed one).
                    
                dx  = dx /nseg

                !! ===== FIND ROOT FUNCTION ===== [**must be a function**]

                pi4 = pi/4
                a1 = 2./allom3
                a2 = 1. + a1
                a3 = allom2**a1


                root1  = a3*((sm +bminc_ind -x1 -((lm +x1 )/ltor)+&
                        &rm +hm )/dwood )/pi4-((sm +bminc_ind -x1 -&
                        &((lm +x1 )/ltor)+rm )/((lm +x1 )*sla *&
                        &(dwood )/latosa))**a2

                ! ======================================================

                !evaluate f(x1) = LHS of eqn (22) at x1

                fx1  = root1 

                !Find approximate location of leftmost root on the interval (x1,x2).
                !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).

                fmid  = fx1 
                xmid  = x1 

                i = 1
                do
                    xmid  = xmid  + dx 

                    root2  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                    &rm +hm )/dwood)/pi4-((sm +bminc_ind -xmid -&
                    &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                    &(dwood)/latosa))**a2

                    fmid  = root2 

                    if ((fmid *fx1 ) .le. 0. .or. xmid  .ge. x2 ) exit  !sign has changed or we are over the upper bound

                    if (i > 50) stop 'Too many iterations allocmod'

                    i = i + 1

                enddo
                print*, "Number of loops: ", i," alloc_typeNormal: ", normal

                !the interval that brackets zero in f(x) becomes the new bounds for the root search
                x1  = xmid  - dx 
                x2  = xmid 

                !Apply bisection method to find root on the new interval (x1,x2)
                fx1  = root1 

                if (fx1  .ge. 0.) then
                    sign  = -1.
                else
                    sign  =  1.
                end if

                rtbis  = x1 
                dx     = x2  - x1 

                !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
                !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment

                i = 1
                do
                    dx    = 0.5 * dx 
                    xmid  = rtbis  + dx 


                    root3  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                    &rm +hm )/dwood)/pi4-((sm +bminc_ind -xmid -&
                    &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                    &dwood/latosa))**a2

                    fmid  = root3 

                    if (fmid  * sign  .le. 0.) rtbis  = xmid 

                    if (dx  < xacc .or. abs(fmid ) <= yacc) exit

                    if (i > 50) stop 'Too many iterations allocmod'

                    i = i + 1
                enddo
                print*, "Number of loops: ", i," alloc_typeNormal: ", normal

                !Now rtbis contains numerical solution for lminc_ind given eqn (22)
                lminc_ind  = rtbis
            endif

            !Calculate increments in other compartments using allometry relationships
            rminc_ind  = (lm  + lminc_ind ) / ltor - rm        !eqn (9)

            sminc_ind  = bminc_ind  - rminc_ind  - lminc_ind   !eqn (1)

            ! print*, 'LEAF_INC (gC/ind)', (lminc_ind /1.D3), 'ROOT_INC (gC/ind)', (rminc_ind /1.D3),&
            ! & 'SAP_INC(gC/ind)', (sminc_ind /1.D3), pls, 'NORMAL'

        else

            !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
            !print*, 'anormal'
            normal = .false.

            !Attempt to distribute this year's production among leaves and roots only

            lminc_ind  = (bminc_ind -lm /ltor+rm )/(1.+1./ltor)  !eqn (33)
            !print*,'anormal', lminc_ind

            if (lminc_ind  > 0.) then
            

                !Positive allocation to leafmass

                rminc_ind  = bminc_ind  - lminc_ind   !eqn (31)
                
                !Add killed roots (if any) to below-ground litter

                if (rminc_ind  < 0.) then

                    lminc_ind  = bminc_ind 
                    rminc_ind  = (lm  + lminc_ind ) / ltor - rm 

                    litter_bg  = litter_bg  + abs(rminc_ind ) * nind 

                end if
                
                i = 1

            else

                !Negative allocation to leaf mass

                rminc_ind  = bminc_ind 
                lminc_ind  = (rm  + rminc_ind ) * ltor - lm   !from eqn (9)
            !  print*, 'lminc', lminc_ind
                !Add killed leaves to litter

                litter_ag_fast  = litter_ag_fast  + abs(lminc_ind ) * nind 
                
                i = 2

            endif

            !Calculate sminc_ind (may be negative)
    
            sminc_ind  = (lm  + lminc_ind ) * sla  /&
            & latosa * dwood * height  - sm   !eqn (35)

            !Convert killed sapwood to heartwood

            hm_new  = hm  + abs(sminc_ind )

        !  print*, pls, 'ANNORMAL'


        endif !normal/abnormal allocation
        print*, "I value final: ", i," alloc_typeNormal: ", normal

        !Increment on C compartments 
        cl_inc = lminc_ind
        !print*, 'clinc', cl_inc !GOTOLITTER???
        cr_inc = rminc_ind
        !print*, 'cRinc', cr_inc
        cs_inc = sminc_ind
        !print*, 'csinc', cs_inc
        ch_inc = hm_new
        !print*, 'chinc', hm
        cw_inc = cs_inc + ch_inc

        
        !  print*, 'ctotal_inc', ctotal_inc
        if(cl_inc.le.0.)then
        cl_inc = 0. !this leaf goes to litter
        endif

        if(cs_inc.le.0.) then
        cs_inc = 0. !this sap goes to heart (I think it already done in the code)
        endif
        ctotal_inc = cl_inc + cr_inc + cs_inc + ch_inc

        !Compartments plus increment
        
        lm_2 = lm  + cl_inc

        ch_2 = hm_new + ch_inc

        cs_2 = sm + cs_inc

        rm_2 = rm + cr_inc

        cw_2 = cs_2 + ch_2
    end subroutine allocate

end module allocation
