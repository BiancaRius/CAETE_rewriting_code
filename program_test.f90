program test2
    use establishment
    use types
    use FPC
    
    implicit none
    

    real(r_8) :: est_pls2
    real(r_8) :: FPC_total_gc2
    real(r_8) :: FPC_ind2
    real(r_8) :: FPC_pls2
    real(r_8) :: npls_alive2 = 1.0D0

    real(r_8) :: csap_new2
    real(r_8) :: cleaf_new2
    real(r_8) :: cheart_new2
    real(r_8) :: croot_new2
    real(r_8) :: cwood_new2
    real(r_8) :: dens_new2

    real(r_8) :: csap_sapl2   
    real(r_8) :: cleaf_sapl2
    real(r_8) :: cheart_sapl2
    real(r_8) :: croot_sapl2

    real(r_8) :: csap_old2 =1.0D0
    real(r_8) :: cleaf_old2=1.0D0
    real(r_8) :: cheart_old2=1.0D0
    real(r_8) :: croot_old2=1.0D0
    real(r_8) :: dens_old2=1.0D0

    real(r_8) :: x1_1 =1.0D0
    
   

        call gc_occupation(x1_1,FPC_total_gc2, FPC_ind2, FPC_pls2)
    
        call establish(npls_alive2,FPC_total_gc2,est_pls2)

        ! print*, est_pls2

        call sapling_allometry(cleaf_sapl2, csap_sapl2, cheart_sapl2, croot_sapl2)

        ! print*, cleaf_sapl2

        call shrink(csap_old2,cleaf_old2,cheart_old2, croot_old2,&
                    est_pls2, dens_old2, csap_sapl2, cleaf_sapl2, cheart_sapl2, croot_sapl2,&
                    csap_new2, cleaf_new2,cheart_new2, croot_new2, cwood_new2,dens_new2)

        

        ! print*, csap_new2,cleaf_new2,cheart_new2, croot_new2, cwood_new2

        !update variables dfh

        

end program test2