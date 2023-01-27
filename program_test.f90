program test2
    use establishment
    use types
    use FPC
    
    
    implicit none
    

    real(r_8) :: est_pls2
    real(r_8) :: FPC_total_gc2
    real(r_8) :: FPC_ind2
    real(r_8) :: FPC_pls2
    
    !!!!!!! ALIVE PLSs
    ! real(r_8) :: npls_alive2 = 1.0D0


    !C that enters in the logic (kgC/m2) for each PLS
    real(r_8) :: cleaf_pls_in  = 1.0D0
    real(r_8) :: csap_pls_in   = 1.0D0
    real(r_8) :: cheart_pls_in = 10.0D0
    real(r_8) :: croot_pls_in  = 1.0D0

    !density that enters in the logic (ind/m2) for each PLS
    real(r_8) :: dens_pls_in = 10. 
    
    !C that is an output from the logic (kgC/m2) for each PLS
    real(r_8) :: cleaf_pls_out  
    real(r_8) :: csap_pls_out   
    real(r_8) :: cheart_pls_out 
    real(r_8) :: croot_pls_out  

    !density that is an output from the logic (ind/m2) for each PLS
    real(r_8) :: dens_pls_out
    
   

        call gc_occupation(dens_pls_in,cleaf_pls_in,csap_pls_in, cheart_pls_in,&
            croot_pls_in, FPC_total_gc2, FPC_ind2, FPC_pls2)
        ! print*, FPC_pls2
        ! call establish(npls_alive2,FPC_total_gc2,est_pls2)

        ! ! print*, est_pls2

        ! call sapling_allometry(cleaf_sapl2, csap_sapl2, cheart_sapl2, croot_sapl2)

        ! ! print*, cleaf_sapl2

        ! call shrink(csap_old2,cleaf_old2,cheart_old2, croot_old2,&
        !             est_pls2, dens_old2,&
        !             csap_new2, cleaf_new2,cheart_new2, croot_new2, cwood_new2,dens_new2)

        

        ! print*, csap_new2,cleaf_new2,cheart_new2, croot_new2, cwood_new2

        !update variables dfh

        

end program test2