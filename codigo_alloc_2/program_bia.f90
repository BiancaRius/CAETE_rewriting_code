program program2
    use constants
    use types
    use allocation
    
    implicit none
    integer(i_4) :: i
    real(r_8) :: leaf_in  = 2. !kgC/m2 initial C leaf input 
    real(r_8) :: root_in  = 2.7 !kgC/m2 initial C root input
    real(r_8) :: wood_in  = 10.  !kgC/m2 initial C wood input
    real(r_8) :: sap_in
    real(r_8) :: heart_in 
    real(r_8) :: bminc_in = 3.4 !carbon (NPP) available to be allocated
                                !basically NPPt - NPPt-1. NPP accumulated in the year/month/day
                                !gc/ind/time_step

    real(r_8) :: dens_in = 3 !ind/m2 initial density of individuals 

    real(r_8) :: leaf_out
    real(r_8) :: root_out
    real(r_8) :: sap_out
    real(r_8) :: heart_out



    sap_in   = 87.!0.05*wood_in
    heart_in = 324.!0.95*wood_in

    
        call alloc(leaf_in, root_in, sap_in, heart_in, bminc_in,dens_in,&
            leaf_out, root_out, sap_out, heart_out)

          print*, 2.e-4 * exp(6.15 - 0.46 * log(4. * 12.))

        !update variables
        ! leaf_in = leaf_out
        
       

end program program2