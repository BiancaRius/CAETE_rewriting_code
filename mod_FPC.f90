!Este módulo calcula as variáveis para indivíduos médios de cada PLS
!Designa-se também a calcular o FPC (foliar projective cover) através
!da estruturação alométrica dos PLSs. Além disso, calcula a mortalidade por espaço,
!ou seja, considerando uma área disponível de 1ha, os PLSs lenhosos só podem ocupar
!95% desta área. Caso ultrapassem esse montante, haverá penalização (redução na densidade 
!de indivíduos) que impactará os compartimentos de carbono. Caso não ultrapassem, será 
!possível que novos indivíduos se estabeleçam. Tais indivíduos são chamados de "sapling",
!apresentam alometria própria, mas não há cohort, de modo que ao incorporá-los na população do 
!PLSs, deve haver um ajustamento no balanço de carbono, o "shrink", que representa uma reestruturação
!dos indivíduos médios.
!A grande maioria dos cálculos aqui é feita para o indivíduo médio, portanto, a maioria das variáveis são
!divididas pela densidade (indivíduos/m²). Pools são usados em gC.
!Este código é baseado principalmente no modelo LPJ (Sitch et al 2003 e Smith et al 2001) e no código
!do LPJ-MLFire

module FPC


    use types
    use params
    use establishment
    

implicit none

private

public :: gc_occupation

contains

    subroutine gc_occupation(dens_pls, FPC_total_gc, FPC_ind, FPC_pls)
        !VARIABLE INPUTS
        real(r_8), intent(in) :: dens_pls
        
        !VARIABLES OUTPUTS 
        real(r_8), intent(out) :: FPC_total_gc !(m2) total FPC in a grid cell considering all PLS
        real(r_8), intent(out) :: FPC_ind !(m2)avg individual FPC
        real(r_8), intent(out) :: FPC_pls !(m2) FPC for each PLS

        !INTERNAL VARIABLES
        !!!Structuring variables
        !!!!!!!!!!!!!!!!Must come from allometry, provisorialy will have a specified value
        real(r_8) :: diam_pls
        real(r_8) :: crown_area_pls
        real(r_8) :: lai_pls
        real(r_8) :: height_pls

        !gc occupation variables
        real(r_8) :: exc_area_gc
        

        !mortality variables
        real(r_8) :: nind_kill_FPC
        real(r_8) :: z2
        
        !initializing variables
        FPC_total_gc = 0.0D0
        FPC_ind = 0.0D0
        FPC_pls = 0.0D0
        diam_pls = 0.0D0
        crown_area_pls = 0.0D0
        lai_pls = 0.0D0
        height_pls = 0.0D0
        nind_kill_FPC = 0.0D0

    !Structuring PLSs [diameter, crown area, leaf area index and height]
        diam_pls = 100.
        crown_area_pls = 100. 
        lai_pls = 100.
        height_pls = 100.

    !dens is an input or output?
            

    !------------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Calculating Foliage Projective Cover of average individual(FPC_ind) and Fractional
    !Projective cover for PLS (FPC_pls2)
        FPC_ind = (1 - (exp (-0.5 * lai_pls)))
    
        FPC_pls = (crown_area_pls * dens_pls) * FPC_ind 
        
!!!!!!!!!!!!!!!!!!calculo FPC_total_gc é feito a partir do acúmulo de todos os PLS.
        !!!!!!!!!!por enquanto trabalhando com 1

        FPC_total_gc = 1000.

    !Verify if total FPC is gt or lt 95% of grid cell area
        if (FPC_total_gc.gt.gc_area_95) then
            print*, 'gt'

            call exc_area(FPC_total_gc, FPC_pls, exc_area_gc, nind_kill_FPC)
            ! print*, exc_area_gc

        else

            !call establish(npls_alive, FPC_total_gc, est_pls)
            !call subroutine of establishment (only the c final?)
            !call establishment

            !call shrink
            print*, 'lt'
        
        endif

    !calculates mortality due to gc occupation, greff, and wood density
        call mortality(nind_kill_FPC,z2)

    end subroutine gc_occupation

    subroutine exc_area(FPC_total_gc, FPC_pls, exc_area_gc, nind_kill_FPC)
        !VARIABLES INPUTS
        real(r_8), intent(in) :: FPC_total_gc !(m2) total FPC in a grid cell considering all PLSs
        real(r_8), intent(in) :: FPC_pls !(m2) FPC for each PLS

        !VARIABLES OUTPUTS
        real(r_8), intent(out) :: exc_area_gc !(m2) excendent in area for a gc
        real(r_8), intent(out) :: nind_kill_FPC !number of avg ind. that will die due to (ind/m2) excedent area occupation

        !INTERNAL VARIABLES
        real(r_8) :: FPC_dec_pls
        
        !initializing variables
        exc_area_gc = 0.0D0
        FPC_dec_pls = 0.0D0
        nind_kill_FPC = 0.0D0

        ! Excedent area           
        exc_area_gc = FPC_total_gc - gc_area_95

        !FPC decay to fit the 95% of gc occupation
!!!!!!!! if(FPC_pls_2(j,k).gt.0) then !se a ocupação é maior q zero = PLS vivo.
        !     FPC_dec(j,k) = min(FPC_pls_2(j,k), exc_area(k)*(FPC_pls_2(j,k)/FPC_total_accu_2(k)))
        ! else
        !     FPC_dec(j,k) = 0
        ! endif

        FPC_dec_pls = min(FPC_pls, exc_area_gc * (FPC_pls/FPC_total_gc))
        
        !kill ind due to exc are
        nind_kill_FPC = 10.*3.!(dens1(j,k) * FPC_dec(j,k))/FPC_pls_2(j,k) !NIND_KILL.
        ! print*, nind_kill_FPC

    end subroutine exc_area

    subroutine mortality(nind_kill_FPC, z2)
        !calculates mortality (due to gc occupation, greff, and wood density) and
        !carbon loss due to residence time

        !VARIABLES INPUTS
        real(r_8), intent(in) :: nind_kill_FPC

        !VARIABLES OUTPUTS
        real(r_8), intent(out) :: z2

        !INTERNAL VARIABLES
        ! real(r_8) :: nind_kill_FPC

        !initializing variables
        z2 = 0.0D0

        z2 = 2.**1.

        print*, nind_kill_FPC

        ! greff_pls = carbon_increment(j,k)/(cl2(j,k)*spec_leaf(j,k)) !growth efficiency in m2/gC

        ! mort_wd = exp(-2.66+(0.255/dwood(j,k))) !mort by wd (from Sakschewski et al 2015)

        ! mort_greff = mort_wd/(1+k_mort2*greff_pls)
                    

        !             nind_kill_greff(j,k) = (dens1(j,k) * mort_greff(j,k)) !valores absolutos de ind.
        !             ! print*, 'NIND_KILL', nind_kill_greff(j,k)

        !             !soma nind_kill
        !             nind_kill_total(j,k) = nind_kill_FPC(j,k) + nind_kill_greff(j,k) !em ind/m2

        !             !mort(j,k) = ((dens1(j,k)-nind_kill_total(j,k))/dens1(j,k)) !quanto vai morrer em relação a densidade atual
        !             mort(j,k) = nind_kill_total(j,k)/dens1(j,k) !% de ind que devem morrer em relação ao que tinha anteriormente


    end subroutine mortality
    
end module FPC



