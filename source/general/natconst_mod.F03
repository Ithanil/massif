module natconst
    use kinds
    implicit none

    real(dp), parameter ::  pi=acos(-1.0_dp)

    real(dp), parameter ::  mtoA = 10.0_dp**(-10)
    real(dp), parameter ::  stons = 10.0_dp**(-9)
    real(dp), parameter ::  stops = 10.0_dp**(-12)
    real(dp), parameter ::  stofs = 10.0_dp**(-15)

    !values from http://www.chemie.fu-berlin.de/chemistry/general/constants_mathmode.html
    real(dp), target    ::  hplank_SI=6.6260755e-34_dp      ![J*s]=[kg*m^2/s^2/K]
    real(dp), target    ::  hbar_SI=1.05457266e-34_dp       ![J*s]=[kg*m^2/s^2/K]
    real(dp), target    ::  kboltz_SI= 1.380658e-23_dp      ![J/K]=[kg*m^2/s^2/K]
    real(dp), target    ::  elcharge_SI=1.60217733e-19_dp   ![C]
    real(dp), target    ::  spolight_SI=2.99792458e8_dp     ![m/s]

    !values from above with conversion J*s -> K, m -> A, s -> ps
    real(dp), target    ::  hplank_VMC=43.761232706_dp      ![ K ]
    real(dp), target    ::  hbar_VMC=6.964816501_dp         ![ K ]
    real(dp), target    ::  kboltz_VMC=1.0_dp               !??
    real(dp), target    ::  elcharge_VMC=1.60217733e-19_dp  ![C]
    real(dp), target    ::  spolight_VMC=2.99792458e6_dp    ![A/ps]

    !pointers
    real(dp), pointer   ::  hplank,hbar,kboltz,elcharge,spolight

contains
    subroutine natconst_init(name)
        character(*)    ::  name

        if (name.eq.'SI') then
            hplank   => hplank_SI
            hbar     => hbar_SI
            kboltz   => kboltz_SI
            elcharge => elcharge_SI
            spolight => spolight_SI
        else if (name.eq.'VMC') then
            hplank   => hplank_VMC
            hbar     => hbar_VMC
            kboltz   => kboltz_VMC
            elcharge => elcharge_VMC
            spolight => spolight_VMC
        end if
    end subroutine
end module natconst
