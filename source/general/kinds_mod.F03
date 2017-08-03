module kinds
    implicit none

    integer, parameter :: cdp=kind((0.d0,0.d0))
    integer, parameter :: dp=kind(0.d0)                   ! double precision
    integer, parameter :: sp=kind(0.)                     ! single precision

    contains

    subroutine kinds_list(unit)
    integer, intent(in) ::  unit

    write(unit,*) '      KINDS:'
    write(unit,*) 'kind         size'
    write(unit,*) 'dp',dp
    write(unit,*) 'sp',sp

    end subroutine
end module kinds
