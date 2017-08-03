module perbounds
    use kinds
    implicit none

    !An object of the perbound_box_* type can be instantiated to use the perbound procedures
    !without passing the (most probably constant) box and perbound vectors at every
    !procedure call. This is potentially faster in both usage and execution!
    type perbound_box_all

        real(dp), allocatable   ::   boxl(:)
        integer                 ::   ndim

    contains

        procedure   ::  init_obj_all
        procedure   ::  dist_calc
        procedure   ::  dvecdist_calc
        procedure   ::  dvec_calc
        procedure   ::  pershift => pershift_obj_all
        procedure   ::  minimg_dvec => minimg_dvec_obj_all
        procedure   ::  minimg_dvecdist => minimg_dvecdist_obj_all
    end type

    type, extends(perbound_box_all) ::  perbound_box_gen
        logical,  allocatable   ::   perbound(:)
    contains

        procedure   ::  init_obj_gen
        procedure   ::  pershift => pershift_obj_gen
        procedure   ::  minimg_dvec => minimg_dvec_obj_gen
        procedure   ::  minimg_dvecdist => minimg_dvecdist_obj_gen
    end type
contains


    !perbound_box_* type general creator routine
    subroutine perbound_box_create(this,boxl,perbound,ndim)
        class(perbound_box_all),allocatable     ::   this
        real(dp),dimension(ndim)    ::   boxl
        logical, dimension(ndim)    ::   perbound
        integer ::  it,ndim

        write(6,*) 'perbound_box_create...'

        do it=1,ndim
            if (it.eq.ndim) then
                allocate(perbound_box_all::this)
                call this%init_obj_all(boxl,ndim)
                exit
            else
                if ( perbound(it).neqv.perbound(it+1) ) then
                    allocate(perbound_box_gen::this)
                    select type (this)
                        type is (perbound_box_gen)
                        call this%init_obj_gen(boxl,perbound,ndim)
                    end select
                    exit
                end if
            end if
        enddo

    end subroutine

    !--- Type-bound periodic boundary routines for type perbound_box_all ---!

    subroutine init_obj_all(this,boxl,ndim)
        class(perbound_box_all)  ::   this
        real(dp),dimension(ndim) ::   boxl
        integer :: ndim

        allocate(this%boxl(ndim))
        this%boxl(:)=boxl(:)
        this%ndim=ndim

    end subroutine

    subroutine pershift_obj_all(this,xyzvec,ndim)
        class(perbound_box_all) ::  this
        real(dp), dimension(ndim)  ::   xyzvec
        integer :: ndim

        if (this%ndim.eq.ndim) then
            xyzvec(:)=xyzvec(:)-this%boxl(:)*floor(xyzvec(:)/this%boxl(:))
        else
            write(6,*) 'ERR :   [perbounds]   pershift_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine

    subroutine minimg_dvec_obj_all(this,distvec,ndim)
        class(perbound_box_all) ::  this
        real(dp), dimension(ndim)  ::   distvec
        integer :: ndim

        if (this%ndim.eq.ndim) then
            distvec(:)=distvec(:)-this%boxl(:)*anint(distvec(:)/this%boxl(:))
        else
            write(6,*) 'ERR :   [perbounds]   minimg_dvec_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine

    subroutine  minimg_dvecdist_obj_all(this,distvec,dist,ndim)
        class(perbound_box_all) ::  this
        real(dp),dimension(ndim)   ::  distvec
        real(dp)                ::  dist
        integer :: ndim

        if (this%ndim.eq.ndim) then
            distvec(:)=distvec(:)-this%boxl(:)*anint(distvec(:)/this%boxl(:))
            dist=sqrt(dot_product(distvec,distvec))
        else
            write(6,*) 'ERR :   [perbounds] minimg_dvecdist_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine

    subroutine  dist_calc(this,R1,R2,dist,ndim)
    !ndim is needed as an explicit parameter, just use this%ndim when calling
        class(perbound_box_all)     ::   this
        integer     ::  ndim
        real(dp)    ::  R1(ndim),R2(ndim),dist

        real(dp)    ::  disth(ndim)
        disth(:)=R2(:)-R1(:)
        call this%minimg_dvecdist(disth,dist,this%ndim)
    end subroutine

    subroutine  dvecdist_calc(this,R1,R2,distvec,dist,ndim)
    !ndim is needed as an explicit parameter, just use this%ndim when calling
        class(perbound_box_all)     ::   this
        integer     ::  ndim
        real(dp)    ::  R1(ndim),R2(ndim),dist,distvec(ndim)

        distvec(:)=R2(:)-R1(:)
        call this%minimg_dvecdist(distvec,dist,this%ndim)
    end subroutine

    subroutine  dvec_calc(this,R1,R2,distvec,ndim)
    !ndim is needed as an explicit parameter, just use this%ndim when calling
        class(perbound_box_all)     ::   this
        integer     ::  ndim
        real(dp)    ::  R1(ndim),R2(ndim),distvec(ndim)

        distvec(:)=R2(:)-R1(:)
        call this%minimg_dvec(distvec,this%ndim)
    end subroutine

    !--- Type-bound periodic boundary routines for type perbound_box_gen ---!

    subroutine init_obj_gen(this,boxl,perbound,ndim)
        class(perbound_box_gen)  ::   this
        real(dp),dimension(ndim) ::   boxl
        logical, dimension(ndim) ::   perbound
        integer :: ndim

        allocate(this%boxl(ndim),this%perbound(ndim))
        this%boxl(:)=boxl(:)
        this%perbound(:)=perbound(:)
        this%ndim=ndim

    end subroutine

    subroutine pershift_obj_gen(this,xyzvec,ndim)
        class(perbound_box_gen) ::  this
        real(dp), dimension(ndim)  ::   xyzvec
        integer :: it,ndim

        if (this%ndim.eq.ndim) then
            do it=1,this%ndim
                if (this%perbound(it)) then
                    xyzvec(it)=xyzvec(it)-this%boxl(it)*floor(xyzvec(it)/this%boxl(it))
                end if
            enddo
        else
            write(6,*) 'ERR :   [perbounds]   pershift_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine

    subroutine minimg_dvec_obj_gen(this,distvec,ndim)
        class(perbound_box_gen) ::  this
        real(dp), dimension(ndim)  ::   distvec
        integer :: it,ndim

        if (this%ndim.eq.ndim) then
            do it=1,this%ndim
                if (this%perbound(it)) then
                    distvec(it)=distvec(it)-this%boxl(it)*anint(distvec(it)/this%boxl(it))
                end if
            enddo
        else
            write(6,*) 'ERR :   [perbounds]   minimg_dvec_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine

    subroutine  minimg_dvecdist_obj_gen(this,distvec,dist,ndim)
        class(perbound_box_gen) ::  this
        real(dp),dimension(ndim)   ::  distvec
        real(dp)                ::  dist
        integer :: it,ndim

        if (this%ndim.eq.ndim) then
            do it=1,this%ndim
                if (this%perbound(it)) then
                    distvec(it)=distvec(it)-this%boxl(it)*anint(distvec(it)/this%boxl(it))
                end if
            enddo
            dist=sqrt(dot_product(distvec,distvec))
        else
            write(6,*) 'ERR :   [perbounds] minimg_dvecdist_obj_all:  Passed vector dimensions do not match! Stop.'
            stop
        end if
    end subroutine


    !-------              Alternative non-type-bound routines                 -----!

    !General periodic boundary routines,
    !called with box vector boxl and logical vector perbound.

    subroutine pershift_gen(xyzvec,boxl,perbound,ndim)
        real(dp), dimension(ndim)  ::   xyzvec,boxl
        logical,  dimension(ndim)  ::   perbound

        integer ::  it,ndim

        do it=1,ndim
            if (perbound(it)) then
                xyzvec(it)=xyzvec(it)-boxl(it)*floor(xyzvec(it)/boxl(it))

            end if
        enddo

    end subroutine

    subroutine minimg_dvec_gen(distvec,boxl,perbound,ndim)
        real(dp), dimension(ndim)  ::   distvec,boxl
        logical,  dimension(ndim)  ::   perbound

        integer ::  it,ndim

        do it=1,ndim
            if (perbound(it)) then
                distvec(it)=distvec(it)-boxl(it)*anint(distvec(it)/boxl(it))
            end if
        enddo

    end subroutine

    subroutine minimg_dvecdist_gen(distvec,dist,boxl,perbound,ndim)
        real(dp), dimension(ndim)  ::    distvec,boxl
        real(dp)                ::    dist
        logical,  dimension(ndim)  ::    perbound

        integer ::  it,ndim

        do it=1,ndim
            if (perbound(it)) then
                distvec(it)=distvec(it)-boxl(it)*anint(distvec(it)/boxl(it))
            end if
        enddo
        dist=sqrt(dot_product(distvec,distvec))

    end subroutine

    ! *_all routines are for the common case, that periodic boundaries are applied
    ! to all directions. Otherwise use the general ones above and pass the logical
    ! vector perbound.

    subroutine pershift_all(xyzvec,boxl,ndim)
        real(dp), dimension(ndim)  ::   xyzvec,boxl
        integer     ::   ndim

        xyzvec(:)=xyzvec(:)-boxl(:)*floor(xyzvec(:)/boxl(:))

    end subroutine

    subroutine minimg_dvec_all(distvec,boxl,ndim)
        real(dp), dimension(ndim)  ::    distvec,boxl
        integer     ::   ndim

        distvec(:)=distvec(:)-boxl(:)*anint(distvec(:)/boxl(:))

    end subroutine

    subroutine minimg_dvecdist_all(distvec,boxl,dist,ndim)
        real(dp), dimension(ndim)  ::    distvec,boxl
        real(dp)    ::    dist
        integer     ::    ndim

        distvec(:)=distvec(:)-boxl(:)*anint(distvec(:)/boxl(:))
        dist=sqrt(dot_product(distvec,distvec))

    end subroutine

end module perbounds
