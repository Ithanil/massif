module mcinterfaces
    use kinds
    use mcwalkers

    implicit none

    type, abstract  ::  mcprobdist

    contains
        procedure(calc_prob_int),deferred,pass (this)  ::  calc_prob
    end type

    abstract interface
        subroutine calc_prob_int(this,walkobj,value)
            import dp,walkers,mcprobdist
            class(mcprobdist)   ::  this
            class(walkers)      ::  walkobj
            real(dp)            ::  value
        end subroutine
    end interface

    type, abstract  ::  mcobserv

        integer         ::  nodim,ncsteps,laststep
        real(dp),allocatable    ::  allsum(:),allvals(:,:),lastval(:)
        logical         ::  storeall

    contains
        procedure(init_observ_int),deferred,pass (this)  ::  init_observ
        procedure(calc_observ_int),deferred,pass (this)  ::  calc_observ
        procedure(final_observ_int),deferred,pass (this)  :: final_observ
        procedure   ::  allocate => sample_allocate
        procedure   ::  step => sample_step
        procedure   ::  average => sample_average
        procedure   ::  deallocate => sample_deallocate
    end type

    abstract interface
        subroutine init_observ_int(this)
            import mcobserv
            class(mcobserv) ::   this
        end subroutine
        subroutine calc_observ_int(this,walkobj,value,nodim)
            import dp,walkers,mcobserv
            class(mcobserv) ::   this
            class(walkers)  ::   walkobj
            integer         ::   nodim
            real(dp)        ::   value(nodim)
        end subroutine
        subroutine final_observ_int(this,value,nodim)
            import dp,mcobserv
            class(mcobserv) ::   this
            integer         ::   nodim
            real(dp)        ::   value(nodim)
        end subroutine
    end interface

    type            ::  mcobserv_listel
        class(mcobserv),allocatable ::  element
    end type

contains

    subroutine sample_allocate(this,ncsteps,storeall)
        class(mcobserv) ::  this
        integer     ::  ncsteps
        logical     ::  storeall

        this%ncsteps=ncsteps
        this%storeall=storeall
        allocate(this%allsum(this%nodim),this%lastval(this%nodim))
        if (this%storeall) then
            allocate(this%allvals(this%nodim,this%ncsteps))
            this%allvals(:,:)=0.0_dp
        end if
        this%allsum(:)=0.0_dp
        this%lastval(:)=0.0_dp
        this%laststep=0
    end subroutine

    subroutine sample_step(this,walkobj,value,nodim,calc)

        class(mcobserv) ::   this
        class(walkers)  ::   walkobj
        integer         ::   nodim
        real(dp)        ::   value(nodim)
        logical         ::   calc

        if (calc) then  !if not, then value must be passed
            call this%calc_observ(walkobj,value,nodim)
        end if

        this%laststep=this%laststep+1

        this%allsum(:)=this%allsum(:)+value(:)
        if (this%storeall) then
            this%allvals(:,this%laststep)=value(:)
        end if
        this%lastval(:)=value(:)
    end subroutine

    subroutine sample_average(this)
        class(mcobserv) ::   this

        call this%final_observ(this%allsum(:),this%nodim)

        if (this%laststep.ne.this%ncsteps) then
            write(6,*)  'WARN:  [mcinterfaces]  sample_average:  Premature usage of average function,'// &
                'laststep was',this%laststep,'but ncsteps was set to:',this%ncsteps,' Dividing by laststep.'
        end if

        this%allsum(:)=this%allsum(:)/this%laststep

    end subroutine

    subroutine sample_deallocate(this)
        class(mcobserv) ::  this

        if (this%laststep.gt.0 .and. this%laststep.lt.this%ncsteps) then
            write(6,*)  'WARN:  [mcinterfaces]  sample_deallocate:  Premature usage of deallocate function,'// &
                'laststep was',this%laststep,'but ncsteps was set to:',this%ncsteps,' Will deallocate anyway.'
        end if

        deallocate(this%allsum)
        deallocate(this%lastval)
        if (this%storeall) then
            deallocate(this%allvals)
        end if

    end subroutine
end module mcinterfaces
