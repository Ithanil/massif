module mcwalkers
    use kinds
    use perbounds

    implicit none

    type moveconfig
        real(dp) dummy
    end type
    type physconfig
        real(dp) dummy
    end type
    type, extends(physconfig)   ::   wmass_config
        real(dp), allocatable   ::   mass(:)
    end type
    type, extends(wmass_config)   ::   wcharg_config
        real(dp), allocatable   ::   charge(:)
    end type

    type walkers
        !Extensible walkers type, suitable for Monte Carlo
        !Create with "call walkers_create(W,...)", where W is class(walkers),allocatable.
        !Rerun init routines with W%init(...) or reinit coords with W%init_coords().
        !Move walkers with W%walk().

        integer :: nwalkers,ndim
        real(dp):: nrho
        real(dp), allocatable   ::  R(:,:),boxl(:)
        procedure(init_coords_std), pointer  :: init_coords_ptr => null()
        class(moveconfig), allocatable       :: mconf_dummy
        class(physconfig), allocatable       :: pconf   !can store masses, charges etc.
        class(perbound_box_all), allocatable :: perbounds
    contains

        procedure :: init_coords => init_coords_handle
        procedure :: print => walkers_printout
        !to be replaced by childs:
        procedure :: init => init_move_null
        procedure :: walk => move_null

    end type

    type, extends(moveconfig) :: config_all_primitive
        real(dp)    ::  maxdr
    end type
    type, extends(walkers) :: walkers_all_primitive
        class(config_all_primitive), allocatable :: mconf
    contains
        procedure :: init => init_all_primitive
        procedure :: walk => move_all_primitive
    end type


contains
    subroutine init_coords_handle(this)
        class(walkers) :: this
        call this%init_coords_ptr()
    end subroutine

    !---               General init routines                 ---!

    subroutine walkers_create(this,nwalkers,ndim,boxl,perbound,mconfig,pconfig,init_coords_name)
        !General constructor based on standard parameters and extension parameters in mconfig.
        !The dynamic type of mconfig selects the move routine to be used.
        class(walkers), allocatable  ::  this
        class(moveconfig)   ::  mconfig
        class(physconfig)   ::  pconfig
        integer             ::  nwalkers,ndim
        real(dp)            ::  boxl(ndim)
        character(*)        ::  init_coords_name
        logical             ::  perbound(ndim)

        write(6,*) 'walkers_create...'

        select type (mconfig)

            type is (config_all_primitive)
            allocate(walkers_all_primitive::this)
            select type (this)
                type is (walkers_all_primitive)
                allocate(config_all_primitive::this%mconf)
            end select

            type is (moveconfig)
            allocate(walkers::this)
            allocate(moveconfig::this%mconf_dummy)
            write(6,*)  'WARN: [mcwalkers] walkers_create:  Using moveconfig dummy move config,'// &
                ' not suitable for real calculations!'
        end select

        select type (pconfig)

            type is (physconfig)
            allocate(physconfig::this%pconf)
            select type (pconf_obj => this%pconf)
                type is (physconfig)
                pconf_obj=pconfig
            end select
            write(6,*) 'WARN: [mcwalkers]  walkers_create:  Using physconfig dummy physical config,'// &
                ' not suitable for real calculations!'

            type is (wmass_config)
            allocate(wmass_config::this%pconf)
            select type (pconf_obj => this%pconf)
                type is (wmass_config)
                pconf_obj=pconfig
            end select

            type is (wcharg_config)
            allocate(wcharg_config::this%pconf)
            select type (pconf_obj => this%pconf)
                type is (wcharg_config)
                pconf_obj=pconfig
            end select

        end select

        call this%init(nwalkers,ndim,boxl,perbound,mconfig,init_coords_name)
        call walkers_pershift(this)
    end subroutine

    subroutine init_walkers_std(this,nwalkers,ndim,boxl,perbound,init_coords_name)
        !General init routine that can be called within extension init routines.
        class(walkers)  ::  this
        integer         ::  nwalkers,ndim
        real(dp)        ::  boxl(ndim)
        character(*)    ::  init_coords_name
        logical         ::  perbound(ndim)

        integer         ::  it

        this%nwalkers=nwalkers
        this%ndim=ndim
        allocate(this%R(this%ndim,this%nwalkers),this%boxl(this%ndim))
        this%R(:,:)=0.0_dp
        this%boxl(:)=boxl(:)
        this%nrho=this%nwalkers
        do it=1,this%ndim
            this%nrho=this%nrho/boxl(it)
        enddo
        call perbound_box_create(this%perbounds,this%boxl,perbound,this%ndim)

        select case (init_coords_name)

            case ('std')
                this%init_coords_ptr=>init_coords_std
            case default
                this%init_coords_ptr=>init_coords_std
                write(6,*) 'WARN: [walkers] init_walkers_std: No matching init_coords name! '// &
                    'Defaulting to init_coords_std ("std").'

        end select

        call this%init_coords()

    end subroutine

    !----       Coordinate init routines       ----!

    subroutine init_coords_std(this)
        !Calls cube or non-cube standard init, depending on the box.
        class(walkers)  ::  this

        integer :: it

        do it=1,this%ndim
            if (it.eq.this%ndim) then
                call init_cube_std(this)
                exit
            else
                if ( this%boxl(it).ne.this%boxl(it+1) ) then
                    call init_noncube_std(this)
                    exit
                end if
            end if
        enddo

    end subroutine


    subroutine init_cube_std(this)
        !Standard walker initialization in Francescos diploma thesis

        class(walkers)  ::  this

        real(dp),allocatable :: lattice(:,:),lvec(:),drvecs(:,:)
        logical ,allocatable :: assigned(:)

        real(dp)    ::  wldx
        integer     ::  nlat,nlatpdim
        integer     ::  it,it2,itl,ity,itz
        logical     ::  found

        real(dp),parameter  :: disph=0.05_dp


        !calculate smallest integer nlat, such that nlat**3>=nwalkers
        nlat=int(this%nwalkers**(1.0_dp/this%ndim))+1
        nlatpdim=nlat**this%ndim
        if (nlatpdim.lt.this%nwalkers) then
            nlat=nlat+1
            write(6,*) 'NOTE: [walkers] init_cube_std: nlat estimation via nlat=int(nw^(1/dim))+1 did not result'// &
                'in nlat^dim>=nw (numerical error)! Program continues with nlat=nlat+1.'
        end if

        allocate(lattice(this%ndim,nlatpdim),lvec(nlat),assigned(nlatpdim),drvecs(this%ndim,this%nwalkers))

        !construct a cubic lattice within the box for nlatpdim positions

        wldx=this%boxl(1)/nlat

        lvec(:)=0.5_dp*wldx
        do it=2,nlat
            lvec(it)=lvec(it)+(it-1)*wldx
        enddo

        ity=0
        itz=0
        do it=1,nlatpdim,nlat
            lattice(1,it:it+nlat-1)=lvec(:)

            if (this%ndim.gt.1) then
                if (this%ndim.gt.2) then
                    if (mod(it-1,nlat*nlat).eq.0) then
                        ity=0
                        itz=itz+1

                        lattice(3,it:it+nlat*nlat-1)=lvec(itz)

                    end if
                end if
                ity=ity+1

                lattice(2,it:it+nlat-1)=lvec(ity)

            end if
        enddo

        !Choose nwalkers positions randomly:
        !To have an unbiased random selection with simple code, an iterator itl, which
        !represents lattice positions, iterates randomly in the range from 1 to nlatpdim.
        !If it hits an unassigned position, it is assigned to the current it (which takes walker numbers).

        assigned(:)=.false.
        do it=1,this%nwalkers
            found=.false.
            do while(.not.found)
                itl=int(nlatpdim*rand())+1
                if (.not.assigned(itl)) then
                    this%R(:,it)=lattice(:,itl)
                    assigned(itl)=.true.
                    found=.true.
                    exit
                end if
            enddo
        enddo

        !Random displacement by maximal disph/2 of lattice parameter
        call random_number(drvecs(:,:))
        drvecs(:,:)=disph*wldx*(drvecs(:,:)-0.5_dp)
        this%R(:,:)=this%R(:,:)+drvecs(:,:)
    end subroutine

    subroutine init_noncube_std(this)
        class(walkers)  ::  this

    end subroutine


    !---                Mover init and move routines                 ---!

    !periodic boundary shift
    subroutine walkers_pershift(this)
        class(walkers)  :: this
        integer         :: it

        do it=1,this%nwalkers
            call this%perbounds%pershift(this%R(:,it),this%ndim)
        enddo

    end subroutine

    !Routines for parent type (walkers), not for practical usage (null)

    subroutine init_move_null(this,nwalkers,ndim,boxl,perbound,mconfig,init_coords_name)
        class(walkers)     ::  this
        class(moveconfig)  ::  mconfig
        integer         ::  nwalkers,ndim
        real(dp)        ::  boxl(ndim)
        character(*)    ::  init_coords_name
        logical         ::  perbound(ndim)

        call init_walkers_std(this,nwalkers,ndim,boxl,perbound,init_coords_name)
        write(6,*) 'WARN: [walkers] init_move_null: Move routine not set!!'

        return
    end subroutine

    subroutine move_null(this)
        class(walkers)     ::  this

        write(6,*) 'WARN: [walkers] move_null: Move routine not set!! Doing nothing.'
    end subroutine


    !Routines for child type walkers_all_primitive

    subroutine init_all_primitive(this,nwalkers,ndim,boxl,perbound,mconfig,init_coords_name)
        class(walkers_all_primitive)     ::  this
        class(moveconfig)     ::  mconfig
        integer         ::  nwalkers,ndim
        real(dp)        ::  boxl(ndim)
        character(*)    ::  init_coords_name
        logical         ::  perbound(ndim)

        select type (mconfig)
            class is (config_all_primitive)
            this%mconf%maxdr=mconfig%maxdr
        call init_walkers_std(this,nwalkers,ndim,boxl,perbound,init_coords_name)

        class is (moveconfig)
            write(6,*) 'ERR:  [mcwalkers] init_all_primitive: Entered specialized init, but have'// &
                       ' generic config type! Stop.'
        end select
        return
    end subroutine

    subroutine move_all_primitive(this)
        class(walkers_all_primitive)    ::  this
        real(dp),allocatable    ::  drvecs(:,:)

        allocate(drvecs(this%ndim,this%nwalkers))
        call random_number(drvecs(:,:))
        drvecs(:,:)=2.0_dp*this%mconf%maxdr*(drvecs(:,:)-0.5_dp)
        this%R(:,:)=this%R(:,:)+drvecs(:,:)
        call walkers_pershift(this)
    end subroutine

    !--         Printout routine        --!

    subroutine walkers_printout(this)
        class(walkers)  ::  this
        integer         ::  it

        write(6,*)  '----   Walkers Data Printout   ----'
        write(6,*)
        select type (this)
            type is (walkers)
            write(6,*)  'Dynamic type:  walkers'
            type is (walkers_all_primitive)
            write(6,*)  'Dynamic type:  walkers_all_primitive'
        end select

        write(6,*) 'R(ndim,nwalkers)',this%ndim,this%nwalkers
        do it=1,this%nwalkers
            write(6,*) this%R(:,it)
        enddo
        write(6,*) 'Box:',this%boxl(:)
        write(6,*)  'meanxyz:'
        do it=1,this%ndim
            write(6,*)  sum(this%R(it,:))/this%nwalkers
        enddo

        write(6,*)
        write(6,*)  '----                           ____'
    end subroutine
end module mcwalkers
