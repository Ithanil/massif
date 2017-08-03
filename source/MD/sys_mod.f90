module sysmod
    use configmod
    !----------------------------------------------------------------------------------!
    !Module of Systems (Prefix: sysm_ )                                                             !
    !Author: Jan Kessler                                                               !
    !Description:                                                                      !
    !This module provides a class suitable for general multi-molecule PIMD simulations.!
    !----------------------------------------------------------------------------------!
    implicit none

    type sysm_atom
        !Represents an atom (ring polymer) with nb beads.'Contained' routines
        !yield an interface to PI module and related observables.

        !name               <->     string consisting of the  single
        !                           character 'name' of the atom
        !mass,name,r,v will point to arrays of a system object
        integer, pointer :: nb
        real(8), pointer :: mass,r(:,:),v(:,:)
        character, pointer :: name

    end type sysm_atom

    type sysm_molecule
        !Represents a molecule of napm atoms with a total mass mass.
        !Molecular COM and Dipole moment interface are 'contained'.
        integer, pointer :: napm,nb
        type(sysm_atom), pointer :: atom(:)
        real(8) :: mass

    end type sysm_molecule

    type sysm_system
        !The type system yields objects meant to represent a group of chemically equal molecules
        !consisting of atoms, each sampled by a nb bead ring polymer.
        !General system observable routines are 'contained'.
        !
        !Variables:
        !- moltype,pitype give information for external force and PI routines
        !- r,v,dvdr store the total configurational data ('containing' simple print routines)
        !- r,v,mass,name and nm,na,napm,nb are pointees of the lower molecule and atom types

        integer, pointer ::  nm,na,napm,nb
        type(sysm_molecule), pointer :: mol(:)
        real(8) :: totmass
        character(*) moltype

        integer trajfilenum,nforces,nslforces

        real(8), pointer   :: r(:,:,:),v(:,:,:),dvdr(:,:,:,:),pot(:)
        real(8), pointer   :: mass(:)
        character, pointer :: name(:)

    end type sysm_system

    integer sysm_nsys
    type(sysm_system), allocatable :: sysm_asys(:)
    integer, pointer :: sysm_anm(:),sysm_ana(:),sysm_anapm(:),sysm_anb(:)
    !because 'target' is not allowed inside a type, the targets
    !are stored as pointers in the module (a list of them)
contains

    !-----------------------------Construct and Destruct-----------------------------!

    subroutine sysm_modinit(nms,napms,nbs,nsys)
        !set configurational integer targets
        implicit none
        integer nms(nsys),napms(nsys),nbs(nsys),nsys

        allocate(sysm_nm(nsys),sysm_na(nsys),sysm_napm(nsys),sysm_nb(nsys))
        sysm_anm(:)=nms(:)
        sysm_anapm(:)=napms(:)
        sysm_ana(:)=nms(:)*napms(:)
        sysm_anb(:)=nbs(:)

    end subroutine

    !initializer cascade sys->mol->atom
    subroutine sysm_sysinit(sys,mass,name,sysnum,nforces,nslforces,trajfilenum)
        !initializes the system with number sysnum
        !the numbers nm,napm,nb are set in the modules dummy target arrays
        !!! DONT use without properly initialized module via sysm_modinit !!
        implicit none
        type(sysm_system) sys
        integer nm,napm,nb,na,sysnum,nforces,nslforces,trajfilenum
        real(8) mass(sysm_napm(sysnum)*sysm_nm(sysnum))
        character name(sysm_napm(sysnum)*sysm_nm(sysnum))

        integer im,ia,na
        nm=sysm_nm(sysnum)
        napm=sysm_napm(sysnum)
        nb=sysm_nb(sysnum)
        na=nm*napm

        sys%nm=>sysm_anm(sysnum)
        sys%napm=>sysm_anapm(sysnum)
        sys%na=>sysm_ana(sysnum)
        sys%nb=>sysm_anb(sysnum)
        sys%nforces=nforces
        sys%nslforces=nslforces
        sys%totmass=sum(mass(:))
        sys%trajfilenum=trajfilenum

        allocate(sys%mol(nm),sys%r(3,na,nb),sys%v(3,na,nb),sys%dvdr(3,na,nb,sys%nforces), &
            sys%pot(sys%nforces),sys%mass(na),sys%name(na))

        sys%r(:,:,:)=0.d0
        sys%v(:,:,:)=0.d0
        sys%dvdr(:,:,:,:)=0.d0

        sys%mass(:)=mass(:)
        sys%name(:)=name(:)

        do im=1,nm
            ia=(im-1)*napm
            call sysm_molinit(sys,im)
        enddo

    end subroutine

    subroutine sysm_molinit(sys,im)
        implicit none
        type(sysm_system) sys
        integer im

        integer iam,ia

        allocate(sys%mol(im)%atom(sys%napm))

        sys%mol(im)%napm => sys%napm
        sys%mol(im)%nb => sys%nb
        sys%mol%mass=sum(sys%mass((im-1)*sys%napm+1 : im*sys%napm))

        iam=(im-1)*sys%napm
        do ia=iam+1,iam+sys%napm
            call sysm_atominit(sys,im,ia-iam,ia)
        enddo

    end subroutine

    subroutine sysm_atominit(sys,im,iainm,ia)
        implicit none
        type(sysm_system) sys
        integer im,iainm,ia

        sys%mol(im)%atom(iainm)%nb   => sys%nb
        sys%mol(im)%atom(iainm)%mass => sys%mass(ia)
        sys%mol(im)%atom(iainm)%name => sys%name(ia)
        sys%mol(im)%atom(iainm)%r    => sys%r(1:3,ia,1:sys%nb)
        sys%mol(im)%atom(iainm)%v    => sys%v(1:3,ia,1:sys%nb)

    end subroutine

    subroutine sysm_dealloc(sys)
        !deallocator cascade
        implicit none

        type(sysm_system) sys
        integer i,j,k

        do i=1,sys%nm
            deallocate(sys%mol(i)%atom)
        enddo
        deallocate(sys%mol)
        deallocate(sys%r,sys%v,sys%dvdr,sys%mass,sys%name)

    end subroutine

    !--------------------------------Get Coordinates---------------------------------!

    subroutine sysm_mcoords(mol,ra,mode)
        !save atomic centroids of mol in ra
        !mode == 1 -> give positions, 2 -> velocities
        implicit none
        type(sysm_molecule) mol
        real(8) ra(3,mol%napm)
        integer mode

        integer ia,i,ib

        ra(:,:)=0.d0
        if (mode==1) then
            do ia=1,mol%napm
                do ib=1,mol%nb
                    ra(:,ia)=ra(:,ia)+mol%atom(ia)%r(:,ib)
                enddo
                ra(:,ia)=ra(:,ia)/mol%nb
            enddo
        else
            do ia=1,mol%napm
                do ib=1,mol%nb
                    ra(:,ia)=ra(:,ia)+mol%atom(ia)%v(:,ib)
                enddo
                ra(:,ia)=ra(:,ia)/mol%nb
            enddo
        end if

    end subroutine

    subroutine sysm_mcom(mol,rcom,mode)
        !save molecular com of mol in rcom
        !mode == 1 -> give positions, 2 -> velocities
        implicit none
        type(sysm_molecule) mol
        real(8) rcom(3)
        integer mode

        real(8) ra(3,mol%napm)
        integer ia

        call sysm_mcoords(mol,ra,mode)

        rcom(:)=0.d0
        do ia=1,mol%napm
            rcom(:)=rcom(:)+mol%atom(ia)%mass*ra(:,ia)
        enddo
        rcom(:)=rcom(:)/mol%mass

    end subroutine

    subroutine sysm_allcentroids(sys,ra,mode)
        !saves all atomic centroids of sys in ra
        !mode == 1 -> give positions, 2 -> velocities
        implicit none
        type(sysm_system) sys
        real(8) ra(3,sys%na)
        integer mode

        integer ib,ia,im,i

        ra(:,:)=0.d0
        if (mode==1) then
            do ib=1,sys%nb
                ra(:,:)=ra(:,:)+sys%r(:,:,ib)
            enddo
        else
            do ib=1,sys%nb
                ra(:,:)=ra(:,:)+sys%r(:,:,ib)
            enddo
        end if
        ra(:,:)=ra(:,:)/sys%nb

    end subroutine

    subroutine sysm_allmcoms(sys,rm,mode)
        !saves all molecular com of sys in rm
        !mode == 1 -> give positions, 2 -> velocities
        implicit none
        type(sysm_system) sys
        real(8) rm(3,sys%nm),ra(3,sys%na)
        integer mode

        integer k,im,ia

        call sysm_allcentroids(sys,ra,mode)

        rm(:,:)=0.d0
        ia=0
        do im=1,sys%nm
            do k=1,sys%napm
                ia=ia+1
                rm(:,im)=rm(:,im)+sys%mass(ia)*ra(:,ia)
            enddo
            rm(:,im)=rm(:,im)/sys%mol(im)%mass
        enddo

    end subroutine


    !-------------------------------Reading, Printing and Tools-------------------------------!
    !subroutine sysm_readin(sys)
    !    !reads in a single time step from sys%trajfilenum, if not 0
    !    implicit none
    !    integer nmr,nar,nbr
    !    real(8) boxlxyz_traj(3)
    !
    !    if (sys%trajfilenum.ne.0) then
    !        read(sys%trajfilenum,*) nmr,nar,nbr
    !        read(sys%trajfilenum,*) cfm_boxlxyz(:)
    !
    !
    !        if (sys%nm .ne. nmr) then
    !            write(6,*) 'Error in sysm_readin: sys%na != nar'
    !        else if (sys%na .ne. nar) then
    !            write(6,*) 'Error in sysm_readin: sys%na != nar'
    !        else if (nbr.eq.sys%nb) then
    !
    !            read(sys%trajfilenum,*) sys%r(:,:,:)
    !
    !        else if (nbr.eq.1) then
    !            ! This allows reading of classical equilibrated
    !            ! configurations to start PI trajectories.
    !
    !            write(6,*) ' Note:// all beads will start at centroid'
    !
    !            ! Read first bead coordinates
    !
    !            read(61,*) sys%r(:,:,1)
    !            ! Copy to all other beads
    !
    !            do k = 2,sys%nb
    !                do j = 1,sys%na
    !                    sys%r(:,j,k) = sys%r(:,j,1)
    !                enddo
    !            enddo
    !
    !        else
    !            write(6,*) 'Error in sysm_readin: Only sys%nb==nbr or nbr==1 allowed! Stop.'
    !            stop
    !        end if
    !
    !        if (cfm_doreftraj) then
    !            write(6,*)'* Initialized from file: ',cfm_initfile
    !        end if
    !    else
    !        write(6,*) 'Error in sysm_readin:   trajfilenum not present!. Stop.!'
    !        stop
    !    end if
    !
    !end subroutine

    subroutine sysm_printme(sys,channel,pr,pv)
        implicit none
        type(sysm_system) sys
        real(8) rc(3,sys%na),vc(3,sys%na)
        logical pr,pv !print r? print v?
        integer channel

        integer i,j,k

        if (pr) then
            call sysm_allcentroids(sys,rc,1)
            do i=1,sys%na
                write(channel,*) sys%name(i),rc(:,i)
            enddo
        end if
        if (pv) then
            call sysm_allcentroids(sys,vc,2)
            do i=1,sys%na
                write(channel,*) sys%name(i),vc(:,i)
            enddo
        end if

    end subroutine


    subroutine sysm_centroid(vec,cvec,nb)
        implicit none
        real(8) vec(3,nb),cvec(3)
        integer nb

        integer ib
        cvec(:)=0.d0
        do ib=1,nb
            cvec(:)=cvec(:)+vec(:,ib)
        enddo
        cvec(:)=cvec(:)/nb
    end subroutine


end module sysmod

