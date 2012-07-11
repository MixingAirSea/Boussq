      program threedfluid
       use threedstuff
       implicit none
       integer :: Nx,Ny,Nz,Nme,Ncpu,ierr
       real :: Re,Ri,Pr,hinv,invRo
       character(2) :: thefile
       real, parameter :: pi = 4.0*atan(1.0)
       integer :: reps,i,isave
       real :: tic,toc,cit,cot
       character(14) :: title
       character(15) :: ftitle
!      Input keys
       integer :: ncom
       character(len=32) :: arg
       character(len=2) :: zpoints, ypoints, xpoints
       character(len=2) :: layer, steps, Rich, Roinv
!
        call cpu_time(cit)
        call MPI_INIT(ierr)
        call errcheck

        call MPI_COMM_RANK(MPI_COMM_WORLD,Nme, ierr)
        call errcheck

        call MPI_COMM_SIZE(MPI_COMM_WORLD, Ncpu, ierr)
        call errcheck
!
       layer   = "h0"
       zpoints = "Nz"
       ypoints = "Ny"
       xpoints = "Nx"
       steps   = "Nt"
       Rich    = "Ri"
       Roinv   = "Ro"
!       Default Values
       Nx = 64
       Ny = 64
       Nz = 64
       reps = 360
       hinv = 1.0
       Ri = .25
       invRo = 0.0
!
       ncom = command_argument_count()
       i = 1
       do while(i.lt.ncom)
        call get_command_argument(i, arg)

        if(trim(arg).eq.zpoints) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg,'(i10)') Nz

        elseif(trim(arg).eq.ypoints) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg, '(i10)') Ny

        elseif(trim(arg).eq.xpoints) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg, '(i10)') Nx

        elseif(trim(arg).eq.layer) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg,'(f6.0)') hinv

        elseif(trim(arg).eq.Rich) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg,'(f6.0)') Ri

        elseif(trim(arg).eq.Roinv) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg,'(f6.0)') invRo

        elseif(trim(arg).eq.steps) then
         i = i+1
         call get_command_argument(i, arg)
         read(arg,'(i10)') reps

        end if
        i = i+1

       end do
       !
       thefile = "hi"
!       Ensure hinv >= 1
       hinv = abs(hinv)
       if((hinv<1.0).and.(hinv.ne.0)) then
        hinv=1.0/hinv
       endif
!       Choose Reynolds number so that motion at the Kolmogorov scale
!       is resolved
!       eta = dx/(2*pi)
       Re = 1.0/Nz !  eta
       Re = Re**(-4.0/3.0)
       Pr = 1.0

!       Correct the Reynolds number estimate for Prandtl number effects
       Re = Re/Pr
!       isave = approximate scaled time between snapshots
       isave = 1
!       reps = approximate scaled time of simulation
       reps  = ceiling(reps*5.0*Nz)
       if(reps<4) then
        reps=4
       endif
       isave = ceiling(isave*5.0*Nz)
       ! Create all the working arrays
       call startthreed(thefile,Re,Ri,Pr,invRo,Nx,Ny,Nz)
!       call initdertest(1,2,3)
       call cordwrite(thefile)
!       call stepwrite(thefile)
!       call dertestXX()
!       call stepwrite(thefile)
!       call dertestYY()
!       call stepwrite(thefile)
!1       call dertestZZ()
!       call stepwrite(thefile)
       call initialvals(hinv)
       call stepwrite(thefile)
       call cpu_time(tic)
       do i = 1,reps
        call tstep
        if(modulo(i,isave).eq.0) then
         call stepwrite(thefile)
        end if
       end do
       call cpu_time(toc)
       call stepwrite(thefile)
       call endthreed
       call cpu_time(cot)
       if(Nme.eq.0) print *, Nx,Ny,Nz,reps,Ncpu,tic-cit,toc-cit,cot-cit
      endprogram
