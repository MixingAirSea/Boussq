      module threedstuff
       use threedNetCDF
       implicit none
       include 'mpif.h'
       include 'fftw3.f'
       integer,private :: Nx,Ny,Nz
       integer,private :: Nxf
       integer,private :: Mx,My,Mz
       integer,private :: Mxf
       integer,private :: ax,az,axf,bufsize
       integer,private :: Ncpu,Nme,ierr,rc,req
       integer, dimension(MPI_STATUS_SIZE) :: mpistat
       real, private :: Lx,Ly,Lz,Re,Ri,invRo,Pr,dt,T
       real, private, allocatable, dimension(:) :: x,y,z,k,l,msin,mcos
       complex, private, allocatable, dimension(:,:,:) :: ufro,vfro,wfro,bfro,zfro
       complex, private, allocatable, dimension(:,:,:) :: dufro,dvfro,dwfro,dbfro
       complex, private, allocatable, dimension(:,:,:) :: uvfro, vvfro, vwfro, wwfro
       complex, private, allocatable, dimension(:,:,:) :: vbfro, wbfro
       complex, private, allocatable, dimension(:,:,:) :: uzufro, uzvfro
       complex, private, allocatable, dimension(:,:,:) :: uzwfro, uzbfro
       real, private, allocatable, dimension(:,:,:) :: esinfro, ecosfro,eboufro
       real, private, allocatable, dimension(:,:,:) :: ucol,vcol,wcol,bcol
       real, private, allocatable, dimension(:) :: sendbuf, recvbuf
       complex, private, allocatable, dimension(:) :: slabin,slabout
       real, private, dimension(:), allocatable :: zfft
       real, private, dimension(:,:), allocatable :: xyfft
       complex, private, dimension(:,:), allocatable :: klfft
       real, private, dimension(:,:,:), allocatable :: ksinsq, kcossq
       integer(8), private :: fwddftxyfft, invdftxyfft
       integer(8), private :: fwdsinzfft, invsinzfft
       integer(8), private :: fwdcoszfft, invcoszfft
       contains

       subroutine startthreed(thefile,Rein,Riin,Prin,Roin,Nxin,Nyin,Nzin)
        character(*), intent(in) :: thefile
        integer, intent(in) :: Nxin,Nyin,Nzin
        real, intent(in) :: Rein,Riin,Prin
        real, intent(in) :: Roin
        integer :: i

!        call MPI_INIT(ierr)
        call errcheck

        call MPI_COMM_RANK(MPI_COMM_WORLD,Nme, ierr)
        call errcheck

        call MPI_COMM_SIZE(MPI_COMM_WORLD, Ncpu, ierr)
        call errcheck

        call srand(Nme)

        call arraysizes(Nxin,Nyin,Nzin)

        Lz = 2.0

        Lx = Nx*Lz/Nz

        Ly = Ny*Lz/Nz

        Re = Rein

        Ri = Riin

        Pr = Prin

        invRo = Roin

        T = 0.0

        dt = .1*Lz/Nz
        call makeworkingarrays
        if(Nme.eq.0) then
         call newfile(thefile,Re,Ri,Pr,invRo,Nx,Ny,Nz,Lz)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck

        call fillcords

        call startfftw
       end subroutine

       subroutine initdertest(xmode,ymode,zmode)
        integer :: h,i,j
        integer, intent(in):: xmode,ymode,zmode
        real :: pi
        pi = 4*atan(1.0)
        vcol = 0
        wcol = 0
        ucol = 0
        bcol = 0
        do h=1,Nz
         do i=1,Ny
          do j=1,ax
           ucol(h,i,j) = cos(xmode*x(j)*pi*2/Lx)
           ucol(h,i,j) = ucol(h,i,j)*cos(ymode*y(i)*pi*2/Ly)
           bcol(h,i,j) = ucol(h,i,j)*sin(zmode*pi*(.5+z(h)/Lz))
           ucol(h,i,j) = ucol(h,i,j)*cos(zmode*pi*(.5+z(h)/Lz))
           !+(rand()-.5)/(Nx*Ny*Nz)
           end do
         end do
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        call fwdfft
        T = T-dt/2
       endsubroutine

       subroutine dertestXX()
        call ddx(ufro,vfro)
        call ddx(bfro,wfro)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        call invfft
        T = T+dt/8
       endsubroutine

       subroutine dertestYY()
        call ddy(ufro,vfro)
        call ddy(bfro,wfro)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        call invfft
        T = T+dt/8
       endsubroutine

       subroutine dertestZZ()
        call Fddcosz(ufro,wfro)
        call Fddsinz(bfro,vfro)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        call invfft
        T = T+dt/8
       end subroutine

       subroutine tstep

         call ddtfwd

        ufro = ecosfro*(ufro+dufro*dt)
        vfro = ecosfro*(vfro+dvfro*dt)
        wfro = esinfro*(wfro+dwfro*dt)
        bfro = eboufro*(bfro+dbfro*dt)

        call nodiv

        call ddtmod(dufro,dvfro,dwfro,dbfro)

        ufro = ufro+ecosfro*dufro*dt
        vfro = vfro+ecosfro*dvfro*dt
        wfro = wfro+esinfro*dwfro*dt
        bfro = bfro+eboufro*dbfro*dt

        call nodiv

        T = T+dt

       end subroutine

       subroutine ddtfwd()
        complex, allocatable,dimension(:,:,:) :: workfro
        allocate(workfro(Mxf,Ny,az))
        call invfft
        call fwdauxfft

        call ddx(uzufro,workfro)
        dufro = -workfro
        call ddy(uvfro,workfro)
        dufro = dufro-workfro
        call Fddsinz(uzwfro,workfro)
        dufro = dufro-workfro
        dufro = dufro+invRo*vfro

        call ddx(uzvfro,workfro)
        dvfro = -workfro
        call ddy(vvfro,workfro)
        dvfro = dvfro-workfro
        call Fddsinz(vwfro,workfro)
        dvfro = dvfro-workfro
        dvfro = dvfro-invRo*(ufro+zfro)

        call ddx(uzwfro,workfro)
        dwfro = -workfro
        call ddy(vwfro,workfro)
        dwfro = dwfro-workfro
        call Fddcosz(wwfro,workfro)
        dwfro = dwfro-workfro
        dwfro = dwfro+Ri*bfro

        call ddx(uzbfro,workfro)
        dbfro = -workfro
        call ddy(vbfro,workfro)
        dbfro = dbfro-workfro
        call Fddcosz(wbfro,workfro)
        dbfro = dbfro-workfro
        dbfro = dbfro- wfro

        deallocate(workfro)

       end subroutine

       subroutine ddtmod(duin,dvin,dwin,dbin)
        complex,dimension(:,:,:),intent(inout)::duin,dvin,dwin,dbin
        complex, allocatable, dimension(:,:,:) :: dtold, ddnext
        call invfft
        call fwdauxfft
        allocate(dtold(Mxf,Ny,az))
        allocate(ddnext(Mxf,Ny,az))

        dtold = duin
        call ddx(uzufro,ddnext)
        duin = -ddnext
        call ddy(uvfro,ddnext)
        duin = duin-ddnext
        call Fddsinz(uzwfro,ddnext)
        duin = duin-ddnext
        duin = duin+invRo*vfro
        duin = (duin-dtold)/2

        dtold = dvin
        call ddx(uzvfro,ddnext)
        dvin = -ddnext
        call ddy(vvfro,ddnext)
        dvin = dvin-ddnext
        call Fddsinz(vwfro,ddnext)
        dvin = dvin-ddnext
        dvin = dvin-invRo*(ufro+zfro)
        dvin = (dvin-dtold)/2

        dtold = dwin
        call ddx(uzwfro,ddnext)
        dwin = -ddnext
        call ddy(vwfro,ddnext)
        dwin = dwin-ddnext
        call Fddcosz(wwfro,ddnext)
        dwin = dwin-ddnext
        dwin = dwin+Ri*bfro
        dwin = (dwin-dtold)/2

        dtold = dbin
        call ddx(uzbfro,ddnext)
        dbin = -ddnext
        call ddy(vbfro,ddnext)
        dbin = dbin-ddnext
        call Fddcosz(wbfro,ddnext)
        dbin = dbin-ddnext
        dbin = dbin-wfro
        dbin = (dbin-dtold)/2

        deallocate(dtold)
        deallocate(ddnext)
       end subroutine

       subroutine cordwrite(filename)
        integer :: i,j,recv,f,g,h
        character(*), intent(in) :: filename
        call invfft
        j = modulo(Nme+1,Ncpu)
        i = modulo(Nme-1,Ncpu)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.ne.0) then
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
        end if
        call cordsave(filename,X,Y,Z,Nme*ax+1,ax)
        call MPI_SEND(Nme,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.eq.0) then
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       end subroutine


       subroutine stepwrite(filename)
        integer :: i,j,recv,f,g,h
        character(*), intent(in) :: filename
        call invfft
        j = modulo(Nme+1,Ncpu)
        i = modulo(Nme-1,Ncpu)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.ne.0) then
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
        end if
        call stepsave(filename,Ucol,Vcol,Wcol,bcol,t,Nme*ax+1)
        call MPI_SEND(Nme,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.eq.0) then
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       end subroutine


       subroutine altwrite(filename)
        integer :: i,j,recv,f,g,h
        character(*), intent(in) :: filename
        call invfft
        j = modulo(Nme+1,Ncpu)
        i = modulo(Nme-1,Ncpu)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.eq.0) then
         open(1,file=filename,action='write',status='replace')
        else
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
         open(1,file=filename,action='write',status='old',position='append')
        end if
        do g = 1,min(ax,Nx-Nme*ax)
         do h = 1,Ny
          do f = 1,Nz
           write(1,*) x(g),y(h),z(f),ucol(f,h,g),vcol(f,h,g),wcol(f,h,g),bcol(f,h,g),2*z(f)/Lz
          end do
          write(1,*) ""
         end do
         write(1,*) ""
        end do
        close(1)
        call MPI_SEND(Nme,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,ierr)
        call errcheck
        if(Nme.eq.0) then
         call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
         call errcheck
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       end subroutine

       subroutine Fringwrite(filename)
        integer :: i,j,recv,f,g,h
        character(*), intent(in) :: filename
        call invfft
         j = modulo(Nme+1,Ncpu)
         i = modulo(Nme-1,Ncpu)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call errcheck
         if(Nme.eq.0) then
          open(1,file=filename,action='write',status='replace')
          do g = 1,az
           do h = 1,Ny
            do f = 1,Nxf
             write(1,*) k(f),l(h),mcos(g),msin(g),abs(ufro(f,h,g)),abs(vfro(f,h,g)),abs(wfro(f,h,g)),abs(bfro(f,h,g))
            end do
            write(1,*) ""
           end do
           write(1,*) ""
          end do
          close(1)
          call MPI_ISEND(Nme,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,req,ierr)
          call errcheck
          call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          call MPI_WAIT(req,mpistat,ierr)
          call errcheck
         else
          call MPI_RECV(recv,1,MPI_INTEGER,i,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          open(1,file=filename,action='write',status='old',position='append')
          do g = 1,min(az,Nz-az*Nme)
           do h = 1,Ny
            do f = 1,Nxf
             write(1,*) k(f),l(h),mcos(g),msin(g),abs(ufro(f,h,g)),abs(vfro(f,h,g)),abs(wfro(f,h,g)),abs(bfro(f,h,g))
            end do
            write(1,*) ""
           end do
           write(1,*) ""
          end do
          close(1)
          call MPI_SEND(Nme,1,MPI_INTEGER,j,j,MPI_COMM_WORLD,ierr)
          call errcheck
         end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       end subroutine

       integer function myrank()
        myrank = Nme
        return
       end function

       subroutine startfftw
        call dfftw_plan_dft_r2c_2d(fwddftxyfft,Nx,Ny,xyfft,klfft,FFTW_PATIENT)
        call dfftw_plan_dft_c2r_2d(invdftxyfft,Nx,Ny,klfft,xyfft,FFTW_PATIENT)
        call dfftw_plan_r2r_1d(fwdsinzfft,Nz,zfft,zfft,FFTW_RODFT10,FFTW_PATIENT)
        call dfftw_plan_r2r_1d(invsinzfft,Nz,zfft,zfft,FFTW_RODFT01,FFTW_PATIENT)
        call dfftw_plan_r2r_1d(fwdcoszfft,Nz,zfft,zfft,FFTW_REDFT10,FFTW_PATIENT)
        call dfftw_plan_r2r_1d(invcoszfft,Nz,zfft,zfft,FFTW_REDFT01,FFTW_PATIENT)
       end subroutine

!> Converts (U,V,W,B) data into Fourier space
       subroutine fwdfft
        real, allocatable, dimension(:,:,:) :: workrow
        allocate(workrow(Mx,Ny,az))

        call fwdcosz(ucol)
        call col2row(ucol,workrow)
        call fwdxyfft(workrow,ufro)

        call fwdcosz(vcol)
        call col2row(vcol,workrow)
        call fwdxyfft(workrow,vfro)

        call fwdsinz(wcol)
        call col2row(wcol,workrow)
        call fwdxyfft(workrow,wfro)

        call fwdsinz(bcol)
        call col2row(bcol,workrow)
        call fwdxyfft(workrow,bfro)
        deallocate(workrow)
       end subroutine
!> Converts Z into Fourier space
       subroutine makezfro
        integer :: i,j
        real, allocatable, dimension(:,:,:) :: workcol
        real, allocatable, dimension(:,:,:) :: workrow
        allocate(workrow(Mx,Ny,az))
        allocate(workcol(Mz,Ny,ax))

        do i=1,Ny
         do j=1,ax
          workcol(1:Nz,i,j) = z(1:Nz)
         enddo
        enddo
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,zfro)
        deallocate(workcol)
        deallocate(workrow)
       end subroutine

!> Converts data into Fourier space
       subroutine fwdauxfft
        integer :: i,j
        real, allocatable, dimension(:,:,:) :: workcol
        real, allocatable, dimension(:,:,:) :: workrow
        allocate(workrow(Mx,Ny,az))
        allocate(workcol(Mz,Ny,ax))

        do i=1,Ny
         do j=1,ax
          workcol(1:Nz,i,j) = (ucol(1:Nz,i,j)+z(1:Nz))
         enddo
        enddo
        workcol = workcol*workcol
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,uzufro)

        do i=1,Ny
         do j=1,ax
          workcol(1:Nz,i,j) = vcol(1:Nz,i,j)*(ucol(1:Nz,i,j)+z(1:Nz))
         enddo
        enddo
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,uzvfro)

        workcol = vcol*vcol
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,vvfro)

        do i=1,Ny
         do j=1,ax
          workcol(1:Nz,i,j) = wcol(1:Nz,i,j)*(ucol(1:Nz,i,j)+z(1:Nz))
         enddo
        enddo
        call fwdsinz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,uzwfro)

        workcol = wcol*wcol
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,wwfro)

        workcol = vcol*ucol
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,uvfro)

        workcol = vcol*wcol
        call fwdsinz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,vwfro)

        do i=1,Ny
         do j=1,ax
          workcol(1:Nz,i,j) = bcol(1:Nz,i,j)*(ucol(1:Nz,i,j)+z(1:Nz))
         enddo
        enddo
        call fwdsinz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,uzbfro)

        workcol = vcol*bcol
        call fwdsinz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,vbfro)

        workcol = wcol*bcol
        call fwdcosz(workcol)
        call col2row(workcol,workrow)
        call fwdxyfft(workrow,wbfro)

        deallocate(workcol)
        deallocate(workrow)
       end subroutine

!> Convert data into real space
       subroutine invfft
        real, allocatable, dimension(:,:,:) :: workrow
        allocate(workrow(Mx,Ny,az))
        call invxyfft(ufro,workrow)
        call row2col(workrow,ucol)
        call invcosz(ucol)

        call invxyfft(vfro,workrow)
        call row2col(workrow,vcol)
        call invcosz(vcol)

        call invxyfft(wfro,workrow)
        call row2col(workrow,wcol)
        call invsinz(wcol)

        call invxyfft(bfro,workrow)
        call row2col(workrow,bcol)
        call invsinz(bcol)

        deallocate(workrow)
       end subroutine

       subroutine fwdsinz(A)
        real, intent(inout), dimension(:,:,:) :: A
        integer :: i,j
        do i=1,ax
         do j=1,Ny
          zfft = A(1:Nz,j,i)
          call dfftw_execute(fwdsinzfft)
          A(1:Nz,j,i) = zfft/(2*Nz)
         end do
        end do
       end subroutine

       subroutine invsinz(A)
        real, intent(inout), dimension(:,:,:) :: A
        integer :: i,j
        do i=1,ax
         do j=1,Ny
          zfft = A(1:Nz,j,i)
          call dfftw_execute(invsinzfft)
          A(1:Nz,j,i) = zfft
         end do
        end do
       end subroutine

       subroutine fwdcosz(A)
        real, intent(inout), dimension(:,:,:) :: A
        integer :: i,j
        do i=1,ax
         do j=1,Ny
          zfft = A(1:Nz,j,i)
          call dfftw_execute(fwdcoszfft)
          A(1:Nz,j,i) = zfft/(2*Nz)
         end do
        end do
       end subroutine

       subroutine invcosz(A)
        real, intent(inout), dimension(:,:,:) :: A
        integer :: i,j
        do i=1,ax
         do j=1,Ny
          zfft = A(1:Nz,j,i)
          call dfftw_execute(invcoszfft)
          A(1:Nz,j,i) = zfft
         end do
        end do
       end subroutine

       subroutine fwdxyfft(R,F)
        real, intent(in), dimension(:,:,:) :: R
        complex, intent(out), dimension(:,:,:) :: F
        integer :: i
         do i=1,az
          xyfft = R(1:Nx,:,i)
          call dfftw_execute(fwddftxyfft)
          F(1:Nxf,:,i) = klfft/(Nx*Ny)
         end do
       end subroutine

       subroutine invxyfft(F,R)
        real, intent(out), dimension(:,:,:) :: R
        complex, intent(in), dimension(:,:,:) :: F
        integer :: i
         do i=1,az
          klfft = F(1:Nxf,:,i)
          call dfftw_execute(invdftxyfft)
          R(1:Nx,:,i) = xyfft
         end do
       end subroutine

       subroutine ddx(fin,fout)
        complex, intent(in), dimension(:,:,:) :: fin
        complex, intent(out), dimension(:,:,:) :: fout
        integer :: i,j
        forall(i=1:Ny,j=1:az)
         fout(:,i,j) = (0,1)*k(:)*fin(:,i,j)
        end forall
       end subroutine

       subroutine ddy(fin,fout)
        complex, intent(in), dimension(:,:,:) :: fin
        complex, intent(out), dimension(:,:,:) :: fout
        integer :: i,j
        forall(i=1:Mxf,j=1:az)
         fout(i,:,j) = (0,1)*l(:)*fin(i,:,j)
        end forall
       end subroutine

       subroutine Fddsinz(fin,fout)
        complex, intent(in), dimension(:,:,:) :: fin
        complex, intent(out), dimension(:,:,:) :: fout
        integer :: j, send,recv, nyq
        send = modulo(Nme+1,Ncpu)
        recv = modulo(Nme-1,Ncpu)
        nyq = min(az,Nz-Nme*az)
        fout = (0.0,0.0)
        slabout = reshape(fin(:,:,nyq),[Mxf*Ny])
        call MPI_ISEND(slabout,Mxf*Ny,MPI_DOUBLE_COMPLEX,send,send,MPI_COMM_WORLD,req,ierr)
        call errcheck
        forall(j=2:az)
         fout(:,:,j) = mcos(j)*fin(:,:,j-1)
        end forall
        call MPI_RECV(slabin,Mxf*Ny,MPI_DOUBLE_COMPLEX,recv,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
        call errcheck
        fout(:,:,1) = mcos(1)*reshape(slabin,[Mxf,Ny])
        call MPI_WAIT(req,mpistat,ierr)
        call errcheck
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       endsubroutine

       subroutine Fddcosz(fin,fout)
        complex, intent(in), dimension(:,:,:) :: fin
        complex, intent(out), dimension(:,:,:) :: fout
        integer :: j, send,recv, nyq
        send = modulo(Nme-1,Ncpu)
        recv = modulo(Nme+1,Ncpu)
        nyq = min(az,Nz-Nme*az)
        fout = 0
        slabout = -mcos(1)*reshape(fin(:,:,1),[Mxf*Ny])
        call MPI_ISEND(slabout,Mxf*Ny,MPI_DOUBLE_COMPLEX,send,send,MPI_COMM_WORLD,req,ierr)
        call errcheck
        forall(j=2:az)
         fout(:,:,j-1) = -mcos(j)*fin(:,:,j)
        end forall
        call MPI_RECV(slabin,Mxf*Ny,MPI_DOUBLE_COMPLEX,recv,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
        call errcheck
        fout(:,:,nyq) = reshape(slabin,[Mxf,Ny])
        call MPI_WAIT(req,mpistat,ierr)
        call errcheck
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
       endsubroutine

!> Takes as argument the inverse length scale
!> of the stratification and velocity profiles
       subroutine initialvals(hinv)
        real, intent(in) :: hinv
        integer :: h,i,j
        real :: pi,dbdz
        pi = 4*atan(1.0)
        vcol = 0
        wcol = 0
        ucol = 0
        bcol = 0
        if(hinv.eq.0) then
         do h=1,Nz
          do i=1,Ny
           do j=1,ax
            ucol(h,i,j) = -z(h)
            if(Ri.ne.0) then
             bcol(h,i,j) =-z(h)+(rand()-.5)/Nz
            else
             wcol(h,i,j) = (rand()-.5)/Nz
            end if
           end do
          end do
         end do

        elseif(hinv.eq.1.0) then
         do h=1,Nz
          do i=1,Ny
           do j=1,ax
            ucol(h,i,j) = 0
            if(Ri.ne.0) then
             bcol(h,i,j) = (rand()-.5)/Nz
            else
             wcol(h,i,j) = (rand()-.5)/Nz
            end if
           end do
          end do
         end do

        else
         do h=1,Nz
          do i=1,Ny
           do j=1,ax
            ucol(h,i,j) = -z(h)+tanh(z(h)*hinv)
            if(Ri.ne.0) then
             bcol(h,i,j) = -z(h)+tanh(z(h)*hinv)
             bcol(h,i,j) = bcol(h,i,j)+(rand()-.5)*hinv/(Nz*cosh(z(h)*hinv)**2)
            else
             wcol(h,i,j) = (rand()-.5)/Nz
            end if
           end do
          end do
         end do
        endif
        call fwdfft
        call nodiv
        T=0
       end subroutine

       subroutine dumpcols
        integer :: i,j
         do i=1,Ny
          do j=1,ax
           print *, ucol(:,i,j)
           print *, ""
           print *, vcol(:,i,j)
           print *, ""
           print *, wcol(:,i,j)
           print *, ""
           print *, ""
          end do
         end do
       end subroutine

        subroutine row2col(row,col)
         real, intent(in), dimension(:,:,:) :: row
         real, intent(out), dimension(:,:,:) :: col
         integer :: recv,send,i
         do i=1,Ncpu
          recv = modulo(Nme+i,Ncpu)
          send = modulo(Nme-i,Ncpu)
          call fillbufR2C(row,sendbuf,recv)
          call MPI_ISEND(sendbuf,bufsize,MPI_DOUBLE_PRECISION,recv,recv,MPI_COMM_WORLD,req,ierr)
          call errcheck
          call MPI_RECV(recvbuf,bufsize,MPI_DOUBLE_PRECISION,send,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          call readbufR2C(col,recvbuf,send)
          call MPI_WAIT(req,mpistat,ierr)
          call errcheck
         end do
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call errcheck
        end subroutine

        subroutine col2row(col,row)
         real, intent(in), dimension(:,:,:) :: col
         real, intent(out), dimension(:,:,:) :: row
         integer :: recv,send,i
         do i=1,Ncpu
          recv = modulo(Nme+i,Ncpu)
          send = modulo(Nme-i,Ncpu)
          call fillbufC2R(col,sendbuf,recv)
          call MPI_ISEND(sendbuf,bufsize,MPI_DOUBLE_PRECISION,recv,recv,MPI_COMM_WORLD,req,ierr)
          call errcheck
          call MPI_RECV(recvbuf,bufsize,MPI_DOUBLE_PRECISION,send,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          call readbufC2R(row,recvbuf,send)
          call MPI_WAIT(req,mpistat,ierr)
          call errcheck
         end do
        end subroutine

       subroutine fillbufR2C(A,chunk,r)
        real, intent(in), dimension(:,:,:) :: A
        real, intent(out), dimension(:) :: chunk
        integer, intent(in) :: r
        chunk = reshape(A(1+r*ax:(r+1)*ax,:,:),[ax*Ny*az])
       end subroutine

       subroutine readbufR2C(A,chunk,s)
        real, intent(inout), dimension(:,:,:) :: A
        real, intent(in), dimension(:) :: chunk
        integer, intent(in) :: s
        integer :: i,j
          A(s*az+1:(s+1)*az,:,:) = reshape(chunk,[az,Ny,ax],[0.0],[3,2,1])
       end subroutine

       subroutine fillbufC2R(A,chunk,r)
        real, intent(in), dimension(:,:,:) :: A
        real, intent(out), dimension(:) :: chunk
        integer, intent(in) :: r
        chunk = reshape(A(1+r*az:(r+1)*az,:,:),[ax*Ny*az])
       end subroutine

       subroutine readbufC2R(A,chunk,s)
        real, intent(inout), dimension(:,:,:) :: A
        real, intent(in), dimension(:) :: chunk
        integer, intent(in) :: s
        integer :: i,j
        A((s*ax+1):(s+1)*ax,:,:) = reshape(chunk,[ax,Ny,az],[0.0],[3,2,1])
       end subroutine

       subroutine makecols
        allocate(ucol(Mz,Ny,ax))
        allocate(vcol(Mz,Ny,ax))
        allocate(wcol(Mz,Ny,ax))
        allocate(bcol(Mz,Ny,ax))
       end subroutine

       subroutine killcols
        deallocate(ucol)
        deallocate(vcol)
        deallocate(wcol)
        deallocate(bcol)
       end subroutine

       subroutine makeFcols
        allocate(ufro(Mxf,Ny,az))
        allocate(vfro(Mxf,Ny,az))
        allocate(wfro(Mxf,Ny,az))
        allocate(bfro(Mxf,Ny,az))
        allocate(zfro(Mxf,Ny,az))

!        allocate(uwfro(Mxf,Ny,az))
!        allocate(wflipfro(Mxf,Ny,az))

        allocate(dufro(Mxf,Ny,az))
        allocate(dvfro(Mxf,Ny,az))
        allocate(dwfro(Mxf,Ny,az))
        allocate(dbfro(Mxf,Ny,az))

        allocate(vvfro(Mxf,Ny,az))
        allocate(wwfro(Mxf,Ny,az))
        allocate(uvfro(Mxf,Ny,az))
        allocate(vwfro(Mxf,Ny,az))

        allocate(vbfro(Mxf,Ny,az))
        allocate(wbfro(Mxf,Ny,az))

        allocate(uzufro(Mxf,Ny,az))
        allocate(uzvfro(Mxf,Ny,az))
        allocate(uzwfro(Mxf,Ny,az))
        allocate(uzbfro(Mxf,Ny,az))

       end subroutine

       subroutine killFcols
        deallocate(ufro)
        deallocate(vfro)
        deallocate(wfro)
        deallocate(bfro)

!        deallocate(uwfro)
!        deallocate(wflipfro)

        deallocate(dufro)
        deallocate(dvfro)
        deallocate(dwfro)
        deallocate(dbfro)

        deallocate(vvfro)
        deallocate(wwfro)
        deallocate(uvfro)
        deallocate(vwfro)

        deallocate(vbfro)
        deallocate(wbfro)

        deallocate(uzufro)
        deallocate(uzvfro)
        deallocate(uzwfro)
        deallocate(uzbfro)

        deallocate(esinfro)
        deallocate(eboufro)
        deallocate(ecosfro)
       end subroutine


       subroutine makeworkingarrays


        !Columns
        call makecols
        !Fourier Space Rows
        call makeFcols

        ! Send and Receive Buffers
        allocate(sendbuf(ax*Ny*az))
        allocate(recvbuf(ax*Ny*az))
        allocate(slabin(Mxf*Ny))
        allocate(slabout(Mxf*Ny))

        !FFTW working arrays
        allocate(xyfft(Nx,Ny))
        allocate(klfft(Nxf,Ny))
        allocate(zfft(Nz))

       end subroutine

       subroutine endthreed
        call MPI_FINALIZE(ierr)
        call errcheck
       end subroutine


       subroutine arraysizes(Nxin,Nyin,Nzin)
        integer, intent(in) :: Nxin,Nyin,Nzin

        Nx = Nxin
        Ny = Nyin
        Nz = Nzin
        Nxf = Nx/2+1

        if (modulo(Nx,Ncpu).eq.0) then
         Mx = Nx
        else
         Mx = (1+Nx/Ncpu)*Ncpu
        end if

        ax = Mx/Ncpu

        if (modulo(Nz,Ncpu).eq.0) then
         Mz = Nz
        else
         Mz = (1+Nz/Ncpu)*Ncpu
        end if

        az = Mz/Ncpu

!       Calculate the dimension of the Fourier trasform
        if(modulo(Nxf,Ncpu).eq.0) then
         Mxf = Nxf
        else
         Mxf = (1+Nxf/Ncpu)*Ncpu
        end if

        axf = Mxf/Ncpu

        bufsize = ax*Ny*az

       end subroutine

       subroutine fillcords
        allocate(x(ax))
        call fillx

        allocate(y(Ny))
        call filly

        allocate(z(Mz))
        call fillz

        allocate(k(Mxf))
        call fillk()

        allocate(l(Ny))
        call filll()

        allocate(msin(az))
        allocate(mcos(az))
        call fillm()

        allocate(ksinsq(Nxf,Ny,az))
        allocate(kcossq(Nxf,Ny,az))
        call fillksq()

        allocate(esinfro(Mxf,Ny,az))
        allocate(eboufro(Mxf,Ny,az))
        allocate(ecosfro(Mxf,Ny,az))
        call fillfactor()

       end subroutine

       subroutine fillfactor
        integer :: h,i,j
        esinfro = 0
        eboufro = 0
        ecosfro = 0
        do h = 1,Nxf
         do i = 1,Ny
          do j = 1,az
            esinfro(h,i,j) = exp(-dt*ksinsq(h,i,j)/Re)
            eboufro(h,i,j) = exp(-dt*ksinsq(h,i,j)/(Pr*Re))
            ecosfro(h,i,j) = exp(-dt*ksinsq(h,i,j)/Re)
          end do
         end do
        end do
       end subroutine

       subroutine fillx
        integer :: i
        x = 0
        do i=1,ax
         x(i) = Lx*((Nme*ax+i-0.5)/Nx-0.5)
        end do
       end subroutine

       subroutine filly
        integer :: i
        do i=1,Ny
         y(i) = Ly*((i-0.5)/Ny-0.5)
        end do
       end subroutine

       subroutine fillz
        integer :: i
        z = 0
        do i=1,Nz
         z(i) =  Lz*((i-0.5)/Nz-0.5)
        end do
       end subroutine

       subroutine fillk()
        real :: pi
        integer :: i
        pi = 4*atan(1.0)
        k = 0
        do i = 1,Nxf
         k(i) = (i-1)
        end do
        k = k*2*pi/Lx
       end subroutine

       subroutine filll()
        real :: pi
        integer :: i
        pi = 4*atan(1.0)
        l = 0
        do i = 1,Ny
         l(i) = i-1-Ny*floor(2.0*(i-1.0)/Ny)
        end do
         l = l*2*pi/Ly
       end subroutine

       subroutine fillm()
        real :: pi
        integer :: i
        pi = 4*atan(1.0)
        do i = 1,az
         if(i+Nme*az.le.Nz) then
          msin(i) = i+az*Nme
          mcos(i) = i+az*Nme-1
         else
          msin(i) = 0
          mcos(i) = 0
         end if
        end do
          msin = msin*pi/Lz
          mcos = mcos*pi/Lz
       end subroutine

       subroutine fillksq()
        integer :: h,i,j,nyq
         nyq = min(az,Nz-az*Nme)
         ksinsq = 0
         kcossq = 0
        forall(h=1:Nxf,i=1:Ny,j=1:nyq)
         ksinsq(h,i,j) = k(h)**2+l(i)**2+msin(j)**2
         kcossq(h,i,j) = k(h)**2+l(i)**2+mcos(j)**2
        end forall
       end subroutine

       subroutine div(fout)
        complex,intent(out),dimension(:,:,:) :: fout
        complex, allocatable, dimension(:,:,:) :: thediv
        fout = 0
        allocate(thediv(Mxf,Ny,az))
        call ddx(ufro,fout)
        call ddy(vfro,thediv)
        fout = fout+thediv
        call Fddsinz(wfro,thediv)
        fout = fout+thediv
        deallocate(thediv)
       end subroutine

       subroutine showdiv
        integer :: h,i,j
        integer :: hmax,imax,jmax
        real :: testdiv,maxdiv
        complex, allocatable, dimension(:,:,:) :: thediv
        allocate(thediv(Mxf,Ny,az))
        call div(thediv)
        maxdiv = 0
        hmax = 0
        imax = 0
        jmax = 0
        do j = 1,az
         do i = 1,Ny
          do h=1,Nxf
           testdiv = abs(thediv(h,i,j))
           if(maxdiv.le.testdiv) then
            maxdiv = testdiv
            hmax = h
            imax = i
            jmax = j
           end if
          end do
         end do
        end do
         print *, hmax,imax,jmax,maxdiv
         maxdiv = 0
         hmax = 0
         imax = 0
         jmax = 0
         deallocate(thediv)
       end subroutine



       subroutine nodiv
        complex, allocatable, dimension(:,:,:) :: thediv
        complex, allocatable, dimension(:,:,:) :: invdiv

        allocate(thediv(Mxf,Ny,az))

        call div(thediv)

        allocate(invdiv(Mxf,Ny,az))

        call invdelsq(thediv,invdiv)

        call ddx(invdiv,thediv)

        ufro = ufro-thediv

        call ddy(invdiv,thediv)

        vfro = vfro-thediv

        call Fddcosz(invdiv,thediv)

        wfro = wfro-thediv
!        if(Ncpu-Nme.eq.1) wfro(:,:,Nz-az*Nme)=0

        deallocate(thediv)
        deallocate(invdiv)
       end subroutine

       subroutine invdelsq(fin,fout)
        complex, intent(in), dimension(:,:,:) :: fin
        complex, intent(out), dimension(:,:,:) :: fout
        integer :: h,i,j
        fout = 0
        forall(h=1:Nxf,i=1:Ny,j=1:az,kcossq(h,i,j).ne.0)
         fout(h,i,j) = -fin(h,i,j)/kcossq(h,i,j)
        end forall
       end subroutine

       subroutine errcheck()
        if (ierr .ne. MPI_SUCCESS) then
         print *,'MPI program Error', ierr
         print *, 'Terminating'
         call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        end if
       end subroutine

      end module
