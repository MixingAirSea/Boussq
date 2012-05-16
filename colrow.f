! Creates and Transposes a matrix of size equal to the number of processes squared
! This is acomplished by subdividing the matrix into
      module colrow
       implicit none
       include 'mpif.h'
       include 'fftw3.f'
       integer :: numtasks, rank, ierr, rc,req
       integer, dimension(MPI_STATUS_SIZE) :: mpistat
       integer :: chtype
       real, allocatable, dimension(:,:) :: myrows, mycols
       complex, allocatable, dimension(:,:) :: myimagrows
       real, allocatable, dimension(:) :: chin,chout
       integer :: rows,cols,chsize,N,M,nch,mch,Nxf
       real, allocatable, dimension(:) :: zfft,xfft
       complex, allocatable, dimension(:) :: kfft
       integer(8), private :: fwddftxfft, invdftxfft
       integer(8), private :: fwdsinzfft, invsinzfft
       integer(8), private :: fwdcoszfft, invcoszfft

       contains


       subroutine initcolrow(Nx,Nz)
       integer, intent(in) :: Nx,Nz
       call MPI_INIT(ierr)
       call errcheck

       call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
       call errcheck

       call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
       call errcheck

       rows = Nx
       Nxf = 1+Nx/2
       cols = Nz


       if(modulo(rows,numtasks).eq.0) then
        N = rows
       else
        N = (1+rows/numtasks)*numtasks
       endif
       nch = N/numtasks

       if(modulo(cols,numtasks).eq.0) then
        M = cols
       else
        M = (1+cols/numtasks)*numtasks
       endif
       mch = M/numtasks

       chsize = mch*nch

       allocate(xfft(cols))
       allocate(kfft(1+cols/2))
       allocate(zfft(rows))
       allocate(chin(chsize))
       allocate(chout(chsize))
       allocate(myrows(M,nch))
       allocate(myimagrows(1+cols/2,nch))
       allocate(mycols(N,mch))
       mycols = -1
       myrows = 0
       call dfftw_plan_dft_r2c_1d(fwddftxfft,cols,xfft,kfft,FFTW_PATIENT)
       call dfftw_plan_dft_c2r_1d(invdftxfft,cols,kfft,xfft,FFTW_PATIENT)
       call dfftw_plan_r2r_1d(fwdsinzfft,rows,zfft,zfft,FFTW_RODFT10,FFTW_PATIENT)
       call dfftw_plan_r2r_1d(invsinzfft,rows,zfft,zfft,FFTW_RODFT01,FFTW_PATIENT)
       call dfftw_plan_r2r_1d(fwdcoszfft,rows,zfft,zfft,FFTW_REDFT10,FFTW_PATIENT)
       call dfftw_plan_r2r_1d(invcoszfft,rows,zfft,zfft,FFTW_REDFT01,FFTW_PATIENT)

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call errcheck
       end subroutine

       subroutine killcolrow()

        if(fwddftxfft.ne.0) call dfftw_destroy_plan(fwddftxfft)
        if(invdftxfft.ne.0) call dfftw_destroy_plan(invdftxfft)
        fwddftxfft = 0
        invdftxfft = 0
        if(fwdsinzfft.ne.0) call dfftw_destroy_plan(fwdsinzfft)
        if(invsinzfft.ne.0) call dfftw_destroy_plan(invsinzfft)
        fwdsinzfft = 0
        invsinzfft = 0
        if(fwdcoszfft.ne.0) call dfftw_destroy_plan(fwdcoszfft)
        if(invcoszfft.ne.0) call dfftw_destroy_plan(invcoszfft)
        fwdcoszfft = 0
        invcoszfft = 0

!        call MPI_TYPE_FREE(chtype,ierr)

        if(allocated(xfft))   deallocate(xfft)
        if(allocated(kfft))   deallocate(kfft)
        if(allocated(zfft))   deallocate(zfft)
        if(allocated(chin))   deallocate(chin)
        if(allocated(chout))  deallocate(chout)
        if(allocated(myrows)) deallocate(myrows)
        if(allocated(myimagrows)) deallocate(myimagrows)
        if(allocated(mycols)) deallocate(mycols)
        call MPI_FINALIZE(ierr)
        call errcheck
       end subroutine

       subroutine fillmyrows()
       integer :: i,j

       ! Fill the initial matrix
       ! negative numbers indicate data "outside" the matrix
       do i=1,nch
        do j=1,M
         if(((rank*nch+i).le.rows).and.(j.le.cols)) then
          myrows(j,i) = rank*nch+i+100*j
!          myrows(j,i) = cos(2*(j-.5)*4*atan(1.0)/cols)
!          myrows(j,i) = myrows(j,i)*cos((rank*nch+i-.5)*4*atan(1.0)/rows)
!          myrows(j,i) = 10*myrows(j,i)
         else
!          myrows(j,i) = 0
          myrows(j,i) = -rank*nch-i-100*j
         end if
        end do
       end do

       end subroutine

       subroutine dumpmyrows()
       integer :: i
       print *, "My rows"
        do i=1,nch
         print *, nint(myrows(:,i))
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        end subroutine

       subroutine dumpmyimagrows()
       integer :: i
       print *, "My imaginary rows"
        do i=1,nch
         print *, nint(abs(myimagrows(:,i)))
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call errcheck
        end subroutine

        subroutine row2col()
         integer :: recv,send,i
         do i=1,numtasks
          recv = modulo(rank+i,numtasks)
          send = modulo(rank-i,numtasks)
          call makechunkR2C(myrows,chout,recv)
          call MPI_ISEND(chout,chsize,MPI_DOUBLE_PRECISION,recv,recv,MPI_COMM_WORLD,req,ierr)
          call errcheck
          call MPI_RECV(chin,chsize,MPI_DOUBLE_PRECISION,send,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          call putchunkR2C(mycols,chin,send)
          call MPI_WAIT(req,mpistat,ierr)
          call errcheck
         end do
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call errcheck
        end subroutine

        subroutine col2row
         integer :: recv,send,i
         do i=1,numtasks
          recv = modulo(rank+i,numtasks)
          send = modulo(rank-i,numtasks)
          call makechunkC2R(mycols,chout,recv)
          call MPI_ISEND(chout,chsize,MPI_DOUBLE_PRECISION,recv,recv,MPI_COMM_WORLD,req,ierr)
          call errcheck
          call MPI_RECV(chin,chsize,MPI_DOUBLE_PRECISION,send,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierr)
          call errcheck
          call putchunkC2R(myrows,chin,send)
          call MPI_WAIT(req,mpistat,ierr)
          call errcheck
         end do
        end subroutine

        subroutine dumpmycols()
         integer :: i
         print *,"my columns"
         do i=1,mch
          print *, nint(mycols(:,i))
         end do
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call errcheck
        end subroutine

       subroutine makechunkR2C(A,chunk,r)
        real, intent(in), dimension(:,:) :: A
        real, intent(out), dimension(:) :: chunk
        integer, intent(in) :: r
        integer :: i,j,chidx
        chunk = reshape(A(1+r*mch:(r+1)*mch,:),[chsize])
       end subroutine

       subroutine putchunkR2C(A,chunk,s)
        real, intent(inout), dimension(:,:) :: A
        real, intent(in), dimension(:) :: chunk
        integer, intent(in) :: s
        integer :: i,j
          A(s*nch+1:s*(nch+1),1:mch) = reshape(chunk,[nch,mch],[0.0],[2,1])
       end subroutine

       subroutine makechunkC2R(A,chunk,r)
        real, intent(in), dimension(:,:) :: A
        real, intent(out), dimension(:) :: chunk
        integer, intent(in) :: r
        chunk = reshape(A(1+r*nch:(r+1)*nch,:),[chsize])
       end subroutine

       subroutine putchunkC2R(A,chunk,s)
        real, intent(inout), dimension(:,:) :: A
        real, intent(in), dimension(:) :: chunk
        integer, intent(in) :: s
        integer :: i,j
        A((s*mch+1):s*(mch+1),1:nch) = reshape(chunk,[mch,nch],[0.0],[2,1])
       end subroutine

       subroutine fwdxfft()
        integer :: i
        do i=1,nch
         if ((rank*mch+i).le.rows) then
          xfft = myrows(1:cols,i)
          call dfftw_execute(fwddftxfft)
          myimagrows(:,i) = kfft/cols
         else
          myimagrows(:,i) = 0
         end if
        end do
       end subroutine

       subroutine fwdcosz()
        integer :: i
         do i=1,mch
          if((rank*nch+i).le.cols) then
           zfft = mycols(1:rows,i)
           call dfftw_execute(fwdcoszfft)
           mycols(1:rows,i)  = zfft/(2*rows)
          end if
         end do
       end subroutine


       subroutine invcosz()
        integer :: i
         do i=1,mch
          if((rank*nch+i).le.cols) then
           zfft = mycols(1:rows,i)
           call dfftw_execute(invcoszfft)
           mycols(1:rows,i)  = zfft
          end if
         end do
       end subroutine

       subroutine errcheck()
       if (ierr .ne. MPI_SUCCESS) then
        print *,'MPI program Error. Terminating.'
        call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
       end if
       end subroutine

       end module
