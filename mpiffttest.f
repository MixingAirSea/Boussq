         program mpiffttest
         use colrow
         implicit none
         integer :: i, reps, Nx,Nz
         real :: tic,toc
         Nx = 768
         Nz = 512
         reps = 1024
         call initcolrow(Nx,Nz)
         call fillmyrows()
!         call dumpmyrows
         call cpu_time(toc)
         do i = 1,reps
          call row2col()
          call fwdcosz()
!          call dumpmycols
          call col2row()
!          call dumpmyrows
          call fwdxfft()
         end do
         call cpu_time(tic)
         print *, (tic-toc)/(reps*Nx*Nz)
!         call dumpmyimagrows
         call killcolrow()
         end program
