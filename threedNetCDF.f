      module threedNetCDF
       use netcdf
       implicit none
       logical, private :: ncexists, ncopen
       integer, parameter :: dims=4,vars=8,kdims=5
       integer, parameter :: xd=3,yd=2,zd=1,td=4
!       integer, parameter :: kd=1,ld=2,md=3,nd=3,rd=4,tdk=5
       integer, parameter :: xv=3,yv=2,zv=1,tv=4
       integer, parameter :: uv=5,vv=6,wv=7,bv=8
!       integer, parameter :: kv=9,lv=10,mv=11,nv=12
!       integer, parameter :: fuv=13,fvv=14,fwv=15,fbv=16
       contains

       subroutine newfile(title,Re,Ri,Pr,invRo,Nx,Ny,Nz,Lz)
        character(*), intent(in) :: title
        integer, intent(in) :: Nx,Ny,Nz
        real, intent(in) :: Re,Ri,Pr,Lz
        real, intent(in) :: invRo
        integer :: ncid,ncstat
        integer :: Reid,Riid,Prid,Lzid
        integer :: Roid
        integer, dimension(dims) :: dimids
        integer, dimension(vars) :: varids
!        integer, dimension(kdims) :: mdimids
!        integer, dimension(kdims) :: ndimids
!        character(16) :: fulltitle
!*    Append the .nc suffix to the file name
!        fulltitle = trim(adjustl(title))//".nc"
        ncstat = NF90_CREATE(trim(adjustl(title))//".nc",cmode=NF90_CLASSIC_MODEL,ncid=ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_DIM(ncid,"Z",Nz,dimids(zd))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Z",NF90_DOUBLE,dimids(zd),varids(zv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(zv),"units","m")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(zv),"long_name","vertical co-ordinate")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(zv),"axis","z")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

       ncstat = NF90_PUT_ATT(ncid,varids(zv),"positive","up")
       if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_DIM(ncid,"Y",Ny,dimids(yd))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Y",NF90_DOUBLE,dimids(yd),varids(yv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(yv),"units","m")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(yv),"long_name","cross stream co-ordinate")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(yv),"axis","y")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_DIM(ncid,"X",Nx,dimids(xd))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"X",NF90_DOUBLE,dimids(xd),varids(xv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(xv),"units","m")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(xv),"long_name","streamwise co-ordinate")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(xv),"axis","x")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!       ncstat = NF90_DEF_DIM(ncid,"k",Nx/2+1,ndimids(kd))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        mdimids(kd)=ndimids(kd)

!        ncstat = NF90_DEF_VAR(ncid,"k",NF90_DOUBLE,ndimids(kd),varids(kv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(kv),"units","rad/m")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(kv),"long_name","streamwise wavenumber")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_DIM(ncid,"l",Ny,ndimids(ld))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        mdimids(ld)=ndimids(ld)

!        ncstat = NF90_DEF_VAR(ncid,"l",NF90_DOUBLE,ndimids(ld),varids(lv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(lv),"units","rad/m")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(lv),"long_name","cross stream wavenumber")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_DIM(ncid,"n",Nz,ndimids(nd))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"n",NF90_DOUBLE,ndimids(nd),varids(nv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(nv),"units","rad/m")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(nv),"long_name","vertical wavenumber")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_DIM(ncid,"m",Nz,mdimids(md))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"m",NF90_DOUBLE,mdimids(md),varids(mv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(mv),"units","rad/m")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(mv),"long_name","vertical wavenumber")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_DIM(ncid,"part",2,mdimids(rd))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ndimids(rd)=mdimids(rd)

        ncstat = NF90_DEF_DIM(ncid,"t",nf90_unlimited,dimids(td))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"t",NF90_DOUBLE,dimids(td),varids(tv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        mdimids(tdk) = dimids(td)
!        ndimids(tdk) = dimids(td)

        ncstat = NF90_PUT_ATT(ncid,varids(tv),"units","s")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(tv),"long_name","time")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(tv),"axis","t")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"u",NF90_DOUBLE,dimids,varids(uv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(uv),"units","m/s")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(uv),"long_name","streamwise velocity anomally")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(uv),"Formula","U=u+Z")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"V",NF90_DOUBLE,dimids,varids(vv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(vv),"units","m/s")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(vv),"long_name","cross-stream velocity")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"W",NF90_DOUBLE,dimids,varids(wv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(wv),"units","m/s")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(wv),"long_name","vertical velocity")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"b",NF90_DOUBLE,dimids,varids(bv))
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(bv),"units","m/s^2")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(bv),"long_name","Reduced buoyancy anomaly")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,varids(bv),"Formula","B=b+Z")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"Fu",NF90_DOUBLE,mdimids,varids(fuv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fuv),"units","m^2/s rad")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fuv),"long_name","streamwise velocity")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"Fv",NF90_DOUBLE,mdimids,varids(fvv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fvv),"units","m^2/s rad")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fvv),"long_name","streamwise velocity")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"Fw",NF90_DOUBLE,ndimids,varids(fwv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fwv),"units","m^2/s rad")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fwv),"long_name","vertical velocity")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_DEF_VAR(ncid,"Fb",NF90_DOUBLE,ndimids,varids(fbv))
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fbv),"units","m^2/s^2 rad")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_ATT(ncid,varids(fbv),"long_name","Reduced buoyancy anomaly")
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Re",NF90_DOUBLE,Reid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Reid,"long_name","Reynolds number")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Reid,"units","")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Pr",NF90_DOUBLE,Prid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Prid,"long_name","Prandtl number")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Prid,"units","")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Ri",NF90_DOUBLE,Riid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Riid,"long_name","Richardson number")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Riid,"units","")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Ro",NF90_DOUBLE,Roid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Roid,"long_name","inverse Rossby number")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Roid,"units","")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_DEF_VAR(ncid,"Lz",NF90_DOUBLE,Lzid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Lzid,"units","m")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,Lzid,"long_name","vertical dimension")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_ATT(ncid,NF90_GLOBAL,"institution","Oregon State University")
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_ENDDEF(ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_VAR(ncid,Reid,Re)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_VAR(ncid,Prid,Pr)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_VAR(ncid,Riid,Ri)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_VAR(ncid,Roid,invRo)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_PUT_VAR(ncid,Lzid,Lz)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_CLOSE(ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
        ncopen = .false.
        ncexists = .true.
       end subroutine

       subroutine cordsave(title,X,Y,Z,Nxp,ax)!k,l,m,n,Nf
        character(*), intent(in) :: title
        real, intent(in), dimension(:) :: X,Y,Z!,k,l,m,o
        integer, intent(in) :: Nxp,ax!,Nf
        integer, dimension(dims) :: dimids
        integer, dimension(vars) :: varids
        integer :: ncid,ncstat,nvar
        ncstat = NF90_OPEN(trim(adjustl(title))//".nc",NF90_WRITE,ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
        ncopen = .true.

        ncstat = NF90_INQ_VARIDS(ncid,nvar,varids)

        ncstat = NF90_PUT_VAR(ncid,varids(zv),Z)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), " saving Z", Z

        ncstat = NF90_PUT_VAR(ncid,varids(yv),Y)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), " saving Y",Y

        ncstat = NF90_PUT_VAR(ncid,varids(xv),X,start=[Nxp])
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), " saving X", Nxp, size(X),X

!        ncstat = NF90_PUT_VAR(ncid,kvarids(kv),k,start=[N])
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
!
!        ncstat = NF90_PUT_VAR(ncid,varids(lv),l)
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_VAR(ncid,kvarids(mv),m)
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

!        ncstat = NF90_PUT_VAR(ncid,kvarids(nv),n)
!        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_CLOSE(ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
        ncopen = .false.

       end subroutine

       subroutine stepsave(title,U,V,W,b,t,Nx)!fu,fv,fw,fb,Nf
        character(*), intent(in) :: title
        real,dimension(:,:,:), intent(in) :: U,V,W,b
        real, intent(in) :: t
        integer, intent(in) :: Nx!,Nf
        integer, dimension(dims) :: dimids, start
        integer, dimension(vars) :: varids
        integer :: ncid,ncstat,Nt
        ncstat = NF90_OPEN(trim(adjustl(title))//".nc",NF90_WRITE,ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
        ncopen = .true.

        ncstat = NF90_INQ_VARIDS(ncid,Nt,varids)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))

        ncstat = NF90_INQUIRE_DIMENSION(ncid,varids(tv),len=Nt)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))


        if(Nx.eq.1) then
         print *, t
         Nt = Nt+1
         ncstat = NF90_PUT_VAR(ncid,varids(tv),t,start=[Nt])
         if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), "saving t"
        end if

        start=[1,1,Nx,Nt]

        ncstat = NF90_PUT_VAR(ncid,varids(uv),U,start=start)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), "saving U",Nx

        ncstat = NF90_PUT_VAR(ncid,varids(vv),V,start=start)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), "saving V",Nx

        ncstat = NF90_PUT_VAR(ncid,varids(wv),W,start=start)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), "saving W",Nx

        ncstat = NF90_PUT_VAR(ncid,varids(bv),b,start=start)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat)), "saving b",Nx

        ncstat = NF90_CLOSE(ncid)
        if(ncstat.ne.nf90_noerr) print *, trim(nf90_strerror(ncstat))
        ncopen = .false.

       end subroutine
      end module
