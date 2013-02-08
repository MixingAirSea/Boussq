#
#
#
 function t=cdiff(ncfile,tstart,tstop)
 tmpdir = '/home/mhoecker/tmp/Cdiff';
 nc = netcdf(ncfile,'r');
 nc = netcdf(ncfile,'r');
 t = nc{'t'}(:);
 Z = nc{'Z'}(:);
 Nt = length(t)
 NZ = length(Z)
 U = nc{'Uz'}(:);
 B = nc{'Bz'}(:);
 NU = size(U)
 NB = size(B)
 ncclose(nc);
 if(nargin>1)
  idx1 = find(t>tstart,1);
  if((nargin>2)&&max(t)>tstop)
   idx2 = find(t>tstop,1)
  else
   idx2 = Nt;
  endif
  else
   idx1 = 1;
   idx2 = Nt;
 endif
 U = U(idx1:idx2,:);
 B = B(idx1:idx2,:);
 [ZZ,tt] = meshgrid(Z,t(idx1:idx2));
 dU = dothediff(U,NZ,1);
 dB = dothediff(B,NZ,0);
#
 U_t = mean(U,1);
 B_t = mean(B,1);
#
 dU_t = dothediff(U_t,NZ,1);
 dB_t = dothediff(B_t,NZ,0);
#
 figure(1)
 pcolor(tt,ZZ,dU); shading flat; colorbar 
#
 figure(2)
 subplot(2,3,1)
 plot(mean(dU,1),Z,'r.-;U;')
 title( "< d < U >_{x,y} / dz >_t")
 subplot(2,3,2)
 plot(U_t,Z,'r.-')
 title( "< U >_{x,y,t}")
 subplot(2,3,3)
 plot(dU_t,Z,'r.-')
 title( "d < U >_{x,y,t} / dz")
 subplot(2,3,4)
 plot(mean(dB,1),Z,'b.-;B;')
 title( "< d <B>_{x,y} / dz >_t")
 subplot(2,3,5)
 plot(B_t,Z,'b.-')
 title( "< B >_{x,y,t}")
 subplot(2,3,6)
 plot(dB_t,Z,'b.-')
 title( "d< B >_{x,y,t} / dz")

endfunction

function dA = dothediff(A,NZ,BC)
 N = size(A);
 if(N(2)>1)
  dA = diff(A,1,2);
  dA = (NZ/2)*([dA,zeros(N(1),1)]+[zeros(N(1),1),dA])/2;
  switch(BC)
   case { 0 }
    dA(:,1)  = (NZ/2)*(A(:,2)+1)/2;
    dA(:,NZ) = (NZ/2)*(1-A(:,NZ-1))/2;
   case { 1 }
    dA(:,1)  = (1+dA(:,2))/2;
    dA(:,NZ) = (1+dA(:,NZ-1))/2;
   endswitch
 else
  dA = diff(U_t);
  dA = (NZ/2)*([dA,0]+[0,dA])/2;
  switch(BC)
   case { 0 }
    dA(1)  = (NZ/2)*(A(2)+1)/2;
    dA(NZ) = (NZ/2)*(1-A(NZ-1))/2;
   case { 1 }
    dA(1)  = (1+dA(2))/2;
    dA(NZ) = (1+dA(NZ-1))/2;
   endswitch
 endif
endfunction

