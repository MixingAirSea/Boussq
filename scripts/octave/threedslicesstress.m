# usage t = threedslicesstress(ncfile,loc,points)
#
# Returns the times ploted
#
# points = (integer) maximum number of points to plot in any dimension
#
# yzslices = (integer) number of slices in Y-Z to plot
#
# ncfile = (string) name of the netCDF file to be ploted;
#
# tstart = (integer) Starting time index to plot;
#
# tstop = (integer) Last time index to plot,
#         zero implies to plot the entire time series;
#
#
 function t=threedslicesstress(ncfile,loc,points,wait)
 tmpdir = '/home/mhoecker/tmp/';
 nctmp = [tmpdir"out.nc"];
 pltfile = [tmpdir 'color.plt'];
 field = ['b','u','W','V'];
 Ntold = 0;
 unix(["touch " loc "0.png"]);
 unix(["rm " loc "*.png"]);
 unix(["touch " loc "0.gif"]);
 unix(["rm " loc "*.gif"]);
 while(true)
 nc = netcdf(ncfile,'r');
 t = nc{'t'}(:);
 Nt = length(t);
 if(Nt==Ntold);
  ncclose(nc);
 else
  tstart=Ntold;
  tstop=Nt-1;
  Ntold=Nt;
  X = nc{'X'}(:);
  Y = nc{'Y'}(:);
  Z = nc{'Z'}(:);
  Lz = nc{'Lz'}(:);
  Re = nc{'Re'}(:);
  Ri = nc{'Ri'}(:);
  Pr = nc{'Pr'}(:);
  Nx = length(X);
  Ny = length(Y);
  Nz = length(Z);
  Nt = length(t);

  Xmax = max(X);
  Ymax = max(Y);
  Lx = 2*(Nx+1)*Xmax/Nx;
  Ly = 2*(Ny+1)*Ymax/Ny;


 idxfile = [loc 'index.html'];
 idxid = fopen(idxfile,'w');
 fprintf(idxid,'<TITLE>Simulation parameters</TITLE>\n');
 fprintf(idxid,'<H1>Simulation Parameters</H1>\n');
 fprintf(idxid,'Re = %g<br>\n',Re);
 fprintf(idxid,'Ri = %g<br>\n',Ri);
 fprintf(idxid,'Pr = %g<br>\n',Pr);
 fprintf(idxid,'Nx = %i<br>\n',Nx);
 fprintf(idxid,'Ny = %i<br>\n',Ny);
 fprintf(idxid,'Nz = %i<br>\n',Nz);
 fprintf(idxid,'Nt = %i<br>\n',Nt);
 fprintf(idxid,'Lx = %g<br>\n',Lx);
 fprintf(idxid,'Ly = %g<br>\n',Ly);
 fprintf(idxid,'Lz = %g<br>\n',Lz);
 fprintf(idxid,'T = %g to %g<br>\n',t(1),t(Nt));
 fprintf(idxid,'<H1>Fields</H1>\n');
 for fi=1:4
 fprintf(idxid,'<a href="%slist.html">%s</a><br>\n',field(fi),field(fi));
 endfor
 fclose(idxid);



  yzslices = floor(Nx*1.0/Ny);
  xyslices = floor(Nz*1.0/Ny);
  h = [0,1,2]*Nz;
  w = [0,1,2]*Nx;
  AR = w(3)/h(3);
  h = h./(2*Nz);
  w = w./(2*Nx);
  scale = .875;
  h = scale*h;
  w = scale*w;
  wpix = 1536;
  hpix = 1024;
  if(AR>1)
   hpix = floor(wpix/AR);
  else
   wpix=hpix*AR;
  endif
  w = w+.25*(1-w(3));
  h = h+.5*(1-h(3));
  stride = 1;
  stridex = ceil(Nx./points);
   stridey = ceil(Ny./points);
   stridez = ceil(Nz./points);
   midx = find(X>=0,1);
   midy = find(Y>=0,1);
   midz = find(Z>=0,1);
   X = X(1:stridex:Nx);
   Y = Y(1:stridey:Ny);
   Z = Z(1:stridez:Nz);
   ncclose(nc);
  for fi=1:4
   idxfile = [loc field(fi) 'list.html'];

   ncfield = ["-v " field(fi) " "];
   stream = ["-d X,,," int2str(stridex) " "];
   cross  = ["-d Y,,," int2str(stridey) " "];
   vert   = ["-d Z,,," int2str(stridez) " "];

   idxid = fopen(idxfile,'w');
   fprintf(idxid,'<TITLE>Plots of the %s field</TITLE>\n',field(fi));
   fprintf(idxid,'<H1>Simulation Parameters</H1>\n');
   fprintf(idxid,'Re = %g<br>\n',Re);
   fprintf(idxid,'Ri = %g<br>\n',Ri);
   fprintf(idxid,'Pr = %g<br>\n',Pr);
   fprintf(idxid,'Nx = %i<br>\n',Nx);
   fprintf(idxid,'Ny = %i<br>\n',Ny);
   fprintf(idxid,'Nz = %i<br>\n',Nz);
   fprintf(idxid,'Nt = %i<br>\n',Nt);
   fprintf(idxid,'Lx = %g<br>\n',Lx);
   fprintf(idxid,'Ly = %g<br>\n',Ly);
   fprintf(idxid,'Lz = %g<br>\n',Lz);
   fprintf(idxid,'T = %g to %g<br>\n',t(1),t(Nt));
   fprintf(idxid,'<H1>Snapshots</H1>\n');
   fclose(idxid);

   for tstep=tstart:tstop
    time = ["-d t," int2str(tstep) " "];

    %
    % Make X-Z plots
    %
    cblims = [0,0];
    Ytics = [];
    for i=1:2
     yidx = ["-d Y," int2str(floor((i-.5)*(Ny/2))) " "];
     unix(["/home/mhoecker/local/bin/ncks" " " time stream vert yidx ncfield ncfile " " nctmp]);
     nc = netcdf(nctmp,'r');
     Ytics = [Ytics,nc{'Y'}(:)];
     b = nc{field(fi)}(:);
     ncclose(nc);
     unix(["rm " nctmp]);

     xzdatfile = [tmpdir 'colorxz' int2str(i) '.dat'];
     datid = fopen(xzdatfile,'w');
     fwrite(datid,length(X),'float');
     fwrite(datid,X,'float');
     for  i=1:length(Z)
      fwrite(datid,Z(i),'float');
      if(field(fi)=='b')
       b(1,:,1,i)=sign(Ri)*(b(1,:,1,i)+Z(i));
      elseif(field(fi)=='u')
       b(1,:,1,i)=b(1,:,1,i)+Z(i);
      endif
      fwrite(datid,b(1,:,1,i),'float');
     endfor
    fclose(datid);

    if(field(fi)!='b')
     cblims(2)=max(cblims(2),max(max(abs(b))));
    endif
   endfor


   %
   %Make X,Y plots
   %
   "X-Y plots"
   Ztics = [];
   for i=1:xyslices
    zidx = ["-d Z," int2str(floor((i-.5)*(Nz/xyslices))) " "];
    unix(["/home/mhoecker/local/bin/ncks" " " time stream cross zidx ncfield ncfile " " nctmp]);
    nc = netcdf(nctmp,'r');
    Ztics = [Ztics,nc{'Z'}(:)];
    b = nc{field(fi)}(:);
    ncclose(nc);
    unix(["rm " nctmp]);
    xydatfile = [tmpdir 'colorxy' int2str(i) '.dat'];
    datid = fopen(xydatfile,'w');
    fwrite(datid,length(X),'float');
    fwrite(datid,X,'float');
    for  j=1:length(Y)
     fwrite(datid,Y(j),'float');
      if(field(fi)=='b')
       b(1,:,j,1)=sign(Ri)*(b(1,:,j,1)+Ztics(i));
      elseif(field(fi)=='u')
       b(1,:,j,1)=b(1,:,j,1)+Ztics(i);
      endif
      fwrite(datid,b(1,:,j,1),'float');
     endfor
    fclose(datid);

    if(field(fi)!='b')
     cblims(2)=max(max(max(abs(b))),cblims(2));
    endif
   endfor

   %
   % Make the Y-Z plots
   %
  Xtics = [];
  for i=1:yzslices
   yzdatfile = [tmpdir 'coloryz' int2str(i) '.dat'];
   xidx = ["-d X," int2str(floor((i-.5)*(Nx/yzslices))) " "];
   unix(["/home/mhoecker/local/bin/ncks" " " time xidx cross vert ncfield ncfile " " nctmp]);
   nc = netcdf(nctmp,'r');
   b = nc{field(fi)}(:);
   Xtics = [Xtics,nc{'X'}(:)];
   ncclose(nc);
   unix(["rm " nctmp]);
   datid = fopen(yzdatfile,'w');
   fwrite(datid,length(Y),'float');
   fwrite(datid,Y,'float');
   for  i=1:length(Z)
    fwrite(datid,Z(i),'float');
    if(field(fi)=='b')
     b(1,1,:,i)=sign(Ri)*(b(1,1,:,i)+Z(i));
    elseif(field(fi)=='u')
     b(1,1,:,i)=b(1,1,:,i)+Z(i);
    endif
    fwrite(datid,squeeze(b(1,1,:,i)),'float');
   endfor
   fclose(datid);

  if(field(fi)!='b')
   cblims(2)=max(max(max(abs(b))),cblims(2));
  endif
  endfor
 % Set color bar limits
 if(cblims(2)<=1e-16|field(fi)=='b')
  cblims(1)=-1;
  cblims(2)=+1;
 else
  cblims(1)=-cblims(2);
 endif

 %
 %Make the gnuplot script
 %
 cblimsid = fopen(pltfile,'w');
 idxid = fopen(idxfile,'a');
 fprintf(cblimsid,'set cbrange [%g:%g]\n',cblims(1),cblims(2));
 pngname = int2str(tstep);
 if(tstep<10)
  pngname = ["0" pngname];
 endif
 if(tstep<100)
  pngname = ["0" pngname];
 endif
 if(tstep<1000)
  pngname = ["0" pngname];
 endif
  pngname = [field(fi) pngname];
 fprintf(idxid,'<a href=%s.png>t = %g</a><br>\n',pngname,t(tstep+1));
 fprintf(cblimsid,'name = "%s"."%s".".png"\n',loc,pngname);
 fprintf(cblimsid,'set output name\n');

 fprintf(cblimsid,'set view map\n');
 fprintf(cblimsid,'set border lw 1\n');

 fprintf(cblimsid,'set pm3d\n');
 fprintf(cblimsid,'set pm3d corners2color c1\n');
%gnuplot> saw1(x,N) = N*x-floor(N*x)
%gnuplot> saw2(x,N) = 1-(N*x-floor(N*x))

 fprintf(cblimsid,'set palette mode HSV\n');
 fprintf(cblimsid,'h(x)=.8*(.5-.5*sgn(x-.5)*(1-(1-2*x+sgn(x-.5))**2)**.25)\n');
% fprintf(cblimsid,'h(x)=.8*x\n');
 fprintf(cblimsid,'s(x)=.75*(1-sgn(x-.5))+.5*(1+sgn(x-.5))*2*(2-2*x)\n');
 fprintf(cblimsid,'v(x)=s(1-x)\n');
 fprintf(cblimsid,'saw2(x,N) = .5+.5*sqrt(N*x-floor(N*x))\n');
 fprintf(cblimsid,'saw1(x,N) = .5+.5*sqrt(1-(N*x-floor(N*x)))\n');
 fprintf(cblimsid,'set palette function h(gray),s(gray)*saw1(gray,12),v(gray)*saw2(gray,12)\n');

 fprintf(cblimsid,'unset surface\n');
 fprintf(cblimsid,'unset colorbox\n');


 fprintf(cblimsid,'set term png rounded truecolor giant lw 2 size %g,%g\n',wpix,hpix);

 fprintf(cblimsid,'set multiplot \n');

% make the X-Y plots
 fprintf(cblimsid,'set colorbox user origin screen %g,%g size screen %g,%g \n',w(3)+.25*(1-w(3)),h(1),.0625*(1-w(3)),h(3)-h(1));
 fprintf(cblimsid,'set cbtics offset -1,0\n');
 fprintf(cblimsid,'set format cb "%s"\n',["% g"]);
 fprintf(cblimsid,'set label 1 "X" at screen %g,%g offset -.5,.5 \n',w(2),h(3));
 fprintf(cblimsid,'set label 2 "Z" at screen %g,%g offset -1.5,0 \n',w(1),h(2));
 fprintf(cblimsid,'set label 3 "Y" at screen %g,%g offset 0,-.5 \n',.5*(w(1)+w(2)),h(1));
 fprintf(cblimsid,'set label 4 "X" at screen %g,%g offset 0,-.5 \n',.5*(w(2)+w(3)),h(1));
 fprintf(cblimsid,'set label 5 "Y" at screen %g,%g offset .5,0 \n',w(3),(h(1)+h(2))/2);
 fprintf(cblimsid,'set label 6 "Z" at screen %g,%g offset .5,0 \n',w(3),(h(2)+h(3))/2);
 for i=1:xyslices
  xydatfile = [tmpdir 'colorxy' int2str(i) '.dat'];
  fprintf(cblimsid,'set label "%s at t=%3.1f" right at screen 1, screen %g offset -1,1 \n',field(fi),t(tstep+1),h(3));
  fprintf(cblimsid,'set lmargin at screen %g\n',w(2));
  fprintf(cblimsid,'set rmargin at screen %g\n',w(3));
  fprintf(cblimsid,'set bmargin at screen %g\n',h(1)+(h(2)-h(1))*(+(i-1)/xyslices));
  fprintf(cblimsid,'set tmargin at screen %g\n',h(1)+(h(2)-h(1))*(+(i+0.0)/xyslices));
  fprintf(cblimsid,'set xlabel " "\n');
  fprintf(cblimsid,'set ylabel " "\n');
  fprintf(cblimsid,'set yrange [%g:%g]\n',-Ymax,Ymax);
  fprintf(cblimsid,'set xrange [%g:%g]\n',-Xmax,Xmax);
  fprintf(cblimsid,'set xtics %g,%g,%g format " "\n',-Lx/2.+Lx/(2.*yzslices),Lx/(1.*yzslices),Lx/2.);
  fprintf(cblimsid,'set ytics %g,%g,%g format " "\n',-Ly/4.,Ly/2.,Ly/4.);
  fprintf(cblimsid,'splot "%s" binary matrix notitle\n',xydatfile);
  if(1==1)
   fprintf(cblimsid,'unset colorbox\n');
   fprintf(cblimsid,'unset label\n');
  endif
 endfor

 fprintf(cblimsid,'unset xlabel\n');
 fprintf(cblimsid,'unset xtics\n');
 fprintf(cblimsid,'unset ylabel\n');
 fprintf(cblimsid,'unset colorbox\n');
 fprintf(cblimsid,'unset label\n')

 % Make X-Z plots
 for i=1:2
  xzdatfile = [tmpdir 'colorxz' int2str(i) '.dat'];
 fprintf(cblimsid,'set lmargin at screen %g\n',w(i));
 fprintf(cblimsid,'set rmargin at screen %g\n',w(i+1));
 fprintf(cblimsid,'set bmargin at screen %g\n',h(2));
 fprintf(cblimsid,'set tmargin at screen %g\n',h(3));
 fprintf(cblimsid,'set xrange [%g:%g] noreverse\n',-Xmax,Xmax);
 fprintf(cblimsid,'set yrange [%g:%g] noreverse\n',-Lz/2,Lz/2);
 fprintf(cblimsid,'set xtics %g,%g,%g format " "\n',-Lx/2.+Lx/(2.*yzslices),Lx/(1.*yzslices),Lx/2.);
 fprintf(cblimsid,'set xtics in mirror\n');
 fprintf(cblimsid,'set ylabel ""\n');
 fprintf(cblimsid,'set ytics %g,%g,%g format " "\n',-Lz/2.+Lz/(2.*xyslices),Lz/(1.*xyslices),Lz/2.);
 fprintf(cblimsid,'splot "%s" binary matrix notitle\n',xzdatfile);
 endfor

 fprintf(cblimsid,'unset xtics\n');
 fprintf(cblimsid,'unset ytics\n');
 fprintf(cblimsid,'unset ylabel\n');

 fprintf(cblimsid,'set tmargin at screen %g\n',h(2));
 fprintf(cblimsid,'set bmargin at screen %g\n',h(1));
 fprintf(cblimsid,'set xrange [-%g:%g]\n',Ymax,Ymax);
 fprintf(cblimsid,'set xtics format " "\n');
 fprintf(cblimsid,'set xtics %g,%g,%g format " "\n',-Ly/4,Ly/2,Ly/4);
 fprintf(cblimsid,'set ylabel " "\n');
 fprintf(cblimsid,'set yrange  [%g:%g] noreverse\n',-Lz/2,Lz/2);
 fprintf(cblimsid,'set ytics %g,%g,%g format " "\n',-Lz/2.+Lz/(2.*xyslices),Lz/(1.*xyslices),Lz/2.);
 for i=1:yzslices
  yzdatfile = [tmpdir 'coloryz' int2str(i) '.dat'];
  fprintf(cblimsid,'set lmargin at screen %g\n',w(1)+(w(2)-w(1))*(-Ny*.5/Nx+(i-.5)/yzslices));
  fprintf(cblimsid,'set rmargin at screen %g\n',w(1)+(w(2)-w(1))*(+Ny*.5/Nx+(i-.5)/yzslices));
  fprintf(cblimsid,'splot "%s" binary matrix notitle\n',yzdatfile);
  if(i==1)
  fprintf(cblimsid,'set ytics format " "\n');
  fprintf(cblimsid,'unset ylabel\n');
  endif
 endfor
 fprintf(cblimsid,'unset xtics\n');
 fprintf(cblimsid,'unset xlabel\n');
 fprintf(cblimsid,'unset ytics\n');
 fprintf(cblimsid,'unset multiplot\n');
 fclose(cblimsid);
 fclose(idxid);
 unix(["gnuplot " pltfile]);
 unix(["rm " pltfile]);

      for i=1:2
       xzdatfile = [tmpdir 'colorxz' int2str(i) '.dat'];
       unix(["rm " xzdatfile]);
      endfor
      for i=1:xyslices
       xydatfile = [tmpdir 'colorxy' int2str(i) '.dat'];
       unix(["rm " xydatfile]);
      endfor
      for i=1:yzslices
       yzdatfile = [tmpdir 'coloryz' int2str(i) '.dat'];
       unix(["rm " yzdatfile]);
      endfor
     endfor
   endfor
   if(wait!=0)
    for fi=1:4
     idxfile = [loc field(fi) 'list.html'];
     idxid = fopen(idxfile,'a');
     typename = 'ogv';
     unix(['/home/mhoecker/bin/pngmovie.sh  -v 2048 -l ' loc field(fi) ' -n ' loc field(fi) ' -t ' typename]);
     fprintf(idxid,'<IMG SRC="%s">\n',[field(fi) '.' typename]);
     fclose(idxid);
    endfor
   endif
  endif
  pause(max([wait,1]))
 endwhile
endfunction
