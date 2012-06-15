# usage t = threedslicesstressavg(ncfile,loc,tstop,tstart)
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
 function t=threedslicesstressavg(ncfile,loc,tstop,tstart)
 tmpdir = '/home/mhoecker/tmp/avg';
 nctstep = [tmpdir"tstep.nc"];
 nctmp = [tmpdir"out.nc"];
 nctmp2 = [tmpdir"out2.nc"];
 pltfile = [tmpdir '.plt'];
 field = ['b','u','W','V'];
 FIELD = ['B','U','W','V'];
 nc = netcdf(ncfile,'r');
 t = nc{'t'}(:);
 X = {nc{'X'}(:),nc{'Y'}(:),nc{'Z'}(:)};
 Lz = nc{'Lz'}(:);
 Re = nc{'Re'}(:);
 Ri = nc{'Ri'}(:);
 Pr = nc{'Pr'}(:);
 Ro = nc{'Ro'}(:);
 ncclose(nc);
 if(Ri==0)
  Ri=1;
 endif
 if(Ro>0)
  Ro=num2str(1./Ro);
 else
  Ro = "Infinity";
 endif
 Nt = length(t);
 NX = [];
 LX = [];
 strideX = [];
 for i=1:3
  NX = [NX,size(X{i})(1)];
  LX = [LX,size(X{3})(1)*Lz./NX(i)];
 endfor
 axid = ["X","Y","Z"];
 #
 [NA,idxX] = sort(NX);
 axid = axid(idxX);
 A = {X{idxX}};
 LA = LX(idxX);
#
 axorder = [axid(1),",",axid(2),",",axid(3)];
#
 bcslices = floor(NA(1)./NA(3));
 abslices = floor(NA(2)./NA(3));
 slicenc = "nice /home/mhoecker/local/bin/ncks -O ";
 ncap2 = ["nice ncap2 -O "];
 avgnc   = ["nice ncwa -O "];
 avgncBB = [avgnc " -a " axid(2) " -b "];
 avgncAB = [avgnc " -a " axid(1) "," axid(2) " -b "];
 reordernc = ["nice /home/mhoecker/local/bin/ncpdq -O -a " axorder " "];
 Aidx = [" -d " axid(1) ",1 "];
 Bidx = [" -d " axid(2) ",1 "];
 Cidx = [" -d " axid(3) ",1 "];
 tstop = min([Nt,tstop]);
for ti=tstart:tstop
 time = [" -d t," int2str(ti) " "];
 timeslice = [slicenc time " -v "];
 for fi=1:3
  timeslice = [timeslice field(fi) ","];
 endfor
 timeslice = [timeslice field(4) " "];
 timeslice = [timeslice ncfile " " nctstep];
 unix(timeslice)
 addlaminar = [ncap2 " -s 'B=b+Z' -s 'U=u+Z' " nctstep " " nctmp];
 unix(addlaminar)
 reorder = [reordernc nctmp " " nctmp2];
 unix(reorder);
 unix(["rm " nctmp]);
 k=0;
 for i=1:3
  for j=i+1:3
   if(i==j)
    avgncij = [avgnc " -a " axid(i) " -b "];
   else
    avgncij = [avgnc " -a " axid(i) ","  axid(j) " -b "];
   endif
   avgncij = [avgncij nctmp2 " " nctmp];
   unix(avgncij);
   iax = find(([1,2,3]!=i)&([1,2,3]!=j)) ;
   nc = netcdf(nctmp,'r');
   k=k+1;
#   figure(k)
   for fi = 1:4
    b = squeeze(nc{FIELD(fi)}(:));
    if(size(b)(2)>size(b)(1))
    b = b';
    endif
    brange = max(abs(b));
    if brange==0
     brange = 1;
    endif
    Arange = max(abs(A{iax}));
#
    datfile = [tmpdir FIELD(fi) 'avg'  axid(i) axid(j) int2str(ti) '.dat'];
    datid = fopen(datfile,'w');
    fwrite(datid,[b';A{iax}'],'float');
    fclose(datid);
    datform = "%float%float";
    pltid = fopen(pltfile,'w');
    fprintf(pltid,'set term png\n');
    fprintf(pltid,'set output "%s"\n',[loc FIELD(fi) "_" axid(i) axid(j) "avg-" int2str(ti) ".png"]);
    fprintf(pltid,'set ylabel "%s"\n',axid(iax));
    fprintf(pltid,'set yrange [%g:%g]\n',-Arange,Arange);
    fprintf(pltid,'set xlabel "%s"\n',FIELD(fi));
    fprintf(pltid,'set xrange [%g:%g]\n',-brange,brange);
    fprintf(pltid,'plot "%s" binary format="%s" u 1:2 w lines title "%s"\n',datfile,datform,["<" FIELD(fi) ">_{" axid(i) axid(j) "}"]);
    fclose(pltid);
    unix(["gnuplot " pltfile]);
    unix(["rm " datfile]);
    unix(["rm " pltfile]);
#
#    subplot(2,2,fi)
#    plot(b,A{iax},[";<" FIELD(fi) ">_{" axid(i) axid(j) "};"])
#    axis([-brange,brange,-Arange,Arange]);
#    xlabel(FIELD(fi))
#    ylabel(axid(iax))
#    print([int2str(ti) "-" axid(iax) "avg.png"],"-dpng")
   endfor
   ncclose(nc);
  endfor
 endfor
 unix(["rm " nctmp]);
endfor

# abdatfile = [tmpdir 'AB' int2str(i) '.dat'];
# datid = fopen(abdatfile,'w');
# fclose(datid);

# unix(["rm " tmpdir "*"]);
endfunction
