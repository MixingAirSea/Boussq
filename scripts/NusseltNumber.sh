#!/bin/bash
TMPDIR="/home/mhoecker/tmp"
if [ -f $2 ] 
 then
  mv $2 $2.old
fi
for (( i = 0; i <= 1201; i++ ))
 do
  ncks -d t,$i,$i -v u,W,b,Re $1 $TMPDIR/last.nc
  PADi=`printf "%04d" $i`
  echo $PADi
  ncap2 -t 4 -v -s 'uz=u.avg($X,$Y);Uz=uz+Z;bz=b.avg($X,$Y);Bz=bz+Z;Wb=Re*(W*(b-bz)).avg($X,$Y);Wu=Re*(W*(u-uz)).avg($X,$Y);Wbz=Wb.avg($Z);Wuz=Wu.avg($Z)' $TMPDIR/last.nc $TMPDIR/Nuss$PADi.nc
  rm $TMPDIR/last.nc
 done
ncrcat $TMPDIR/Nuss*.nc $2
rm $TMPDIR/Nuss*.nc
