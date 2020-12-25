#!/bin/bash


#information on the system
NGEOM=1000  #total number of fragmenys
ATOMS=40

#memory 
maxcore="2000"

for (( c=1;c<$NGEOM; c++ ))
do
  index=$(($c*10))
  time=$(($index*10))
  grep -A40 -m 1 "MD iter: ${index}" geom.xyz | tail -n +2 > geom.tmp
  echo "start DSB_rt_all" >> DSB_rt_all.nw
  echo "memory 2000 mb" >> DSB_rt_all.nw 
  echo "charge 0" >> DSB_rt_all.nw
  echo "geometry noautosym noautoz" >> DSB_rt_all.nw
  tail -$ATOMS geom.tmp >> DSB_rt_all.nw
  echo "end" >> DSB_rt_all.nw
  echo "   " >> DSB_rt_all.nw
  echo "semiemp" >>DSB_rt_all.nw
  echo "mult 1.0" >> DSB_rt_all.nw
  echo "apx INDO/1" >> DSB_rt_all.nw
  echo "scftype rhf" >> DSB_rt_all.nw
  echo "maxiter 100" >> DSB_rt_all.nw
  echo "inttyp 1" >> DSB_rt_all.nw
  echo "end" >> DSB_rt_all.nw
  
  echo "rt_semiemp" >> DSB_rt_all.nw
  echo "nrestarts 10" >> DSB_rt_all.nw
  echo "tmax ${time}" >> DSB_rt_all.nw
  echo "dt 0.1"  >> DSB_rt_all.nw
  echo "norb 12" >> DSB_rt_all.nw
  echo "numorb 48 49 50 51 52 53 54 55 56 57 58 59" >> DSB_rt_all.nw
  if [[ $c -eq 1 ]] 
  then	  
  echo "occuorb 2.0 2.0 2.0 2.0 1.59 1.41 0.41 0.59 0.0 0.0 0.0 0.0" >> DSB_rt_all.nw
  fi
  echo "propagator chebyshev" >> DSB_rt_all.nw
  echo "tag "kick_y" " >> DSB_rt_all.nw
  echo "field_type delta" >> DSB_rt_all.nw
  echo "field_max 0.0000" >> DSB_rt_all.nw
  echo "polarization y" >> DSB_rt_all.nw
  echo "print moocc energy" >> DSB_rt_all.nw
  if [[ $c -ne 1 ]]
  then
  echo "load restart" >> DSB_rt_all.nw
  fi
  echo "end" >> DSB_rt_all.nw
 
  echo "task semiemp rt_semiemp" >> DSB_rt_all.nw

done 
