mkdir ./dipoleDAT
mkdir ./dipoleDAT/x_dat
mkdir ./dipoleDAT/y_dat
mkdir ./dipoleDAT/z_dat

grep "kick" *.out > output.txt
grep system output.txt | grep kick_xtotal  | grep Dipole | awk '{print $2 "   "$3 "   " $4 "   " $5}' > ./dipoleDAT/x_dat/dip.dat
grep system output.txt | grep kick_ytotal  | grep Dipole | awk '{print $2 "   "$3 "   " $4 "   " $5}' > ./dipoleDAT/y_dat/dip.dat
grep system output.txt | grep kick_ztotal  | grep Dipole | awk '{print $2 "   "$3 "   " $4 "   " $5}' > ./dipoleDAT/z_dat/dip.dat

cd ./dipoleDAT/x_dat/ 
~/RT-semiemp/rttddftARG_spec 200 0.0001
cd ../../

cd ./dipoleDAT/y_dat/ 
~/RT-semiemp/rttddftARG_spec 200 0.0001
cd ../../

cd ./dipoleDAT/z_dat/ 
~/RT-semiemp/rttddftARG_spec 200 0.0001
cd ../../


awk '{print $1}' ./dipoleDAT/x_dat/fft.dat > ./dipoleDAT/tTMP.dat
awk '{print $3}' ./dipoleDAT/x_dat/fft.dat > ./dipoleDAT/xTMP.dat
awk '{print $4}' ./dipoleDAT/y_dat/fft.dat > ./dipoleDAT/yTMP.dat
awk '{print $5}' ./dipoleDAT/z_dat/fft.dat > ./dipoleDAT/zTMP.dat

paste ./dipoleDAT/tTMP.dat ./dipoleDAT/xTMP.dat ./dipoleDAT/yTMP.dat ./dipoleDAT/zTMP.dat  > ./dipoleDAT/fftAll.dat
rm -f ./dipoleDAT/tTMP.dat
rm -f ./dipoleDAT/xTMP.dat
rm -f ./dipoleDAT/yTMP.dat
rm -f ./dipoleDAT/zTMP.dat




