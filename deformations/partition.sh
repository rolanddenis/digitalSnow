#!/bin/bash

S=10 # Display&check step
S_max=10000 # Maximum number of step
D=3 # Dimension
P=8 # Number of phase
N_list=(32 32 64 64 128 128 256 256 512 512 1024 1024)   # List of domain size
EPS_list=(2 1 2  1  2   1   2   1   2   1   2    1)     # List of epsilon

if [[ $D == 2 ]]
then
    EXT="pgm"
else
    EXT="vol"
fi

N_last=0
for ((i=0;i<${#N_list[@]};++i))
do
    N=${N_list[i]}
    eps=${EPS_list[i]}
    work_dir=`printf "phase_%02d_N%d_eps%d" $i $N $eps`
    rm -r $work_dir
    mkdir -p $work_dir
    cd $work_dir

    if [[ $i == 0 ]]
    then
        ../deformation${D}d -d $N -s rand --noDist --displayStep $S -e $eps -t $(($eps*$eps)) -n $S_max -p $P
    else
        if [[ $N_last == $N ]]
        then
            ../deformation${D}d -i ../last_result.${EXT} --displayStep $S -e $eps -t $(($eps*$eps)) -n $S_max
        else
            ../deformation${D}d -i ../last_result.${EXT} --overSample $(($N/$N_last)) --displayStep $S -e $eps -t $(($eps*$eps)) -n $S_max
        fi
    fi

    cp `ls -1 *.${EXT} | tail -n 1` ../last_result.${EXT}
    cd ..

    N_last=$N
done
