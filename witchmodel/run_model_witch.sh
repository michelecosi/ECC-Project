#!/bin/bash
cd "C:/Users/modin/Desktop/Ettore/UNIVERSITA/ECC_GAMS/WITCH/witchmodel"


# This is to create a tider nameout 
name=()
for c in {r5,witch13}\_{bau,bau-impacts}  ; do
    name+=($c)
done

# this is to launch the different flags
arg=()
for c in  --n={r5,witch13}\,--policy={bau,bau-impacts} ; do
    arg+=($c)
done


## DON'T MODIFY
for ((i=0 ; i < ${#arg[@]} ; ++i)) ; do
    echo ${arg[i]} --nameout=${name[i]} 
    #C:/GAMS/42/gams run_witch.gms ${arg[i]} --nameout=${name[i]} 
done

