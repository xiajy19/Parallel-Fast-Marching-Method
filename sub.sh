#!/bin/zsh

fmm="fm3dt_parallel"
subv="subinv"

itn=subiter.in
ing=gridi.vtx
sng=gridc.vtx
reft=rtimes.dat
ttim=rtravel.out

typeset -i COUNT
typeset -i ITER

cp $ing $sng

COUNT=1
exec 4>$itn
	print -u4 $COUNT
exec 4<&-

ITER=1
echo "Running Parallel FMM..."
mpirun -np 9 ./$fmm
echo "Parallel FMM Done"
cp $ttim $reft
while [[ $ITER -le 5 ]]
do
	echo "Running Inversion..."
	./$subv
	echo "Inversion Done" 
	echo "Running Parallel FMM..."
	mpirun -np 9 ./$fmm
	echo "Parallel FMM Done"
	ITER=ITER+1
	COUNT=ITER
	exec 4>$itn
		print -u4 $COUNT
	exec 4<&-
done
