#/bin/zsh
zmodload -i zsh/mathfunc

ETA_MAX=0.60
ETA=0.30
STEP=0.0005
while (( $ETA < $ETA_MAX ))
do
	echo "Starting sfere3D with ETA= ${ETA} "
	time ./main2 ${ETA}
	ETA=$((ETA+STEP))
done

