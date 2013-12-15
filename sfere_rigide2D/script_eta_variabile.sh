#/usr/bin/zsh
zmodload -i zsh/mathfunc

ETA_MAX=0.75
ETA=0.3
STEP=0.01
while (( $ETA<$ETA_MAX ))
do
	echo "Starting sfere2D with ETA= ${ETA} "
	time ./main ${ETA}
	ETA=$((ETA+STEP))
done

