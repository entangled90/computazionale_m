#/usr/bin/python

import os

ETA_MAX=0.76
ETA=0.4
STEP=0.1
while ( ETA < ETA_MAX )
	print("Starting sfere2D with ETA= ${ETA} ")
	
	call ["time", "./main","ETA"] 
	ETA=ETA+STEP

