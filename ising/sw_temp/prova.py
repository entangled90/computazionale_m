import random as rd
import numpy as np



N = 32

def init( m ,s):
	for i in range(N):
		for j in range(N):
			if (rd.random()<0.5):
				m[i][j] = 1
			else:
				m[i][j] = -1
			s[i][j]=-1


def get_neighbours (m,s, i,j):
	c = set(range(-1,1))
	n = set()
	for x in c:
		for y in c:
			xx = 
			yy = (j+y+N)%N
			if s[(i+x+N)%N][yy]== -1:
				n.add(s[xx][yy])
	return n


def clusterize (m):
	clusters=[]





if __name__ == '__main__':
	matrix_spin = np.empty((1,1), dtype=int,order='C')
	matrix_cluster = np.empty((1,1),dtype=int, order='C')
	init(matrix_spin,matrix_cluster)
