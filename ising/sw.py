import random as rd
import numpy as np



N = 32

class spin:
	def __init__(self, s=None, i=None,j=None,cluster = None):
		self.s = s
		self.i = i
		self.j = j
		self.cluster = cluster




def init( m ):
	for i in range(N):
		for j in range(N):
			if (rd.random()<0.5):
				m[i][j] = 1
			else:
				m[i][j] = -1


def get_neighbours ()

def clusterize (m):
	clusters=[]





if __name__ == '__main__':
	matrix = 
	init(matrix)
