#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt



def force (r):
	t=1/r;
	return (24*(2*(t)**13 -(t)**7) -0.039);
def potential(r):
	tmp=(1/r)**6
	R_LIM=2.5
	u_R_LIM=4*(1/(R_LIM**12)-1/(R_LIM**6))
	return np.where(r<R_LIM, 4*((tmp*tmp)-(tmp)) - u_R_LIM +0.039*(r-R_LIM),0)
def pot_2(r):
	tmp=(1/r)**6
	return 4*(tmp*tmp -tmp)
x =np.linspace(1,2,50000)

plt.plot(x,force(x),label='force')
plt.plot(x,pot_2(x),label='nonscalato')
plt.plot(x,potential(x),label='potential')
plt.grid()
plt.legend()
plt.show()


