import numpy as np

# multiply two quaternions

def mult_quat(x,y):
	a = x[0]
	b = x[1]
	c = x[2]
	d = x[3]
	e = y[0]
	f = y[1]
	g = y[2]
	h = y[3]

	i = a*e - b*f - c*g - d*h
	j = a*f + b*e + c*h - d*g
	k = a*g - b*h + c*e + d*f
	l = a*h + b*g - c*f + d*e

	return np.array([i,j,k,l])

# using quaternions, rotate a point r by angle theta around an axis u (given by a vector)

def rotate_quat(p, theta, u):
	u = np.array(u)*(1/np.linalg.norm(u))
	#u = np.multiply(u,1/np.linalg.norm(u))
	cth = np.cos(theta/2)
	sth = np.sin(theta/2)
	q = [cth, sth*u[0], sth*u[1], sth*u[2]]
	q_inv = [cth, -sth*u[0], -sth*u[1], -sth*u[2]]

	p_quat = [0,p[0],p[1],p[2]]
	out_quat = mult_quat(mult_quat(q,p_quat),q_inv)
	return out_quat[1:4]