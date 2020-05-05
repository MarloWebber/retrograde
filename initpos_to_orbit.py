import numpy as np
import quat

class Orbit():
	# a 		Semi Major Axis
	# e 		Eccentricity
	# i 		Inclination
	# Omega 	Longitude of Ascending Node
	# omega 	Argument of Periapsis
	# nu 		True Anomaly
	def __init__(self,a,e,i,Omega,omega,nu):
		self.a = a
		self.e = e
		self.i = i
		self.Omega = Omega
		self.omega = omega
		self.nu = nu

# def initpos_to_orbit(r_vec,v_vec,mu):
# 	#r_vec is a vector giving cartesian coordinates of orbiting object,
# 	#v_vec is velocity
# 	#mu is gravitational parameter of object being orbited (=G * mass of object)

# 	#angular momentum vector, points perpendicular to plane of orbit
# 	h_vec = np.cross(r_vec, v_vec)
# 	#node vector, points in direction of ascending node (= h_vec cross z unit vector)
# 	n_vec = [-h_vec[1],h_vec[0],0]
# 	n = np.linalg.norm(n_vec)
# 	#velocity squared
# 	v_sqr = v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2
# 	#r scalar
# 	r = np.linalg.norm(r_vec)
# 	#eccentricity vector, points from apoapsis to periapsis with magnitude of the scalar eccentricity
# 	e_vec = np.multiply(v_sqr/mu - 1/r, r_vec) - np.multiply(np.dot(r_vec, v_vec)/mu, v_vec)
# 	#eccentricity scalar
# 	e = np.linalg.norm(e_vec)
# 	#specific mechanical energy
# 	sme = 0.5*v_sqr - mu/r

# 	#find semi major axis a and semi latus rectum p, use which one?

# 	#eps???
# 	if e == 1:
# 		#parabolic orbit
# 		a = np.inf
# 		p = (h_vec[0]**2 + h_vec[1]**2 + h_vec[2]**2)/mu

# 	else:
# 		#hyperbolic, elliptical or circular orbits
# 		a = -mu/(2*sme)
# 		p = a*(1 - e**2)

# 	#inclination, ranges from 0 (flat prograde) to pi (flat retrograde) 
# 	i = np.arctan2(np.sqrt(h_vec[0]**2 + h_vec[1]**2), h_vec[2])

# 	#longitude of ascending node, ranges from 0 to 2*pi
# 	if e_vec[2] == 0: #equatorial orbit, set to 0 for convenience
# 		Omega = 0
# 	else:
# 		Omega = np.arctan2(h_vec[0],h_vec[1]) % (2*np.pi)

# 	#argument of periapsis
# 	if e_vec[2] == 0: #equatorial orbit
# 		omega = np.sign(h_vec[2]) * np.arctan2(e_vec[1], e_vec[0])
# 	else:
# 		omega = np.sign(e_vec[2]) * np.arccos(np.dot(n_vec,e_vec) / (n*e))
# 	omega = omega % (2*np.pi)

# 	#true anomaly
# 	if e == 0: #circular orbit
# 		if e_vec[2] == 0: #equatorial orbit
# 			nu = -np.sign(v_vec[0]) * np.arccos(r_vec[0] / r)
# 		else:
# 			nu = np.sign(r_vec[2]) * np.arccos(np.dot(n_vec,r_vec) / (n*r))
# 	else:
# 		nu = np.sign(np.dot(r_vec,v_vec)) * np.arccos(np.dot(e_vec,r_vec) / (e*r))
# 	nu = nu % (2*np.pi)

# 	return Orbit(a,e,i,Omega,omega,nu)

# r_vec1 = np.array([0.004,0.0005])
# v_vec1 = np.array([0.005, 0.004])
# mu1 = 0.001
# print initpos_to_orbit(r_vec1,v_vec1,mu1)




#----------------

# import numpy as np
# import quat

def initpos_to_orbit(r_vec,v_vec,mu):
	#r_vec is a vector giving cartesian coordinates of orbiting object,
	#v_vec is velocity
	#mu is gravitational parameter of object being orbited (=G * mass of object)

	#angular momentum vector, points perpendicular to plane of orbit
	h_vec = np.cross(r_vec, v_vec)
	#node vector, points in direction of ascending node (= h_vec cross z unit vector)
	n_vec = [-h_vec[1],h_vec[0],0]
	n = np.linalg.norm(n_vec)
	#velocity squared
	v_sqr = v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2
	#r scalar
	r = np.linalg.norm(r_vec)
	#eccentricity vector, points from apoapsis to periapsis with magnitude of the scalar eccentricity
	e_vec = np.multiply(v_sqr/mu - 1/r, r_vec) - np.multiply(np.dot(r_vec, v_vec)/mu, v_vec)
	#eccentricity scalar
	e = np.linalg.norm(e_vec)
	#specific mechanical energy
	sme = 0.5*v_sqr - mu/r

	#find semi major axis a and semi latus rectum p, use which one?

	#eps???
	if e == 1:
		#parabolic orbit
		a = np.inf
		p = (h_vec[0]**2 + h_vec[1]**2 + h_vec[2]**2)/mu

	else:
		#hyperbolic, elliptical or circular orbits
		a = -mu/(2*sme)
		p = a*(1 - e**2)

	#inclination, ranges from 0 (flat prograde) to pi (flat retrograde) 
	i = np.arctan2(np.sqrt(h_vec[0]**2 + h_vec[1]**2), h_vec[2])

	#longitude of ascending node, ranges from 0 to 2*pi
	if e_vec[2] == 0: #equatorial orbit, set to 0 for convenience
		Omega = 0
	else:
		Omega = np.arctan2(h_vec[0],h_vec[1]) % (2*np.pi)

	#argument of periapsis
	if e_vec[2] == 0: #equatorial orbit
		omega = np.sign(h_vec[2]) * np.arctan2(e_vec[1], e_vec[0])
	else:
		omega = np.sign(e_vec[2]) * np.arccos(np.dot(n_vec,e_vec) / (n*e))
	omega = omega % (2*np.pi)

	#true anomaly
	if e == 0: #circular orbit
		if e_vec[2] == 0: #equatorial orbit
			nu = -np.sign(v_vec[0]) * np.arccos(r_vec[0] / r)
		else:
			nu = np.sign(r_vec[2]) * np.arccos(np.dot(n_vec,r_vec) / (n*r))
	else:
		nu = np.sign(np.dot(r_vec,v_vec)) * np.arccos(np.dot(e_vec,r_vec) / (e*r))
	nu = nu % (2*np.pi)

	#eccentric anomaly
	E = np.arctan2(np.sqrt(1 - e**2) * np.sin(nu), e + np.cos(nu)) % (2*np.pi)
	#mean anomaly
	M = E - e * np.sin(E)

	# return np.array([a,e,i,Omega,omega,M])
	return Orbit(a,e,i,Omega,omega,M)

def orbit_to_position(a,e,i,Omega,omega,M_0,t_0,t,mu):
	#mean anomaly
	M = M_0 + (t - t_0)*np.sqrt(mu/a**3)

	#solve keplers equation using newtons method to find eccentric anomaly
	E_0 = M
	tol = 1e-10
	n = 0
	while n < 1000:
		E_1 = E_0 - (E_0 - e*np.sin(E_0) - M)/(1 - e*np.cos(E_0))
		#print(n, E_0, E_1)
		if np.abs(E_1 - E_0) < tol:
			n = 1000
			E = E_1
		else:
			E_0 = E_1
			n += 1


	#true anomaly
	nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
	#distance to central body
	r = a*(1 - e*np.cos(E))

	#position in orbital plane
	o_vec = np.multiply(r, [np.cos(nu), np.sin(nu), 0])

	#velocity in orbital plane
	ov_vec = np.multiply(np.sqrt(mu*a)/r, [-np.sin(E), np.sqrt(1-e**2)*np.cos(E), 0])

	#sin and cos of each angle to be rotated by
	# somega = np.sin(omega)
	# comega = np.cos(omega)
	# sOmega = np.sin(Omega)
	# cOmega = np.cos(Omega)
	# si = np.sin(i)
	# ci = np.cos(i)

	#rotation matrix
	#rotate = np.array([
	#[comega*cOmega - somega*ci*sOmega, -somega*cOmega - comega*ci*sOmega, 0],
	#[comega*sOmega + somega*ci*cOmega, -somega*sOmega + comega*ci*cOmega, 0],
	#[somega*si, comega*si, 0]])

	#r_vec = np.matmul(rotate, o_vec)
	#v_vec = np.matmul(rotate, ov_vec)

	r_vec = quat.rotate_quat(quat.rotate_quat(quat.rotate_quat(o_vec,-omega,[0,0,1]),-i,[1,0,0]),-Omega,[0,0,1])
	v_vec = quat.rotate_quat(quat.rotate_quat(quat.rotate_quat(ov_vec,-omega,[0,0,1]),-i,[1,0,0]),-Omega,[0,0,1])
	
	y = np.array([r_vec, v_vec])
	return np.reshape(y,6)

# r1 = np.array([1,0,0])
# v1 = np.array([0,1,0])
# mu = 100
# [a,e,i,Omega,omega,M] = initpos_to_orbit(r1,v1,mu)
# print([a,e,i,Omega,omega,M])

# y = orbit_to_position(a,e,i,Omega,omega,M,0,0,mu)
# print(y)

def initpos_to_orbit_2d(r_vec,v_vec,mu):
	r_vec_3 = np.array([r_vec[0], r_vec[1], 0])
	v_vec_3 = np.array([v_vec[0], v_vec[1], 0])

	return initpos_to_orbit(r_vec_3,v_vec_3,mu)