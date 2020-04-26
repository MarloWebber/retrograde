import numpy as np

def func_initpos_to_orbit(r_vec,v_vec,mu):
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

	return [a,e,i,Omega,omega,nu]

# r_vec1 = np.array([0.004,0.0005])
# v_vec1 = np.array([0.005, 0.004])
# mu1 = 0.001
# print initpos_to_orbit(r_vec1,v_vec1,mu1)