import ..retrograde

solarSystems = []

# Earth and Moon system
earthAtmosphere = retrograde.Atmosphere(
		640000,													# planet radius
		15000,													# atmosphere depth
		1,														# density at the bottom
		0,														# density at the top
		[70,145,220,255],										# color at the bottom
		[0,0,0,255]												# color at the top
		)
earth = retrograde.Attractor(
		640000,													# radius
		1,														# density
		0.9,													# surface friction
		[180,170,145,255],										# color
		[200,190,155,255],										# outline color
		earthAtmosphere											# atmosphere object
		[1,1]													# position
		)
moon = retrograde.Attractor(
		160000,		
		1,
		0.9,
		[45,45,45,255],
		[145,145,145,255],
		None,
		[3000000,-3000000]
		)
sol_III = retrograde.SolarSystem(
		'Sol III',												# name
		[earth, moon],											# contents
		[0,0],													# map position
		[10,50,200,255],										# map color
		[255,255,255,255],										# map outline color
		['Sol IV', 'Procyon'],									# hyperdrive links
		4000000,												# hyperdrive threshold
		)
galaxy.append(sol_III)

# Mars, Phobos and Deimos system
marsAtmosphere = retrograde.Atmosphere(
		560000,
		10000,
		0.1,
		0,
		[50,50,50,255],
		[0,0,0,255],
		)
mars = retrograde.Attractor(
		560000
		0.8
		0.9
		[125,45,25,255]
		[155,75,55,255]
		marsAtmosphere
		)
sol_IV = retrograde.SolarSystem(
		'Sol IV',												
		[mars],											
		[0,0],													
		[200,50,10,255],										
		[255,255,255,255],									
		['Sol III', 'Procyon'],									
		3000000,												
		)
galaxy.append(sol_IV)

# Yhivi system
yhiviAtmosphere = retrograde.Atmosphere(
		100000,
		20000,
		0.5,
		0,
		[210,255,120,255],
		[0,0,0,255],
		)
yhivi = retrograde.Attractor(
		100000
		5
		0.9
		[130,220,120,255]
		[220,255,180,255]
		yhiviAtmosphere
		)
procyon = retrograde.SolarSystem(
		'Procyon',												
		[yhivi],											
		[0,0],													
		[10,200,50,255],										
		[255,255,255,255],									
		['Sol IV'],									
		2000000,												
		)
galaxy.append(procyon)
