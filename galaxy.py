# import Atmosphere
# import Attractor
# import SolarSystem

from retrograde import *
from ships import *

solarSystems = []


# Earth and Moon system
earthAtmosphere = Atmosphere(
		640000,													# planet radius
		35000,													# atmosphere depth
		1,														# density at the bottom
		0,														# density at the top
		[70,145,220,255],										# color at the bottom
		[0,0,0,255]												# color at the top
		)
earth = Attractor(
		'Chikunipur',
		640000,													# radius
		1,														# density
		0.9,													# surface friction
		[150,150,150,255],										# color
		[200,190,155,255],										# outline color
		earthAtmosphere,										# atmosphere object
		[1,1],													# position
		make_clouds(650000, [1,1], 1000, 300, 30000, 640000),	# clouds
		0
		)
moon = Attractor(
		'Munus',
		160000,		
		1,
		0.9,
		[45,45,45,255],
		[145,145,145,255],
		None,
		[3000000,-3000000]
		)
sol_III = SolarSystem(
		'Sol III',												# name
		[earth, moon],											# contents
		[0,0],													# map position
		[10,50,200,255],										# map color
		[255,255,255,255],										# map outline color
		['Sol IV', 'Procyon'],									# hyperdrive links
		4000000,												# hyperdrive threshold
		)
solarSystems.append(sol_III)

# Mars, Phobos and Deimos system
marsAtmosphere = Atmosphere(
		560000,
		10000,
		0.1,
		0,
		[50,50,50,255],
		[0,0,0,255]
		)
mars = Attractor(
		'Rakshasa',
		560000,
		0.8,
		0.9,
		[125,45,25,255],
		[155,75,55,255],
		marsAtmosphere,
		[1,1]
		)
sol_IV = SolarSystem(
		'Sol IV',												
		[mars],											
		[10,10],													
		[200,50,10,255],										
		[255,255,255,255],									
		['Sol III', 'Procyon'],									
		3000000												
		)
solarSystems.append(sol_IV)

# Yhivi system
yhiviAtmosphere = Atmosphere(
		100000,
		20000,
		0.5,
		0,
		[210,255,120,255],
		[0,0,0,255]
		)
yhivi = Attractor(
		'Yhivi',
		100000,
		5,
		0.9,
		[130,220,120,255],
		[220,255,180,255],
		yhiviAtmosphere,
		[1,1]
		)
procyon = SolarSystem(
		'Procyon',												
		[yhivi],											
		[-20,100],													
		[10,200,50,255],										
		[255,255,255,255],									
		['Sol IV', 'Mehrangarh'],									
		2000000												
		)
solarSystems.append(procyon)


# Yhivi system
moti_mahal = Attractor(
		'Moti Mahal',
		100000,
		50,
		0.1,
		[60,150,255,255],
		[150,200,255,255],
		None,
		[1,1]
		)

n_mahal_units = 32
mahal_array = []
mahal_array_height = 200000
mahal_array_speed = 45000
mehrangarh_contents = []
stepSize = 2*math.pi/n_mahal_units
for i in range(0,n_mahal_units):
	mahal_unit = shipyard('nano orbital section')
	mahal_unit_instance = Actor('nano orbital section', mahal_unit,(mahal_array_height * math.cos(stepSize* i), mahal_array_height * math.sin(stepSize* i)), [mahal_array_speed * math.cos((stepSize* i) + (0.5* math.pi) ),mahal_array_speed * math.sin((stepSize* i) + (0.5 * math.pi))], stepSize*i)
	mehrangarh_contents.append(mahal_unit_instance)

mehrangarh_contents.append(moti_mahal)

mehrangarh = SolarSystem(
		'Mehrangarh',												
		mehrangarh_contents,											
		[1000,-200],													
		[0,120,255,255],										
		[255,255,255,255],									
		['Procyon'],									
		2000000												
		)
solarSystems.append(mehrangarh)


