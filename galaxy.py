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
		'Chikundipur',
		640000,													# radius
		1,														# density
		0.9,													# surface friction
		[200,180,120,255],										# color
		[200,190,155,255],										# outline color
		earthAtmosphere,										# atmosphere object
		[1,1],													# position
		make_clouds(650000, [1,1], 1000, 300, 30000, 640000),	# clouds
		0
		)

chikmagathur_instance = Actor('Chikmagalur', shipyard('space_station_1'),[0,-700000], [51500,0], 0)
chikmagathur_instance.dockMessage = "Chikmagalur\nThe people who live above the desert world grow coffee in these huge metal cans.\nThey are happy to trade with you and to share their delicious food.\nOn your first visit here, you met a domestic cat, the first living animal you had seen in your whole life."

jellostation_instance = Actor('Jelly Pixel', shipyard('jellyStation'), [0,700000], [-51500,0], 0)
jellostation_instance.dockMessage = "Jelly Pixel\nThe inhabitants of this beautiful station live in luxury.\nThey spend their time creating sculptures and music.\n They have a symbiotic relationship with their more industrious neighbours."
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

alpha_centauri = SolarSystem(
		"α Centauri",												
		[earth, moon,chikmagathur_instance, jellostation_instance],											
		[-180,-90],													
		[240,200,120,255],										
		[255,255,255,255],									
		["Barnard's Star", "Wolf 359", "Lalande 21185","Sirius"],									
		1000000												
		)
solarSystems.append(alpha_centauri)
# sol_III = SolarSystem(
# 		'Sol III',												# name
# 		[earth, moon,chikmagathur_instance, jellostation_instance ],											# contents
# 		[0,0],													# map position
# 		[10,50,200,255],										# map color
# 		[255,255,255,255],										# map outline color
# 		['Sol IV', 'Procyon'],									# hyperdrive links
# 		4000000,												# hyperdrive threshold
# 		)
# solarSystems.append(sol_III)

# Mars, Phobos and Deimos system
# marsAtmosphere = Atmosphere(
# 		560000,
# 		10000,
# 		0.1,
# 		0,
# 		[50,50,50,255],
# 		[0,0,0,255]
# 		)
# mars = Attractor(
# 		'Rakshasa',
# 		560000,
# 		0.8,
# 		0.9,
# 		[125,45,25,255],
# 		[155,75,55,255],
# 		marsAtmosphere,
# 		[1,1]
# 		)
# sol_IV = SolarSystem(
# 		'Sol IV',												
# 		[mars],											
# 		[10,10],													
# 		[200,50,10,255],										
# 		[255,255,255,255],									
# 		['Sol III', 'Procyon'],									
# 		3000000												
# 		)
# solarSystems.append(sol_IV)

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
		[-250,340],													
		[110,245,60,255],										
		[255,255,255,255],									
		["Sirius"],									
		2000000												
		)
solarSystems.append(procyon)



# moti_mahal = Attractor(
# 		'Moti Mahal',
# 		100000,
# 		50,
# 		0.1,
# 		[60,150,255,255],
# 		[150,200,255,255],
# 		None,
# 		[1,1]
# 		)
# mehrangarh_contents.append(moti_mahal)
# n_mahal_units = 16
# mahal_array = []
# mahal_array_height = 200000
# mahal_array_speed = 45000
# mehrangarh_contents = []
# stepSize = 2*math.pi/n_mahal_units
# for i in range(0,n_mahal_units):
# 	mahal_unit = shipyard('nano orbital section')
# 	mahal_unit_instance = Actor('nano orbital section', mahal_unit,(mahal_array_height * math.cos(stepSize* i), mahal_array_height * math.sin(stepSize* i)), [mahal_array_speed * math.cos((stepSize* i) + (0.5* math.pi) ),mahal_array_speed * math.sin((stepSize* i) + (0.5 * math.pi))], stepSize*i)
# 	mehrangarh_contents.append(mahal_unit_instance)



# mehrangarh = SolarSystem(
# 		'Noida',												
# 		mehrangarh_contents,											
# 		[1000,-200],													
# 		[0,120,255,255],										
# 		[255,255,255,255],									
# 		['Procyon'],									
# 		2000000												
# 		)
# solarSystems.append(mehrangarh)







barnards_star = SolarSystem(
		"Barnard's Star",												
		[],											
		[8,-200],													
		[110,35,0,255],										
		[255,255,255,255],									
		["α Centauri"],									
		1000000												
		)
solarSystems.append(barnards_star)



wolf_359 = SolarSystem(
		"Wolf 359",												
		[],											
		[-380,80],													
		[70,65,245,255],										
		[255,255,255,255],									
		["α Centauri"],									
		1000000												
		)
solarSystems.append(wolf_359)

lalande_21185 = SolarSystem(
		"Lalande 21185",												
		[],											
		[-405,72],													
		[10,120,50,255],										
		[255,255,255,255],									
		[],									
		1000000												
		)
solarSystems.append(lalande_21185)

sirius = SolarSystem(
		"Sirius",												
		[],											
		[-90,273],													
		[130,250,255,255],										
		[255,255,255,255],									
		["Procyon", "ε Eridani", "α Centauri"],									
		1000000												
		)
solarSystems.append(sirius)

rakshasaAtmosphere = Atmosphere(
		560000,
		10000,
		0.1,
		0,
		[50,50,50,255],
		[0,0,0,255]
		)
rakshasa = Attractor(
		'Rakshasa',
		560000,
		0.8,
		0.9,
		[125,45,25,255],
		[155,75,55,255],
		rakshasaAtmosphere,
		[1,1]
		)
epsilon_eridani = SolarSystem(
		"ε Eridani",												
		[rakshasa],											
		[315,275],													
		[245,105,65,255],										
		[255,255,255,255],									
		["Sirius", "GJ 1061", "τ Ceti"],									
		1000000												
		)
solarSystems.append(epsilon_eridani)

gj_1061_contents = []
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
gj_1061_contents.append(moti_mahal)
n_mahal_units = 16
mahal_array = []
mahal_array_height = 200000
mahal_array_speed = 45000
gj_1061_contents = []
stepSize = 2*math.pi/n_mahal_units
for i in range(0,n_mahal_units):
	mahal_unit = shipyard('nano orbital section')
	mahal_unit_instance = Actor('nano orbital section', mahal_unit,(mahal_array_height * math.cos(stepSize* i), mahal_array_height * math.sin(stepSize* i)), [mahal_array_speed * math.cos((stepSize* i) + (0.5* math.pi) ),mahal_array_speed * math.sin((stepSize* i) + (0.5 * math.pi))], stepSize*i)
	gj_1061_contents.append(mahal_unit_instance)

gj_1061 = SolarSystem(
		"GJ 1061",												
		gj_1061_contents,											
		[348,311],													
		[10,165,255,255],										
		[255,255,255,255],									
		["ε Eridani"],									
		1000000												
		)
solarSystems.append(gj_1061)

tau_ceti = SolarSystem(
		"τ Ceti",												
		[],											
		[531,168],													
		[245,155,255,255],										
		[255,255,255,255],									
		["ε Eridani"],									
		1000000												
		)
solarSystems.append(tau_ceti)
