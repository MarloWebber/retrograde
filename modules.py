from retrograde import *

'''
Modules are the granular parts of ships and buildings and other actors in the game.
They each consume, store, and produce resources, and the resources can be used for many things.

The module class is arranged like a factory that can produce many individual module objects, because there very many in the game world and some are of the same type as each other, but are still different individuals.
'''

class ModuleEffect(): # a ModuleEffect is just a polygon that is visually displayed when a module is active.
	def __init__(self, effectType, offset=[0,0]):
		self.offset = offset

		if effectType == 'engine 10 flame':
			self.radius = 3
			self.points = [[-self.radius, -self.radius], [self.radius, -self.radius], [0,2*self.radius]]
			self.color = [255,255,255,255]
			self.outlineColor = [255,220,150,255]
			self.illuminator = Illuminator(offset, 1500, self.outlineColor, 3)
		if effectType == 'engine 100 flame':
			self.radius = 9
			self.points = [[-self.radius, -self.radius], [self.radius, -self.radius], [0,2*self.radius]]
			self.color = [255,255,255,255]
			self.outlineColor = [255,220,150,255]
			self.illuminator = Illuminator(offset, 2500, self.outlineColor, 3)

		if effectType == 'cannon 10 flash':
			self.radius = 1
			self.points = [[-self.radius, self.radius], [self.radius, self.radius], [0,-2*self.radius]]
			self.color = [255,255,255,255]
			self.outlineColor = [245,250,255,255]
			self.illuminator = Illuminator(offset, 100, self.outlineColor, 20)

class Module():
	def __init__(self, moduleType, offset=[0,0], angle=0):
		self.enabled = True 									# whether or not the module has enough available resources to function.
		self.active = False 									# if the module is specifically turned on, for example, an engine firing when the player pushes the up arrow key.
		self.offset = offset 									# x,y position on the ship
		self.moduleType = moduleType
		self.angle = angle										# the angle by which the module is slanted compared to the ships forward direction.
		self.effect = None										# a polygon drawn when the module is active, also can come with lights.
		self.mass = 1											# how much it weighs. The actor's mass is a sum of all its modules.

		if self.moduleType == 'generator':
			self.mass = 1
			self.active = True
			self.quiescent = {}
			self.resources = {
				'electricity': 1,								# positive valued resources are produced by the module.
				'fuel': -0.001									# negative valued resources are consumed by the module.
			}
			self.stores = {
				'fuel': 500,									# the stores field is the maximum amount the module can hold.
				'electricity': 500
			}
			self.initialStores = {
				'fuel': 500,									# initialStores is how much the module comes with when it is created.
				'electricity': 50
			}
			self.radius = 5										# radius is used to simplify physical calculations, for example with lighting, and as a baseline for the module's points.
			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			self.color = [75,10,10,255]							# the color fill of the module.
			self.outlineColor = [200,70,70,255]					# the outline color of the module, the outlines are also what receive illumination.

		elif self.moduleType == 'engine 10':
			self.mass = 1
			self.quiescent = {
				'electricity':0.001
			}
			self.resources = {
				'thrust': 100,
				'fuel': -1,
				'electricity':1
			}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			self.points = [[-self.radius, -self.radius*2], [-self.radius, self.radius*2], [self.radius,self.radius*2], [self.radius, -self.radius*2]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.effect = ModuleEffect('engine 10 flame', [0,self.radius*3])

		elif self.moduleType == 'engine 100':
			self.mass = 10
			self.quiescent = {
				'electricity':0.01
			}
			self.resources = {
				'thrust': 1000,
				'fuel': -10,
				'electricity':10
			}
			self.stores = {}
			self.initialStores = {}
			self.radius = 15
			self.points = [[-self.radius, -self.radius*2], [-self.radius, self.radius*2], [self.radius,self.radius*2], [self.radius, -self.radius*2]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.effect = ModuleEffect('engine 100 flame', [0,self.radius*3])

		elif self.moduleType == 'RCS':
			self.mass = 0.2
			self.quiescent = {
				'electricity':0.001
			}
			self.resources = {
				'torque': 5,
				'electricity': -0.2
			}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			self.color = [120,100,100,255]
			self.outlineColor = [170,150,150,255]
			self.momentArm = self.radius						# moment arm is the leverage by which this module can apply torque to the spacecraft. It is a distance.

		elif self.moduleType == 'spar 10':
			self.mass = 1
			self.quiescent = {}
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			self.points = [[-self.radius, -self.radius*10], [-self.radius, self.radius*10], [self.radius,self.radius*10], [self.radius, -self.radius*10]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'box 5':
			self.mass = 2
			self.quiescent = {}
			self.resources = {}
			self.stores = {'cargo':5}
			self.initialStores = {'cargo':0}
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size], [-size, size], [size,size], [size, -size]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'box 10':
			self.mass = 2
			self.quiescent = {}
			self.resources = {}
			self.stores = {'cargo':100}
			self.initialStores = {'cargo':0}
			self.radius = 5
			size = self.radius
			self.points = [[-size*5, -size*5], [-size*5, size*5], [size*5,size*5], [size*5, -size*5]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'box 100':
			self.mass = 100
			self.quiescent = {}
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*100, -size*100], [-size*100, size*100], [size*100,size*100], [size*100, -size*100]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'spar 100':
			self.mass = 1
			self.quiescent = {}
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*10, -size*100], [-size*10, size*100], [size*10,size*100], [size*10, -size*100]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'cannonshell 10':
			self.mass = 0.1
			self.active = True
			self.resources = {}
			self.quiescent = {}
			self.stores = {
				'high explosive': 10
			}
			self.initialStores = {
				'high explosive': 10
			}
			self.radius = 0.5
			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			self.color = [255,230,200,255]
			self.outlineColor = [255,235,225,255]
			self.lifetime = 100 										# the shell lasts for 1000 somethings and then explodes.

		elif self.moduleType == 'cannon 10':
			self.mass = 2
			self.active = False
			self.quiescent = {
				'electricity':0.001
			}
			self.resources = {
				'cannonshell 10': 1,
				'electricity':10
			}
			self.stores = {
				'cannonshell 10': 10
			}
			self.initialStores = {
				'cannonshell 10': 10
			}
			self.radius = 5
			self.points = [[-0.5*self.radius, -self.radius], [-0.5*self.radius, self.radius], [0.5*self.radius,self.radius], [0.5*self.radius, -self.radius]]
			self.color = [30,30,30,255]
			self.outlineColor = [100,100,30,255]

			self.barrelHole = [0,-self.radius + 1.5]					# the location that bullets are emitted from this module.
			self.muzzleVelocity = 250000								# the speed that bullets come out at.
			self.cooldownTime = 100										# how long it takes before the gun fires again.
			self.cooldownValue = 0
			self.effect = ModuleEffect('cannon 10 flash', [0,-self.radius])

		elif self.moduleType == 'hyperdrive 10':
			self.mass = 5
			self.quiescent = {
				'electricity':0.001
			}
			self.resources = {
				'electricity':-1000,
				'warp energy':1000
			}
			self.stores = {
				'warp energy':100
			}
			self.initialStores = {
				'warp energy':0
			}
			self.radius = 5
			size = self.radius
			self.points = [[-size*1.5, -size*2], [-size*1.5, size*2], [size*1.5,size*2], [size*1.5, -size*2]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]
			self.momentArm = self.radius

		elif self.moduleType == 'starbridge armor':
			self.mass = 2
			self.quiescent = {}
			self.resources = {}
			self.stores = {'armor':10}
			self.initialStores = {'armor':10}
			self.radius = 5
			size = self.radius
			self.points = [[-size*2, -size], [-size*2, size*2], [size*2,size*2], [size*2, -size*2]]
			self.color = [200,200,200,255]
			self.outlineColor = [240,240,240,255]
			self.momentArm = self.radius

		elif self.moduleType == 'turret 10':
			self.mass = 5
			self.active = False
			self.quiescent = {
				'electricity':0.01
			}
			self.resources = {
				'cannonshell 10': 1,
				'electricity':20
			}
			self.stores = {
				'cannonshell 10': 20
			}
			self.initialStores = {
				'cannonshell 10': 20
			}
			self.radius = 5
			self.points = [[-1*self.radius, -self.radius], [-2*self.radius, self.radius], [2*self.radius,self.radius], [1*self.radius, -self.radius]]
			self.color = [30,30,30,255]
			self.outlineColor = [100,100,30,255]

			self.barrelHole = [0,-self.radius + 1.5]					# the location that bullets are emitted from this module.
			self.muzzleVelocity = 250000								# the speed that bullets come out at.
			self.cooldownTime = 100										# how long it takes before the gun fires again.
			self.cooldownValue = 0
			self.effect = ModuleEffect('cannon 10 flash', [0,-self.radius])

			self.firingArc = 1/3* 2*math.pi

		elif self.moduleType == 'nano orbital':
			self.mass = 10000000
			self.quiescent = {}
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*1000, -size*5000], [-size*1000, size*5000], [size*1000,size*5000], [size*1000, -size*5000]]
			self.color = [60,150,255,255]
			self.outlineColor = [150,200,255,255]
			self.momentArm = self.radius