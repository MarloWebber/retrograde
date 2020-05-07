import math, sys, random
import pygame
from pygame.locals import *
from pygame.color import *
import pymunk
from pymunk import Vec2d
import pymunk.pygame_util
import numpy
import time

from AstroPynamics.tools.body import Body as APBody
from AstroPynamics.orbits.orbit import Orbit
from astropy.time import Time, TimeDelta
    
global contact
global shape_to_remove

def mag(x):
	return numpy.sqrt(x.dot(x))

def addRadians(a, b):
	return (a + b + math.pi) % (2*math.pi) - math.pi

def rotate_point(point, angle_rad, center_point=(0, 0)):
    """Rotates a point around center_point(origin by default)
    Angle is in degrees.
    Rotation is counter-clockwise
    """
    new_point = (point[0] - center_point[0], point[1] - center_point[1])  # Shift the point so that center_point becomes the origin
    new_point = (new_point[0] * math.cos(angle_rad) - new_point[1] * math.sin(angle_rad),
                 new_point[0] * math.sin(angle_rad) + new_point[1] * math.cos(angle_rad))
    new_point = [new_point[0] + center_point[0], new_point[1] + center_point[1]] # Reverse the shifting we have done
    return new_point

def rotate_polygon(polygon, angle, center_point=(0, 0)):
    """Rotates the given polygon which consists of corners represented as (x,y)
    around center_point (origin by default)
    Rotation is counter-clockwise
    Angle is in degrees """
    rotated_polygon = []
    for corner in polygon:
        rotated_corner = rotate_point(corner, angle, center_point)
        rotated_polygon.append(rotated_corner)
    return rotated_polygon

def area_of_annulus(inner, outer):
	return math.pi * ((inner**2) - (outer**2))

def make_circle(radius, points):
	# returns a list of points defining a circle.
	circle = []
	angleStep = 2*math.pi / points
	for i in range(0,points):
		circle.append([radius * math.cos(i * angleStep),radius * math.sin(i * angleStep) ])
	return circle

def mass_of_a_sphere(density, radius):
	return 4/3 * math.pi * radius**3

def semiMinorAxis(a, e):
	return a * math.sqrt(1 - (math.pow(e,2))) # https://math.stackexchange.com/questions/1259945/calculating-semi-minor-axis-of-an-ellipse

class Atmosphere():
	def __init__(self,radius, planetPosition):
		self.height = 2000
		self.density = 1
		self.color = (0,50,200)
		self.points = make_circle(radius+self.height, 314)
		self.mass = self.density * area_of_annulus(radius+self.height, radius)
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = planetPosition

class ModuleEffect(): # a ModuleEffect is just a polygon that is visually displayed when a module is active.
	def __init__(self, offset=[0,0]):
		self.offset = offset
		self.radius = 3
		self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]

class Module():
	def __init__(self, moduleType, offset=[0,0], angle=0):
		self.enabled = True # whether or not the module has enough available resources to function.
		self.active = False # if the module is specifically turned on, for example, an engine firing when the player pushes the up arrow key.
		self.offset = offset # x,y position on the ship
		self.moduleType = moduleType
		self.angle = angle

		if self.moduleType is 'generator':
			self.mass = 1
			self.active = True
			self.resources = {
				'electricity': 0.01,
				'fuel': -0.0001
			}
			self.stores = {
				'fuel': 5000,
				'electricity': 500
			}
			self.initialStores = {
				'fuel': 5000,
				'electricity': 50
			}
			self.mass = 1.5
			self.radius = 5
			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			self.color = [150,20,20]

		elif self.moduleType is 'engine':
			self.mass = 1
			self.resources = {
				'thrust': 100,
				'fuel': -0.1,
				'electricity':1
			}
			self.stores = {}
			self.initialStores = {}
			self.mass = 1.5
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size*2], [-size, size*2], [size,size*2], [size, -size*2]]
			self.color = [120,100,100]

			self.effect = ModuleEffect([0,0])

		elif self.moduleType is 'RCS':
			self.mass = 0.2
			self.resources = {
				'torque': 2,
				'electricity': -0.1
			}
			self.stores = {}
			self.initialStores = {}
			self.mass = 1.5
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size], [-size, size], [size,size], [size, -size]]
			self.color = [120,100,100]

			self.momentArm = self.radius

		for point in self.points:
				point[0] += self.offset[0]
				point[1] += self.offset[1]


dinghy = [Module('generator',[0,0]), Module('engine',[0,8]), Module('RCS',[0,-10]) ]
lothar = [Module('generator',[0,0]), Module('engine',[-13,8]), Module('engine',[13,8]), Module('RCS',[-13,-10]), Module('RCS',[13,-10]) ]

class Actor():
	def __init__(self, name, modulesList, position, velocity, isPlayer=False):
		self.type = 'actor'
		self.name = name #str(random.randint(0,1000)) # the individual name of the craft. Set to shipType for now.
		self.modules = []

		for module in modulesList:
			self.modules.append(module)

		self.mass = 0
		self.points = []
		self.availableResources = {}
		self.storagePool = {}
		self.maximumStores = {}

		self.isPlayer = isPlayer

		for module in self.modules:
			self.mass += module.mass
			self.points += module.points
			for resource, quantity in list(module.initialStores.items()):
				if resource not in self.storagePool: self.storagePool[resource] = quantity
				else: self.storagePool[resource] += quantity
			for resource, quantity in list(module.initialStores.items()):
				if resource not in self.availableResources: self.availableResources[resource] = quantity
				else: self.availableResources[resource] += quantity
			for resource, quantity in list(module.stores.items()):
				if resource not in self.maximumStores: self.maximumStores[resource] = quantity
				else: self.maximumStores[resource] += quantity
		
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.body.apply_impulse_at_local_point(velocity, [0,0])
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = 0.9
		self.orbit = None #initpos_to_orbit(self.)
		self.freefalling = True
		self.color = (200,50,50)
		self.keyStates = {
			'up': False,
			'down': False,
			'left': False,
			'right': False
		}
		self.orbiting = None
		self.stepsToFreefall = 5
		self.decompVelocity = 500000
		# self.prevVelocity = [velocity[0], velocity[1]]

	def leaveFreefall(self):
		self.stepsToFreefall = 10
		self.freefalling = False
		self.orbit = None
		
	def doResources(self):
		# - tally the amount of stored resources. first, zero everything out
		for module in self.modules:
			for resource, quantity in list(module.resources.items()):
					self.availableResources[resource] = 0

		# add all the unstored resources being produced by active modules
		for module in self.modules:
			for resource, quantity in list(module.resources.items()):
				if quantity > 0 and module.enabled and module.active:
					if resource not in self.availableResources:
						self.availableResources[resource] = quantity
					else:
						self.availableResources[resource] += quantity

		# - turn modules on and off
		for module in self.modules:
			module.enabled = True
			for resource, quantity in list(module.resources.items()):
				if quantity < 0:
					if resource in self.storagePool:
						availableAmount = self.storagePool[resource] + self.availableResources[resource]
					else:
						availableAmount = self.availableResources[resource]
					if availableAmount < abs(quantity):
						module.enabled = False
			if module.moduleType == 'RCS':
				module.active = self.keyStates['left'] or self.keyStates['right']
			elif module.moduleType == 'engine':
				module.active = self.keyStates['up']

		# - consume and produce resources
		for module in self.modules:
			if module.enabled and module.active:
				for resource, quantity in list(module.resources.items()):
					if quantity > 0: # producing resource
						if resource in self.storagePool:
							remainingCapacity = self.maximumStores[resource] - self.storagePool[resource]
							if quantity > remainingCapacity:
								self.storagePool[resource] += remainingCapacity
								self.availableResources[resource] += quantity - remainingCapacity
							else:
								self.storagePool[resource] += quantity

					else: # consuming resource
						self.availableResources[resource] += quantity # adding a negative number is a subtraction
						if self.availableResources[resource] < 0:
							self.storagePool[resource] += self.availableResources[resource]
							self.availableResources[resource] = 0

		for resource, quantity in list(self.storagePool.items()):
			if quantity > self.maximumStores[resource]: quantity = self.maximumStores[resource]			

	def doModuleEffects(self, keyStates, timestepSize):
		ifThrustHasBeenApplied = False
		for module in self.modules:
			if module.enabled and module.active:
				for giveResource, giveQuantity in list(module.resources.items()): #module.produces.items():
					if giveResource == 'thrust':
						if keyStates['up']:
							force = [(giveQuantity * timestepSize * 500 * math.cos(addRadians(module.angle, math.pi * 0.5))), -giveQuantity * timestepSize * 500 * math.sin(addRadians(module.angle, math.pi * 0.5) )]
							self.body.apply_impulse_at_local_point(force, (0,0))
							ifThrustHasBeenApplied = True

					elif giveResource == 'torque':
						if keyStates['left']:
							# apply two impulses, pushing in opposite directions, an equal distance from the center to create torque
							self.body.apply_impulse_at_local_point([-giveQuantity,0], [0,-module.momentArm])
							self.body.apply_impulse_at_local_point([giveQuantity,0], [0,module.momentArm])
						elif keyStates['right']:
							self.body.apply_impulse_at_local_point([giveQuantity,0], [0,-module.momentArm])
							self.body.apply_impulse_at_local_point([-giveQuantity,0], [0,module.momentArm])

		return ifThrustHasBeenApplied

class Attractor():
	def __init__(self, planetType, position, gravitationalConstant):
		self.type = 'attractor' # actor or attractor type. this is very stupid and should be designed out
		self.planetName = planetType # just set the planetName to something easy for now. # the name of the individual instance of this planet type.
		
		if planetType == 'earth':
			self.radius = 320000
			self.density = 0.5
			self.friction = 0.9
			self.color = (190,165,145)
			self.atmosphere = Atmosphere(self.radius, position)

		elif planetType == 'moon':
			self.radius = 80000
			self.density = 0.75
			self.friction = 0.9
			self.color = (145,145,145)
			self.atmosphere = None #Atmosphere(self)

		# create pymunk physical body and shape
		self.mass = mass_of_a_sphere(self.density, self.radius)
		size = self.radius
		self.points = make_circle(self.radius, 120)
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = self.friction

		# create astropynamics orbit-able body
		self.APBody = APBody(self.planetName, self.mass * gravitationalConstant * 0.1655, self.radius)

class World():
	def __init__(self):
		pygame.init()
		self.clock = pygame.time.Clock() # the pygame clock is NOT the same as the simulation clock.
		self.time = 0 # the number of timesteps that have passed in-universe. used for physics and orbital calculations.
		self.space = pymunk.Space()
		self.space.gravity = (0.0, 0.0)
		self.gravitationalConstant = 0.03
		self.actors = []
		self.attractors = []
		self.resolution = (1280,780)
		self.screen = pygame.display.set_mode(self.resolution)
		self.draw_options = pymunk.pygame_util.DrawOptions(self.screen)
		self.ch = self.space.add_collision_handler(0, 0)
		self.ch.data["surface"] = self.screen
		self.ch.post_solve = self.handle_collision
		self.viewpointObject = None
		self.viewpointObjectIndex = 0
		self.player = None
		self.zoom = 1 # the actual applied zoom number.
		self.pan = [0,0]
		self.rotate = 0
		self.timestepSize = 0.2/60.0 #1.0/60.0
		pygame.key.set_repeat(50,50) # holding a key down repeats the instruction. https://www.pygame.org/docs/ref/key.html
		self.font = pygame.font.SysFont('dejavusans', 15)
		self.showHUD = False

	def gravityForce(self, actorPosition, attractorPosition, attractorMass):
		distance = attractorPosition - actorPosition # scalar distance between two bodies
		magnitude = mag(distance)
		gravity = self.gravitationalConstant * attractorMass
		appliedGravity = gravity/(magnitude**2)
		components = numpy.divide(distance, magnitude)
		force = components * appliedGravity * self.timestepSize
		return force

	def gravitate(self, actor, force):
		rotatedForce = Vec2d(force[0], force[1])
		rotatedForce = rotatedForce.rotated(-actor.body.angle)
		actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])

	def inputs(self):
		for event in pygame.event.get():
			if event.type == KEYDOWN and event.key == K_ESCAPE:
				self.running = False
			elif event.type == KEYDOWN and event.key == K_RIGHTBRACKET:
				self.viewpointObjectIndex += 1
				if self.viewpointObjectIndex >= len(self.actors):
					self.viewpointObjectIndex = 0
				self.viewpointObject = self.actors[self.viewpointObjectIndex]
			elif event.type == KEYDOWN and event.key == K_LEFTBRACKET:
				self.viewpointObjectIndex -= 1
				if self.viewpointObjectIndex < 0:
					self.viewpointObjectIndex = len(self.actors) - 1
				self.viewpointObject = self.actors[self.viewpointObjectIndex]
			elif event.type == KEYDOWN and event.key == K_EQUALS:
				self.zoom += self.zoom * 0.5
			elif event.type == KEYDOWN and event.key == K_MINUS:
				self.zoom -= self.zoom * 0.5
			elif event.type == KEYDOWN and event.key == K_COMMA:
				self.timestepSize += self.timestepSize * 0.5
			elif event.type == KEYDOWN and event.key == K_PERIOD:
				self.timestepSize -= self.timestepSize * 0.5
			elif event.type == KEYDOWN and event.key == K_LEFT:
				self.player.keyStates['left'] = True
			elif event.type == KEYUP and event.key == K_LEFT:
				self.player.keyStates['left'] = False
			elif event.type == KEYDOWN and event.key == K_RIGHT:
				self.player.keyStates['right'] = True
			elif event.type == KEYUP and event.key == K_RIGHT:
				self.player.keyStates['right'] = False
			elif event.type == KEYDOWN and event.key == K_UP:
				self.player.keyStates['up'] = True
			elif event.type == KEYUP and event.key == K_UP:
				self.player.keyStates['up'] = False
			elif event.type == KEYDOWN and event.key == K_DOWN:
				self.player.keyStates['down'] = True
			elif event.type == KEYUP and event.key == K_DOWN:
				self.player.keyStates['down'] = False
			elif event.type == KEYDOWN and event.key == K_h:
				self.showHUD = not self.showHUD
        
	def add(self, thing):  
		self.space.add(thing.body, thing.shape)
		if thing.type == 'actor':
			self.actors.append(thing)
		elif thing.type == 'attractor':
			self.attractors.append(thing)

		self.player = self._getPlayer()
		self.viewpointObject = self.player

	def decomposeActor(self, actor, modules):

		if len(actor.modules) == 1:
			return # the actor is already fully decomposed, destroy it if you want
		else:
			# remove the actor from the space.
			self.space.remove(actor.body)
			self.space.remove(actor.shape)

			if actor in self.actors:
				self.actors.remove(actor)

			# create a new actor, minus the module
			for module in modules:
				actor.modules.remove(module)

				# create the module on it's own as a new actor
				fragmentPosition = [actor.body.position[0] + (module.offset[0] * math.cos(actor.body.angle)), actor.body.position[1] +  (module.offset[1] * math.sin(actor.body.angle))]
				self.add(Actor(actor.name + ' fragment', [module], fragmentPosition, actor.body.velocity, False))

			if not actor.modules:
				pass
			else:
				self.add( Actor(actor.name, actor.modules, actor.body.position, actor.body.velocity, actor.isPlayer ) ) # add the remaining parts of the original actor back into the space

	def physics(self):
		for actor in self.actors:

			# actor.prevVelocity = actor.body.velocity

			actor.doResources()

			# figure out which attractor you are orbiting
			strongestForce = None
			strongestAttractor = None
			for attractor in self.attractors:
				force = self.gravityForce(actor.body.position, attractor.body.position, attractor.body.mass)
				if strongestAttractor is None or mag(force) > mag(strongestForce):
					strongestForce = force
					strongestAttractor = attractor

			if strongestAttractor is not actor.orbiting or actor.orbiting is None:
				actor.orbiting = strongestAttractor
				
			# figure out if the actor is freefalling by seeing if any engines or collisions have moved it.
			if actor.doModuleEffects(actor.keyStates, self.timestepSize):
				actor.leaveFreefall()

			# if it is freefalling, move it along the orbital track.
			if actor.freefalling and actor.orbit is not None:
				actor.orbit.updTime(self.timestepSize)
				cartesian = actor.orbit.cartesianCoordinates(actor.orbit.tAn)
				actor.body.position =  [cartesian[0] ,cartesian[1] ]

				# you also must update the actor's velocity, or else when it leaves the track it will have the same velocity it entered with, leading to weird jumps.
				trackSpeed = actor.orbit.getSpeed(actor.orbit.tAn)

				# there is almost definitely a way to figure this out from the ellipse's properties. You would need to find tangent to the ellipse. But I figured it out by computing the position one step into the future, and then finding the angle between positions.
				futureSteptAn = actor.orbit.tAnAtTime(self.timestepSize)
				futureStepCoordinates = actor.orbit.cartesianCoordinates(futureSteptAn)
				angle = math.atan2( futureStepCoordinates[1] - actor.body.position[1], futureStepCoordinates[0] - actor.body.position[0] )

				actor.body.velocity = [trackSpeed * math.cos(angle), trackSpeed * math.sin(angle)]

			else:
				# if it is not freefalling, add some gravity to it so that it moves naturally, and try to recalculate the orbit.
				self.gravitate(actor, strongestForce)
				if actor.stepsToFreefall > 0:
					actor.stepsToFreefall -= 1
				else:
					actor.freefalling = True
					try:
						actor.orbit = Orbit.fromStateVector(numpy.array([actor.body.position[0],actor.body.position[1],1]), numpy.array([actor.body.velocity[0],actor.body.velocity[1],1]), actor.orbiting.APBody, Time('2000-01-01 00:00:00'), actor.name + " orbit around " + actor.orbiting.planetName)
					except:
						pass
				
		self.space.step(self.timestepSize)
		self.time += self.timestepSize

	def rotatePolygon(self, points, angle):
		return Rotate2D(points,(0,0),angle)

	def transformForView(self, position):
		if self.viewpointObject == None:
			return position
		else:
			transformedPosition = position - self.viewpointObject.body.position # offset everything by the position of the viewpointObject, so the viewpoint is at 0,0
			transformedPosition = transformedPosition * self.zoom  # shrink or expand everything around the 0,0 point
			transformedPosition[0] += 0.5 * self.resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
			transformedPosition[1] += 0.5 * self.resolution[1] # add half the height.
			return transformedPosition

	def drawCircle(self,color, position, radius):
		pygame.draw.circle(self.screen, color, [int(position[0]), int(position[1])], int((radius * self.zoom)))

	def drawActor(self, actor):
		rotatedPoints = rotate_polygon(actor.points,actor.body.angle)  # orient the polygon according to the body's current direction in space.
		transformedPoints = []
		for rotatedPoint in rotatedPoints:
			transformedPoints.append(self.transformForView(rotatedPoint + actor.body.position)) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
		pygame.draw.lines(self.screen, actor.color, True, transformedPoints)

	def drawModule(self, actor, module):
		# draw the outline of the module.
		transformedPoints = []
		for rotatedPoint in module.points: 
			transformedPoints.append(self.transformForView(rotatedPoint + actor.body.position )) # zoom and pan to fit the screen.
		rotatedPoints = rotate_polygon(transformedPoints,actor.body.angle, self.transformForView(actor.body.position))  # orient the module according to the actor's current direction in space.
		pygame.draw.lines(self.screen, module.color, True, rotatedPoints)

		# put a circle in the middle if it is enabled, and a smaller red circle in the middle of that, if it is activated.
		if module.enabled:
			activeCircle = self.transformForView(module.offset + actor.body.position)
			activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
			self.drawCircle(module.color, activeCircle, 2)
			if module.active:
				self.drawCircle((255,0,0), activeCircle, 1)

		# rocket engines have a line coming out of them.
		if module.enabled and module.active:
			for giveResource, giveQuantity in list(module.resources.items()): 
				if giveResource == 'thrust':
					if actor.keyStates['up']:
						forceAngle = actor.body.angle
						force = [(giveQuantity * math.cos(addRadians(forceAngle, math.pi * 0.5))), giveQuantity * math.sin(addRadians(forceAngle, math.pi * 0.5) )]
						activeCircle = self.transformForView(module.offset + actor.body.position)
						activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
						ananas = (int(activeCircle[0] + force[0] * self.zoom), int(activeCircle[1]+force[1] * self.zoom ) )
						# print ananas
						pygame.draw.lines(self.screen, (255,255,200), True, [activeCircle,ananas])


	def getActorFromBody(self, body):
		for actor in self.actors:
			if body is actor.body:
				return actor
		return None

	def handle_collision(self, arbiter, space, data):
		# for c in arbiter.contact_point_set.points:
		# 	r = max( 3, abs(c.distance*5) )
		# 	r = int(r)
		# 	p = tuple(map(int, c.point_a))
			# pygame.draw.circle(data["surface"], THECOLORS["red"], p, r, 0)
			# self.drawCircle((255,255,255), self.transformForView(p), 100) # this doesn't even work, i don't know why
			# print( self.transformForView(p))

		shapes = arbiter._get_shapes()
		# for shape in shapes:
		bodyA = shapes[0]._get_body()
		actorA = self.getActorFromBody(bodyA)

		bodyB = shapes[0]._get_body()
		actorB = self.getActorFromBody(bodyB)

		# get difference in velocity
		# deltaV = actorA.prevVelocity - actorB.prevVelocity
		# print deltaV
		ke = arbiter._get_total_ke()
		ke_A = ke / actorA.mass
		ke_B = ke / actorB.mass

		# print(ke_A)

		actorA.leaveFreefall()
		actorB.leaveFreefall()

		if ke_A > actorA.decompVelocity:
			self.decomposeActor(actorA, actorA.modules)
		if ke_B > actorB.decompVelocity:
			self.decomposeActor(actorB, actorB.modules)


			# print(deltaV)
			# if deltaV > actor.decompVelocity:
			# 	self.decomposeActor(actor, actor.modules)

	def _getPlayer(self):
		for actor in self.actors:
			if actor.isPlayer:
				return actor


	def drawAPOrbit(self, orbit, color):

		points = []


		# print("period " + str(orbit.getPeriod()))
		# print("apoapsis " + str(orbit.getApoapsis()))
		# print("periapsis " + str(orbit.getPeriapsis()))
		# print("crashing " + str(orbit.isCrashing()))

		for i in range(0,100):
			temp_vec3d = orbit.cartesianCoordinates(i * (2 * math.pi / 100))
			# print(temp_vec3d)
			point = (temp_vec3d[0], temp_vec3d[1])
			point = self.transformForView(point)
			points.append(point)

		for i in range(0,100):
			if i > 0:
				pygame.draw.lines(self.screen, color, True, (points[i-1], points[i]))


	


		# self.drawCircle((255,0,0),self.transformForView(chickybabe), 5000)


	def drawHUD(self):

		# self.drawEllipse(self.player.orbit, self.attrplayer)

		# show the player what resources are available
		i = 0
		hudList = {}
		for resource, quantity in list(self.player.availableResources.items()):
			hudList[resource] = quantity
		for resource, quantity in list(self.player.storagePool.items()):
			if resource in hudList:
				hudList[resource] += self.player.storagePool[resource]
			else:
				hudList[resource] = self.player.storagePool[resource]

		for availableResource, availableQuantity in list(hudList.items()):
			textsurface = self.font.render(str(availableResource) + ': ' + str(availableQuantity), False, (255, 255, 255))
			self.screen.blit(textsurface,(30,i * 20))
			i += 1

		textsurface = self.font.render('warp: ' + str(self.timestepSize * 3 * 100), False, (255, 255, 255))
		self.screen.blit(textsurface,(30,i * 20))
		i += 1
		textsurface = self.font.render('zoom: ' + str(self.zoom), False, (255, 255, 255))
		self.screen.blit(textsurface,(30,i * 20))
		i += 1
		textsurface = self.font.render('time: ' + str(self.time), False, (255, 255, 255))
		self.screen.blit(textsurface,(30,i * 20))
		i += 1
		textsurface = self.font.render('freefalling: ' + str(self.player.freefalling), False, (255, 255, 255))
		self.screen.blit(textsurface,(30,i * 20))




		# print the navcircle
		n_navcircle_lines = 32
		navcircleLinesLength = 10
		navcircleInnerRadius = 250

		for n in range(0,n_navcircle_lines):
			angle = n * (2 * math.pi / n_navcircle_lines)
			start = ((navcircleInnerRadius * math.cos(angle)) + (self.resolution[0]*0.5) , (navcircleInnerRadius* math.sin(angle)) +( self.resolution[1] * 0.5) )
			end = ((navcircleInnerRadius + navcircleLinesLength) * math.cos(angle)+ (self.resolution[0]*0.5), (navcircleInnerRadius + navcircleLinesLength) * math.sin(angle)+ (self.resolution[1]*0.5))
			# navcircleLines.append([start, end])
			pygame.draw.lines(self.screen, (100,100,100), True, (start,end))

		blipLength = (navcircleInnerRadius-navcircleLinesLength)
		angle = self.player.body.angle - 0.5 * math.pi
		start = ((blipLength * math.cos(angle)) + (self.resolution[0]*0.5) , (blipLength* math.sin(angle)) +( self.resolution[1] * 0.5) )
		end = ((navcircleInnerRadius) * math.cos(angle)+ (self.resolution[0]*0.5), (navcircleInnerRadius) * math.sin(angle)+ (self.resolution[1]*0.5))
		pygame.draw.lines(self.screen, (200,0,10), True, (start,end))

		# draw the trajectory
		# n_trajectory_points = 100
		# trajectory_point_position = self.player.body.position
		# trajectory_point_velocity = self.player.body.velocity
		# prev_trajectory_point = trajectory_point_position
		# for n in xrange(0,n_trajectory_points):
		# 	sumGravity = [0,0]
		# 	for attractor in self.attractors:
		# 		sumGravity += self.gravityForce(self.player.body.position, attractor.body.position, attractor.body.mass) 
		# 	prev_trajectory_point = trajectory_point_position
		# 	trajectory_point_position += trajectory_point_velocity * self.timestepSize
		# 	sumGravity[0] = sumGravity[0] * self.timestepSize
		# 	sumGravity[1] = sumGravity[1] * self.timestepSize
		# 	trajectory_point_velocity += sumGravity 
		# 	pygame.draw.lines(self.screen, (200,0,10), True, (self.transformForView(prev_trajectory_point),self.transformForView(trajectory_point_position)))


	def graphics(self):
		### Clear screen
		self.screen.fill(THECOLORS["black"])

		### Draw stuff
		for attractor in self.attractors:
			if attractor.atmosphere != None:
				self.drawActor(attractor.atmosphere)

				# atmosphere_pos = self.transformForView([-attractor.radius + attractor.body.position[0], -attractor.radius + attractor.body.position[1]])
				# self.screen.blit(gradients.radial( (attractor.radius + attractor.atmosphere.height) * self.zoom  , (0,255,0,100), (100,0,50,255)), atmosphere_pos )

			self.drawActor(attractor)
		for actor in self.actors:
			if actor.orbit is not None:
				self.drawAPOrbit(actor.orbit, (100,100,100))
			for module in actor.modules:
				self.drawModule(actor, module)

		# self.drawAPOrbit(self.calibrationOrbit, (100,100,100))

		self.drawHUD()


		pygame.display.flip()
		self.clock.tick(150)
		pygame.display.set_caption("fps: " + str(self.clock.get_fps()))

	def step(self):
		self.player = self._getPlayer()
		self.inputs()
		self.physics()
		self.graphics()

	def start(self):

		planet_erf = Attractor('earth', [1,1], self.gravitationalConstant)
		planet_moon = Attractor('moon',[-1000000,-1000000], self.gravitationalConstant)
		self.add(planet_erf)
		self.add(planet_moon)

		dinghy_instance = Actor('player Lothar', lothar,(10, -322100), [300000,0], True)
		lothar_instance = Actor('NPC lothar', lothar,(1, -322160), [300000,0])
		self.add(dinghy_instance)
		self.add(lothar_instance)
		
		# self.calibrationOrbit = Orbit.fromStateVector(numpy.array([self.player.body.position[0],self.player.body.position[1],1]), numpy.array([self.player.body.velocity[0],self.player.body.velocity[1],1]), self.attrplayer.APBody, Time('2000-01-01 00:00:00'), "calibration orbit")
				

		self.running = True
		while self.running:
			self.step()

mundus = World()
mundus.start()
