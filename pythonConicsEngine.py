import math, sys, random

import pygame
from pygame.locals import *
from pygame.color import *

import pymunk
from pymunk import Vec2d
import pymunk.pygame_util

import numpy

import initpos_to_orbit as orbmech

# def draw_collision(arbiter, space, data):
#     for c in arbiter.contact_point_set.points:
#         r = max( 3, abs(c.distance*5) )
#         r = int(r)

#         p = pymunk.pygame_util.to_pygame(c.point_a, data["surface"])
#         pygame.draw.circle(data["surface"], THECOLORS["black"], p, r, 1)
    
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
    # angle_rad = math.radians(angle % 360)
    # Shift the point so that center_point becomes the origin
    new_point = (point[0] - center_point[0], point[1] - center_point[1])
    new_point = (new_point[0] * math.cos(angle_rad) - new_point[1] * math.sin(angle_rad),
                 new_point[0] * math.sin(angle_rad) + new_point[1] * math.cos(angle_rad))
    # Reverse the shifting we have done
    new_point = [new_point[0] + center_point[0], new_point[1] + center_point[1]]
    return new_point

def rotate_polygon(polygon, angle, center_point=(0, 0)):
    """Rotates the given polygon which consists of corners represented as (x,y)
    around center_point (origin by default)
    Rotation is counter-clockwise
    Angle is in degrees
    """
    rotated_polygon = []
    for corner in polygon:
        rotated_corner = rotate_point(corner, angle, center_point)
        rotated_polygon.append(rotated_corner)
    return rotated_polygon

def area_of_annulus(inner, outer):
	return math.pi * ((inner*inner) - (outer*outer))

def make_circle(radius, points):
	# returns a list of points defining a circle.
	circle = []
	angleStep = 2*math.pi / points
	for i in xrange(0,points):
		circle.append([radius * math.cos(i * angleStep),radius * math.sin(i * angleStep) ])
	return circle


class Orbit():
	# a 		Semi Major Axis
	# e 		Eccentricity
	# i 		Inclination
	# Omega 	Longitude of Ascending Node
	# omega 	Argument of Periapsis
	# nu 		True Anomaly
	def __init__(a,e,i,Omega,omega,nu):
		self.a = a
		self.e = e
		self.i = i
		self.Omega = Omega
		self.omega = omega
		self.nu = nu

# 	def semiMinorAxis(self):
# 		self.b = self.a * math.sqrt(1 - (math.pow(self.e,2))) # https://math.stackexchange.com/questions/1259945/calculating-semi-minor-axis-of-an-ellipse

# 	def semiLatusRectum(self):
# 		self.l = (math.pow(self.b,2) / self.a) # https://en.wikipedia.org/wiki/Ellipse

class Atmosphere():
	def __init__(self,planet):
		self.height = 2000
		self.density = 1

		self.color = (0,50,200)

		# self.seaLevelDrag = 0.1

		self.points = make_circle(planet.radius+self.height, 314)

		self.mass = self.density * area_of_annulus(planet.radius+self.height, planet.radius)
		# print self.mass
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = planet.body.position

class ModuleEffect(): # a ModuleEffect is just a polygon that is visually displayed when a module is active.
	def __init__(self, offset=[0,0]):
		self.offset = offset
		self.radius = 3
		self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]

class Module():
	def __init__(self, moduleType, offset=[0,0], angle=0):

		self.enabled = False # whether or not the module has enough available resources to function.
		self.active = False # if the module is specifically turned on, for example, an engine firing when the player pushes the up arrow key.

		self.offset = offset # x,y position on the ship
		self.moduleType = moduleType

		self.angle = angle

		if self.moduleType is 'generator':
			self.mass = 1
			self.active = True

			self.produces = {
				'electricity': 0.01
			}
			self.consumes = {
				'fuel': 0.00001
			}
			self.stores = {
				'fuel': 100,
				'electricity': 100
			}
			self.currentStores = {
				'fuel': 100,
				'electricity': 50
			}

			self.mass = 1.5
			self.radius = 5

			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			for point in self.points:
				point[0] += self.offset[0]
				point[1] += self.offset[1]

			self.color = [200,20,20]


		elif self.moduleType is 'engine':
			self.mass = 1

			self.produces = {
				'thrust': 50
			}
			self.consumes = {
				'fuel': 0.1
			}
			self.stores = {
				'thrust': 10
			}
			self.currentStores = {
				'thrust': 1,
			}

			self.mass = 1.5
			self.radius = 5

			size = self.radius
			self.points = [[-size, -size*2], [-size, size*2], [size,size*2], [size, -size*2]]
			for point in self.points:
				point[0] += self.offset[0]
				point[1] += self.offset[1]

			self.color = [120,100,100]

			self.effect = ModuleEffect([0,0])

		elif self.moduleType is 'RCS':
			self.mass = 0.2

			self.produces = {
				'torque': 1
			}
			self.consumes = {
				'electricity': 0.1
			}
			self.stores = {
				'torque': 1
			}
			self.currentStores = {
				'torque': 0,
			}

			self.mass = 1.5
			self.radius = 5

			size = self.radius
			self.points = [[-size, -size], [-size, size], [size,size], [size, -size]]
			for point in self.points:
				point[0] += self.offset[0]
				point[1] += self.offset[1]

			self.color = [120,100,100]

			self.momentArm = self.radius


class Actor():
	def __init__(self, position, velocity):
		self.type = 'actor'

		self.mass = 0
		self.points = []

		self.modules = []
		self.modules.append(Module('generator',[0,0]))
		self.modules.append(Module('engine',[0,8]))
		self.modules.append(Module('RCS',[0,-5.5]))

		for module in self.modules:
			self.mass += module.mass
			self.points += module.points

		# print self.mass
		# print self.points

		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.body.apply_impulse_at_local_point(velocity, [0,0])
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = 0.5

		self.orbit = None #initpos_to_orbit(self.)

		self.freefalling = True
		self.interacting = False

		self.color = (200,50,50)

		# self.image = None #preprocessed image used to 'blit' onto the screen.

		self.availableResources = {}

		# self.isPlayer = False
		self.keyStates = {
			'up': False,
			'down': False,
			'left': False,
			'right': False
		}


	# def enterFreefall(self, attractor):
	# 	self.orbit = func_initpos_to_orbit()
	# 	self.freefalling = True

	# def leaveFreefall(self):
	# 	self.freefalling = False

	def getAvailableResources(self):
		for module in self.modules:
			for availableResource in module.currentStores.keys():
				self.availableResources[availableResource] = 0

		for module in self.modules:
			for availableResource, availableQuantity in module.currentStores.items():
				# print availableResource
				# print availableQuantity
				# print self.availableResources
				if self.availableResources[availableResource] == None:
					self.availableResources[availableResource] = availableQuantity
				else:
					self.availableResources[availableResource] += availableQuantity

		for module in self.modules:		# first, allow all the modules to produce what they need, so you know what is available.
			if module.enabled and module.active:
				for resource, quantity in module.produces.items():
					self.availableResources[resource] += quantity


	def doResourceConsumption(self):
		for module in self.modules:
			if module.active:
				module.enabled = False
				for needResource, needQuantity in module.consumes.items():
					for providerModule in self.modules:
						if needResource in providerModule.currentStores and providerModule.currentStores[needResource] > 0:
							if providerModule.currentStores[needResource] < needQuantity:
								needQuantity -= providerModule.currentStores[needResource]
								providerModule.currentStores[needResource] = 0
							else:
								providerModule.currentStores[needResource] -= needQuantity
								module.enabled = True
								break

	def doResourceStorage(self):
		for module in self.modules:
			if module.active and module.enabled:
				for giveResource, giveQuantity in module.produces.items():
					for accepterModule in self.modules:
						if giveResource in accepterModule.stores:
							remainingCapacity = accepterModule.stores[giveResource] - accepterModule.currentStores[giveResource]
							if remainingCapacity > 0:
								if remainingCapacity > giveQuantity:
									accepterModule.currentStores[giveResource] += giveQuantity
									break
								else:
									accepterModule.currentStores[giveResource] = accepterModule.stores[giveResource]
									giveQuantity -= remainingCapacity

	def doModuleActivation(self):
		for module in self.modules:
			if module.moduleType == 'RCS':
				module.active = self.keyStates['left'] or self.keyStates['right']
			if module.moduleType == 'engine':
				module.active = self.keyStates['up']

	def doResources(self):
		self.doModuleActivation()
		self.getAvailableResources()
		self.doResourceConsumption()
		self.doResourceStorage()

		

		

		
									

	def doModuleEffects(self, keyStates):
		for module in self.modules:
			if module.enabled and module.active:
				for giveResource, giveQuantity in module.stores.items(): #module.produces.items():
					if giveResource == 'thrust':
						if keyStates['up']:
							# force = [(giveQuantity * math.sin(module.angle + self.body.angle)), giveQuantity * math.cos(module.angle + self.body.angle) ]
							# self.body.apply_impulse_at_local_point(force, module.offset)

							# forceAngle = addRadians(actor.body.angle, 0)
							forceAngle = self.body.angle#addRadians(module.angle, actor.body.angle)
							# print forceAngle

							force = [(giveQuantity * math.cos(addRadians(forceAngle, math.pi * 0.5))), -giveQuantity * math.sin(addRadians(forceAngle, math.pi * 0.5) )]
							self.body.apply_impulse_at_local_point(force, (0,0))


					elif giveResource == 'torque':
						if keyStates['left']:
							# apply two impulses, pushing in opposite directions, an equal distance from the center to create torque
							self.body.apply_impulse_at_local_point([-giveQuantity,0], [0,-module.momentArm])
							self.body.apply_impulse_at_local_point([giveQuantity,0], [0,module.momentArm])
						elif keyStates['right']:
							# apply two impulses, pushing in opposite directions, an equal distance from the center to create torque
							self.body.apply_impulse_at_local_point([giveQuantity,0], [0,-module.momentArm])
							self.body.apply_impulse_at_local_point([-giveQuantity,0], [0,module.momentArm])





class Attractor():
	def __init__(self):
		self.type = 'attractor'
		
		self.radius = 320000
		self.density = 0.5

		self.mass = self.density * (math.pi * (self.radius * self.radius))

		# size = self.radius
		self.points = make_circle(self.radius, 314)
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		# self.body.position = position
		self.body.position = 1, 1
		# self.body.apply_impulse_at_local_point(velocity, [0,0])
		self.shape = pymunk.Poly(self.body, self.points)


		#inertia = pymunk.moment_for_circle(self.mass, 0, self.radius, (0,0))
		#self.body = pymunk.Body(self.mass, inertia)
		
		#self.shape = pymunk.Circle(self.body, self.radius, (0,0))
		self.shape.friction = 0.5

		self.color = (190,165,145)

		self.atmosphere = Atmosphere(self)

		

class World():
	def __init__(self):
		pygame.init()

		self.clock = pygame.time.Clock()
		self.space = pymunk.Space()
		self.space.gravity = (0.0, 0.0)
		self.gravitationalConstant = 300

		self.actors = []
		self.attractors = []

		self.resolution = (840,680)
		self.screen = pygame.display.set_mode(self.resolution)
		self.draw_options = pymunk.pygame_util.DrawOptions(self.screen)
		self.draw_options.flags = self.draw_options.flags ^ pymunk.pygame_util.DrawOptions.DRAW_COLLISION_POINTS 

		self.ch = self.space.add_collision_handler(0, 0)
		self.ch.data["surface"] = self.screen
		# self.ch.post_solve = draw_collision

		self.viewpointObject = None
		self.player = None
		self.zoom = 1
		self.pan = [0,0]
		self.rotate = 0

		self.timestepSize = 0.2/60.0 #1.0/60.0

		pygame.key.set_repeat(50,50) # holding a key down repeats the instruction. https://www.pygame.org/docs/ref/key.html

		self.font = pygame.font.SysFont('dejavusans', 15)

		# playerKey states are used to say whether the key is being held down or not. This is because the recommended pygame key event code only provides keyup and keydown events.
		# self.playerKeyStates = {
		# 	'up': False,
		# 	'down': False,
		# 	'left': False,
		# 	'right': False
		# }


	def gravitate(self, actor, attractor):
		distance = attractor.body.position - actor.body.position # scalar distance between two bodies
		magnitude = mag(distance)
		gravity = self.gravitationalConstant * attractor.body.mass
		appliedGravity = gravity/(magnitude * magnitude)

		components = numpy.divide(distance, magnitude)

		force = components * appliedGravity * self.timestepSize

		# print force

		rotatedForce = Vec2d(force[0], force[1])
		rotatedForce = rotatedForce.rotated(-actor.body.angle)

		actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])

	
	def inputs(self):
		for event in pygame.event.get():
			if event.type == KEYDOWN and event.key == K_ESCAPE:
				self.running = False
			elif event.type == KEYDOWN and event.key == K_RIGHTBRACKET:
				if self.viewpointObject == self.actors[0]:
					self.viewpointObject = self.attractors[0]
				else:
					self.viewpointObject = self.actors[0]

			elif event.type == KEYDOWN and event.key == K_EQUALS:
				self.zoom += self.zoom * 0.5
				# print self.zoom

			elif event.type == KEYDOWN and event.key == K_MINUS:
				self.zoom -= self.zoom * 0.5
				# print self.zoom

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


			# elif event.type == KEYDOWN and event.key == K_p:
			# 	pygame.image.save(screen, "contact_with_friction.png")
        
	def add(self, thing):  
		
		self.space.add(thing.body, thing.shape)
		
		if thing.type == 'actor':
			self.actors.append(thing)
			
		elif thing.type == 'attractor':
			self.attractors.append(thing)

	def physics(self):
		### Update physics
		# print 'physics'
		for actor in self.actors:
			# print 'actor'
			actor.doResources()
			actor.doModuleEffects(actor.keyStates)

			for attractor in self.attractors:
				if actor.freefalling:
					if actor.orbit == None:
						self.gravitate(actor, attractor)
						actor.orbit = orbmech.initpos_to_orbit_2d(actor.body.position, actor.body.velocity, self.gravitationalConstant * attractor.body.mass)
						print(actor.orbit)
					else:
						self.gravitate(actor, attractor)
				else:
					self.gravitate(actor, attractor)

		self.space.step(self.timestepSize)

	def rotatePolygon(self, points, angle):
		# def Rotate2D(pts,cnt,ang=pi/4):
		return Rotate2D(points,(0,0),angle)
	

	def transformForView(self, position):
		if self.viewpointObject == None:
			return position
		else:
			transformedPosition = position - self.viewpointObject.body.position # offset everything by the position of the viewpointObject, so the viewpoint is at 0,0

			# #rotate everything around the viewPointObject, so that the current largest gravitational influencer is 'down'.
			# #get angle to influencer
			# aPos = self.attractors[0].body.position
			# vPos = self.viewpointObject.body.position
			# angle = math.atan2(aPos[0] - vPos[0], aPos[1] - vPos[1])
			# print angle
			# print aPos
			# print vPos
			# transformedPosition = rotate_point(transformedPosition, angle)
			transformedPosition = transformedPosition * self.zoom  # shrink or expand everything around the 0,0 point

			transformedPosition[0] += 0.5 * self.resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
			transformedPosition[1] += 0.5 * self.resolution[1] # add half the height.

			return transformedPosition


	def drawCircle(self,color, position, radius):
		# transformedPosition = self.transformForView(position)
		pygame.draw.circle(self.screen, color, [int(position[0]), int(position[1])], int((radius * self.zoom)))

	# def drawEllipse(self, orbit):
	def drawActor(self, actor):
	
		rotatedPoints = rotate_polygon(actor.points,actor.body.angle)  # orient the polygon according to the body's current direction in space.
		# print actor.body.angle

		# aPos = self.attractors[0].body.position
		# vPos = self.actors[0].body.position
		# downAngle = math.atan2(aPos[0] - vPos[0], aPos[1] - vPos[1])

		# print downAngle
		transformedPoints = []

		for rotatedPoint in rotatedPoints:
			transformedPoints.append(self.transformForView(rotatedPoint + actor.body.position))


		pygame.draw.lines(self.screen, actor.color, True, transformedPoints)


	def drawModule(self, actor, module):
	
	
		# print actor.body.angle

		# aPos = self.attractors[0].body.position
		# vPos = self.actors[0].body.position
		# downAngle = math.atan2(aPos[0] - vPos[0], aPos[1] - vPos[1])

		# print downAngle
		transformedPoints = []

		for rotatedPoint in module.points: 
			transformedPoints.append(self.transformForView(rotatedPoint + actor.body.position + module.offset))

		# rotatedPoints = rotate_polygon(transformedPoints,actor.body.angle, [self.resolution[0]*0.5,self.resolution[1]*0.5])  # orient the polygon according to the body's current direction in space.
		rotatedPoints = rotate_polygon(transformedPoints,actor.body.angle, self.transformForView(actor.body.position))  # orient the polygon according to the body's current direction in space.


		pygame.draw.lines(self.screen, module.color, True, rotatedPoints)
		# pygame.draw.polygon(self.screen, actor.color, transformedPoints)
		
		if module.active and module.enabled:
			# pygame.draw.circle(self.screen, module.color, 2, 1)
			# def rotate_point(point, angle_rad, center_point=(0, 0)):
			activeCircle = self.transformForView(module.offset + actor.body.position)
			activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
			self.drawCircle(module.color, activeCircle, 2)

		if module.enabled and module.active:
			for giveResource, giveQuantity in module.stores.items(): #module.produces.items():
				if giveResource == 'thrust':
					if actor.keyStates['up']:
						print actor.body.angle
						print module.angle

						# forceAngle = addRadians(actor.body.angle, 0)
						forceAngle = actor.body.angle#addRadians(module.angle, actor.body.angle)
						print forceAngle

						force = [(giveQuantity * math.cos(addRadians(forceAngle, math.pi * 0.5))), giveQuantity * math.sin(addRadians(forceAngle, math.pi * 0.5) )]
						# self.body.apply_impulse_at_local_point(force, module.offset)
						activeCircle = self.transformForView(module.offset + actor.body.position)
						activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
						ananas = (int(activeCircle[0] + force[0] * self.zoom), int(activeCircle[1]+force[1] * self.zoom ) )
						# print ananas
						pygame.draw.lines(self.screen, module.color, True, [activeCircle,ananas])

	# def blitPlanet(self, attractor):
	# 	# circle_surface = pygame.draw.circle(COLOR, RADIUS, WIDTH)
	# 	screen.blit(attractor.image, attractor.position)
        
	def graphics(self):
		### Clear screen
		self.screen.fill(THECOLORS["black"])

		### Draw stuff
		# self.space.debug_draw(self.draw_options)

	
		for attractor in self.attractors:
			if attractor.atmosphere != None:
				# self.drawCircle(attractor.atmosphere.color, attractor.body.position, attractor.radius + attractor.atmosphere.height)
				self.drawActor(attractor.atmosphere)
			# self.drawCircle(attractor.color, attractor.body.position, attractor.radius)
			self.drawActor(attractor)
		for actor in self.actors:
			for module in actor.modules:
				self.drawModule(actor, module)

		# print self.actors[0].availableResources

		i = 0
		for availableResource, availableQuantity in actor.availableResources.items():
			textsurface = self.font.render(str(availableResource) + ': ' + str(availableQuantity), False, (255, 255, 255))
			self.screen.blit(textsurface,(30,i * 30))
			i += 1

		### Flip screen
		pygame.display.flip()
		self.clock.tick(150)
		pygame.display.set_caption("fps: " + str(self.clock.get_fps()))

		

	def step(self):
		self.inputs()
		self.physics()
		self.graphics()

	def start(self):

		newPlanet = Attractor()
		# newButt = Actor((10, -322100), [14000,0])
		newButt = Actor((10, -322100), [50,0])
		twoButt = Actor((1, -320050), [50,0])
		self.add(newButt)
		self.add(twoButt)
		self.player = self.actors[0]
		self.viewpointObject = self.actors[0]
		self.add(newPlanet)

		self.running = True
		while self.running:
			self.step()


mundus = World()
mundus.start()
