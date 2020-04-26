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

class Module():
	def __init__(self):
		self.type = ''
		self.mass = 1

		self.provides = []
		self.consumes = []

		self.active = False

class Actor():
	def __init__(self):
		self.type = 'actor'
		self.mass = 1
		self.radius = 5
		inertia = pymunk.moment_for_circle(self.mass, 0, self.radius, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = 150, 100
		self.body.apply_impulse_at_local_point([150,0], [0,0])
		self.shape = pymunk.Circle(self.body, self.radius, (0,0))
		self.shape.friction = 0.5

		self.orbit = None #initpos_to_orbit(self.)

		self.freefalling = True
		self.interacting = False

		self.modules = []		
		self.availableResources = []
		self.requiredResources = []
	# def enterFreefall(self, attractor):
	# 	self.orbit = func_initpos_to_orbit()
	# 	self.freefalling = True

	# def leaveFreefall(self):
	# 	self.freefalling = False

	def doResources():
		# self.
		for module in self.modules:
			pass


class Attractor():
	def __init__(self):
		self.type = 'attractor'
		self.mass = 100000
		self.radius = 5
		inertia = pymunk.moment_for_circle(self.mass, 0, self.radius, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = 140, 400
		self.shape = pymunk.Circle(self.body, self.radius, (0,0))
		self.shape.friction = 0.5

class World():
	def __init__(self):
		pygame.init()

		self.clock = pygame.time.Clock()
		self.space = pymunk.Space()
		self.space.gravity = (0.0, 0.0)
		self.gravitationalConstant = 1

		self.actors = []
		self.attractors = []

		self.resolution = (600,600)
		self.screen = pygame.display.set_mode(self.resolution)
		self.draw_options = pymunk.pygame_util.DrawOptions(self.screen)
		self.draw_options.flags = self.draw_options.flags ^ pymunk.pygame_util.DrawOptions.DRAW_COLLISION_POINTS 

		self.ch = self.space.add_collision_handler(0, 0)
		self.ch.data["surface"] = self.screen
		# self.ch.post_solve = draw_collision

		self.viewpointObject = None
		self.zoom = 1
		self.pan = [0,0]
		self.rotate = 0

	def gravitate(self, actor, attractor):
		distance = attractor.body.position - actor.body.position # scalar distance between two bodies
		magnitude = mag(distance)
		gravity = self.gravitationalConstant * attractor.body.mass
		appliedGravity = gravity/(magnitude * magnitude)

		components = numpy.divide(distance, magnitude)

		force = components * appliedGravity

		# print force

		rotatedForce = Vec2d(force[0], force[1])
		rotatedForce = rotatedForce.rotated(-actor.body.angle)

		actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])

	
	def inputs(self):
		for event in pygame.event.get():
			if event.type == KEYDOWN and event.key == K_ESCAPE:
				self.running = False
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

		for actor in self.actors:
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


		dt = 0.2/60.0 #1.0/60.0
		self.space.step(dt)


	def transformForView(self, position):
		if self.viewpointObject == None:
			return position
		else:
			transformedPosition = position - self.viewpointObject.body.position

			transformedPosition[0] += 0.5 * self.resolution[0] # add half the width of the screen
			transformedPosition[1] += 0.5 * self.resolution[1] # add half the height.

			return transformedPosition


	def drawCircle(self,color, position, radius):
		transformedPosition = self.transformForView(position)
		pygame.draw.circle(self.screen, color, [int(transformedPosition[0]), int(transformedPosition[1])], radius)

	# def drawEllipse(self, orbit):
		
        
	def graphics(self):
		### Clear screen
		self.screen.fill(THECOLORS["black"])

		### Draw stuff
		# self.space.debug_draw(self.draw_options)

		for actor in self.actors:
			self.drawCircle((120,100,100), actor.body.position, actor.radius)
		for attractor in self.attractors:
			self.drawCircle((0,30,200), attractor.body.position, attractor.radius)


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
		newButt = Actor()
		self.add(newButt)
		self.add(newPlanet)
		self.viewpointObject = self.actors[0]

		self.running = True
		while self.running:
			self.step()


mundus = World()
mundus.start()
