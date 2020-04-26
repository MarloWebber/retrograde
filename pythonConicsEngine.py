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

def rotate_point(point, angle, center_point=(0, 0)):
    """Rotates a point around center_point(origin by default)
    Angle is in degrees.
    Rotation is counter-clockwise
    """
    angle_rad = math.radians(angle % 360)
    # Shift the point so that center_point becomes the origin
    new_point = (point[0] - center_point[0], point[1] - center_point[1])
    new_point = (new_point[0] * math.cos(angle_rad) - new_point[1] * math.sin(angle_rad),
                 new_point[0] * math.sin(angle_rad) + new_point[1] * math.cos(angle_rad))
    # Reverse the shifting we have done
    new_point = (new_point[0] + center_point[0], new_point[1] + center_point[1])
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
	def __init__(self):
		self.height = 100

		self.color = (0,50,200)

		self.seaLevelDrag = 0.1

class Module():
	def __init__(self):
		self.type = ''
		self.mass = 1

		self.provides = []
		self.consumes = []

		self.active = False

class Actor():
	def __init__(self, position, velocity):
		self.type = 'actor'
		self.mass = 1
		self.radius = 5

		size = self.radius
		self.points = [(-size, -size), (-size, size), (size,size), (size, -size)]
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.body.apply_impulse_at_local_point(velocity, [0,0])
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = 0.5

		self.orbit = None #initpos_to_orbit(self.)

		self.freefalling = True
		self.interacting = False

		self.modules = []		
		self.availableResources = []
		self.requiredResources = []


		self.color = (200,50,50)

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
		self.mass = 50000000
		self.radius = 20000
		inertia = pymunk.moment_for_circle(self.mass, 0, self.radius, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = 1, 1
		self.shape = pymunk.Circle(self.body, self.radius, (0,0))
		self.shape.friction = 0.5

		self.atmosphere = Atmosphere()

		self.color = (190,165,145)

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

		pygame.key.set_repeat(50,50) # holding a key down repeats the instruction. https://www.pygame.org/docs/ref/key.html

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
			if event.type == KEYDOWN and event.key == K_RIGHTBRACKET:

				if self.viewpointObject == self.actors[0]:
					self.viewpointObject = self.attractors[0]
				else:
					self.viewpointObject = self.actors[0]

			if event.type == KEYDOWN and event.key == K_EQUALS:
				self.zoom += 0.1
				print self.zoom

			if event.type == KEYDOWN and event.key == K_MINUS:
				self.zoom -= 0.1
				print self.zoom
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
			transformedPosition = position - self.viewpointObject.body.position # offset everything by the position of the viewpointObject, so it is always in the middle of the screen.

			transformedPosition = transformedPosition * self.zoom  # shrink or expand everything around the 0,0 point

			transformedPosition[0] += 0.5 * self.resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
			transformedPosition[1] += 0.5 * self.resolution[1] # add half the height.

			return transformedPosition

	def rotatePolygon(self, points, angle):
		# def Rotate2D(pts,cnt,ang=pi/4):
		return Rotate2D(points,(0,0),angle)
	

	def drawCircle(self,color, position, radius):
		transformedPosition = self.transformForView(position)
		pygame.draw.circle(self.screen, color, [int(transformedPosition[0]), int(transformedPosition[1])], int((radius * self.zoom)))

	# def drawEllipse(self, orbit):
	def drawActor(self, actor):
	
		rotatedPoints = rotate_polygon(actor.points,actor.body.angle) 
		transformedPoints = []

		for rotatedPoint in rotatedPoints:
			transformedPoints.append(self.transformForView(rotatedPoint + actor.body.position))

		pygame.draw.lines(self.screen, actor.color, True, transformedPoints)
		
        
	def graphics(self):
		### Clear screen
		self.screen.fill(THECOLORS["black"])

		### Draw stuff
		# self.space.debug_draw(self.draw_options)

	
		for attractor in self.attractors:
			if attractor.atmosphere != None:
				self.drawCircle(attractor.atmosphere.color, attractor.body.position, attractor.radius + attractor.atmosphere.height)
			self.drawCircle(attractor.color, attractor.body.position, attractor.radius)
		for actor in self.actors:
			self.drawActor(actor)


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
		newButt = Actor((10, -20060), [100,0])
		twoButt = Actor((1, -20050), [50,0])
		self.add(newButt)
		self.add(twoButt)
		self.add(newPlanet)
		self.viewpointObject = self.actors[0]

		self.running = True
		while self.running:
			self.step()


mundus = World()
mundus.start()
