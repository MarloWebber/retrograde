
import sys, copy, time

import numpy, math, random

# pymunk is used for rigid body dynamics.
import pymunk
from pymunk import Vec2d

# astropynamics is used for orbital mechanics calculations.
from AstroPynamics.tools.body import Body as APBody
from AstroPynamics.orbits.orbit import Orbit
from astropy.time import Time, TimeDelta
import time
    
# global contact
# global shape_to_remove

# pyglet is used primarily for graphics but also handles keyboard input.
import pyglet
from pyglet.gl import *
from pyglet.window import Window
from pyglet.window import key
from pyglet.window import mouse


# cProfile is used along with snakeviz to analyze the game's performance.
import cProfile

# dill is used to save and load information so you can save ships and your game progress.
import pickle
import dill

gravitationalConstant = 0.01

def area_of_annulus(inner, outer):
	return math.pi * ((inner**2) - (outer**2))


def make_circle(radius, points, noise):
	# returns a list of points defining a circle.
	circle = []
	randomnesses = []
	angleStep = 2*math.pi / points

	for i in range(0,points):
		randomnesses.append((random.randint(0,noise)  ))

	# for i in range(0,points):
	# 	if randomnesses[i] > 0:
	# 		randomnesses[i] = randomnesses[i] **2
	# 	else:
	# 		randomnesses[i] = (abs(randomnesses[i]) **2 )* -1

	for i in range(0,points):
		circle.append([((randomnesses[i] + radius)* math.cos(i * angleStep) ),(( randomnesses[i]+ radius)* math.sin(i * angleStep) )  ])
	return circle

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


def mass_of_a_sphere(density, radius):
	return density * 4/3 * math.pi * radius**3

def mass_of_a_circle(density, radius):
	return density * 2 * math.pi * radius**2


def addRadians(a, b):
	return (a + b + math.pi) % (2*math.pi) - math.pi

def mirror_polygon(polygon):
	transformedPoints = []
	for point in polygon:
		transformedPoints.append( [point[0] * -1, point[1]] )
	return transformedPoints


def mag(x):
	return numpy.sqrt(x.dot(x))

def differenceBetweenAngles(a, b):
	return math.pi - abs(abs( a - b) - math.pi); 



class Atmosphere():
	def __init__(self, radius, height, bottomDensity, topDensity, color, outerColor):

		# if atmosphereType is "earthAtmosphere":
		self.height = height
		self.bottomDensity = bottomDensity
		self.topDensity = topDensity
		self.color = color
		self.outerColor = outerColor

		# each atmosphere layer is an annulus.
		# it is rendered as a triangle strip with a gradient between the inner and outer edges.
		self.n_points = 100
		self.innerPoints = make_circle(radius, self.n_points, 0)
		self.outerPoints = make_circle(radius+self.height, self.n_points, 0)
		self.mass = ( (self.bottomDensity + self.topDensity) / 2 ) * area_of_annulus(radius+self.height, radius)

class Illuminator():
	def __init__(self, offset, radius, color, brightness):
		self.radius = 1000
		self.color = color
		self.intensity = 1
		self.offset = offset #(100, -320030)
		self.position = offset # the offset is the illuminators position relative to the parent module. the 'position' field is just computed once per turn, so you don't have to do it once per point.
		self.transformedPosition = [0,0]
		self.isItLightingThisPoint = False
		self.distanceScalar = 0
		self.distanceVector = [0,0]
		self.brightness = brightness

class Maneuver():
	# description of an AI behaviour item.
	def __init__(self, maneuverType, parameter1=None, parameter2=None, parameter3=None):
		self.maneuverType = maneuverType
		self.parameter1 = parameter1 # parameters can be used to mean different things depending on what kind of maneuver it is.
		self.parameter2 = parameter2
		self.parameter3 = parameter3
		self.completed = False
		self.event1 = False
		self.event2 = False
		self.event3 = False
		self.event4 = False

	def perform(self, actor):
		if not self.completed:

			#http://www.braeunig.us/space/orbmech.htm orbital mechanics

			if self.maneuverType == 'attack target':
				pass 

			if self.maneuverType == 'lead target':
				# fly directly into the target, accelerating the whole way. better if you can lead the target.
				if not self.event1:
					print('lead target')
					self.event1 = True 


				if self.parameter1.orbit is not None: # the target is on a rail
					speedAngle = math.atan2(actor.body.velocity[1], actor.body.velocity[0] )
					angleToTarget = math.atan2(self.parameter1.body.position[1]- actor.body.position[1] , self.parameter1.body.position[0]- actor.body.position[0])
					distanceToTarget = mag(actor.body.position - self.parameter1.body.position)
					speedMagnitude = mag(actor.body.velocity)
					timeToTarget = distanceToTarget / speedMagnitude
					futuretAn = self.parameter1.orbit.tAnAtTime(timeToTarget)
					targetFuturePosition = self.parameter1.orbit.cartesianCoordinates(futuretAn)
					color_point[0] = targetFuturePosition[0]
					color_point[1] = targetFuturePosition[1]
					targetLeadAngle = math.atan2(targetFuturePosition[1]- actor.body.position[1] , targetFuturePosition[0]- actor.body.position[0])
					actor.setPoint = targetLeadAngle + 0.5 * math.pi
					


			if self.maneuverType == 'fast intercept':
				# point at the target, fly half way there at full speed. turn around at the halfway point and decelerate the rest of the way.
				if not self.event1:
					print('fast intercept')
					self.event1 = True 


			if self.maneuverType == 'add speed':
				# points the ship in a direction and accelerates until it is comoving with the vector parameter1.
				if not self.event1:
					print('add speed')
					self.event1 = True
					if self.parameter2 is None:
						self.parameter2 = mag(self.parameter1) * 0.1

				speedDifference = self.parameter1 - actor.body.velocity
				resultAngle = math.atan2(speedDifference[1], speedDifference[0])

				actor.setPoint = resultAngle
				if (angleDifference < 0.1 and angleDifference > 0) or angleDifference > (2*math.pi - 0.1): # only fire engines if the ship is pointing vaguely in the right direction
					actor.keyStates['up'] = True
				else:
					actor.keyStates['up'] = False

				if mag(speedDifference) < self.parameter2:
						actor.keyStates['up'] = False
						self.completed = True


			if self.maneuverType == 'change periapsis':
				# event1 prints 'change periapsis'
				# event2 is reaching the apoapsis and starting the burn
				# parameter1 is the new periapsis

				if not self.event1:
					print('change periapsis...')
					self.event1 = True

					if self.parameter2 is None: # if the angle for the new periapsis was passed in as none, you should burn from the apoapsis.
						self.parameter2 = 0

				actor.setPoint = actor.prograde + 0.5 * math.pi
				if actor.orbit is not None:
					if actor.orbit.tAn >  self.parameter2 - 0.1 and actor.orbit.tAn <  self.parameter2 + 0.1: # the ship is near the apoapsis
						self.event2 = True
						print('change periapsis - mark')

					if self.event2:
						angleDifference = addRadians(  actor.setPoint, - actor.body.angle)
						if (angleDifference < 0.1 and angleDifference > 0) or angleDifference > (2*math.pi - 0.1): # only fire engines if the ship is pointing vaguely in the right direction
							actor.keyStates['up'] = True
						else:
							actor.keyStates['up'] = False
					if actor.orbit.getApoapsis() > self.parameter1:
						actor.keyStates['up'] = False
						self.completed = True

			if self.maneuverType == 'change apoapsis':
				if not self.event1:
					print('change apoapsis')
					self.event1 = True
				actor.setPoint = actor.retrograde + 0.5 * math.pi
				if actor.orbit is not None:
					if actor.orbit.tAn > 0 and actor.orbit.tAn < 0.1: # the ship is near the apoapsis
						self.event2 = True
					if self.event2:
						angleDifference = actor.setPoint - actor.body.angle
						if (angleDifference < 0.1 and angleDifference > 0) or angleDifference > (2*math.pi - 0.1): # only fire engines if the ship is pointing vaguely in the right direction
							actor.keyStates['up'] = True
						else:
							actor.keyStates['up'] = False
					if actor.orbit.getApoapsis() < self.parameter1:
						actor.keyStates['up'] = False
						self.completed = True

			if self.maneuverType == 'transfer':

				# parameter1 is the object we're leaving
				# parameter2 is the one we're heading toward

				
				#event1 is announcement
				if not self.event1:
					print('transfer')
					self.event1 = True

					actor.setPoint = actor.prograde + math.pi * 0.5

					# 
				

				# event2 is arrival at burn point to leave orbit.
				if self.event1 and not self.event2 and not self.event3:
					positionDifference = [ self.parameter1.body.position[0]- self.parameter2.body.position[0], self.parameter1.body.position[1]- self.parameter2.body.position[1]]
					ejectionPointAngle = math.atan2(positionDifference[1] , positionDifference[0]) 
					# if ejectionPointAngle > 2 * math.pi:
					# 	ejectionPointAngle -= 2 * math.pi
					ejectionPointAngle = addRadians(actor.orbit.aPe, ejectionPointAngle)
					ejectionPointAngle = addRadians(1.5 * math.pi, ejectionPointAngle)

					print('ejectionPointAngle: '+str(ejectionPointAngle))
					print('addRadians(actor.orbit.tAn ,actor.orbit.aPe): '+str(addRadians(actor.orbit.tAn ,actor.orbit.aPe)))
					# print(addRadians(actor.orbit.tAn ,actor.orbit.aPe))


					actor.setPoint = actor.prograde + math.pi * 0.5

					if (actor.orbit.tAn > ejectionPointAngle -0.1 and actor.orbit.tAn < ejectionPointAngle + 0.1):
						print('reached burn point')
						self.event2 = True

				if self.event2 and not self.event3:
					actor.setPoint = actor.prograde + math.pi * 0.5

					positionDifference = [ self.parameter1.body.position[0]- self.parameter2.body.position[0], self.parameter1.body.position[1]- self.parameter2.body.position[1]]
					desiredPeriapsis = mag(numpy.array(positionDifference)) + self.parameter2.radius
					if self.parameter2.atmosphere is not None:
						desiredPeriapsis += self.parameter2.atmosphere.height

					angleDifference = differenceBetweenAngles(actor.setPoint,actor.body.angle)
					if (angleDifference < 0.5 and angleDifference > -0.5): # only fire engines if the ship is pointing vaguely in the right direction
						actor.keyStates['up'] = True
					else:
						actor.keyStates['up'] = False

					print(actor.orbit.getApoapsis())
					print(desiredPeriapsis)
					if actor.orbit.getApoapsis() > desiredPeriapsis:
						self.event3 = True
						print('raised trajectory to target, ready to cruise')
						actor.keyStates['up'] = False

				#event3 is periapsis reaching desired altitude above target
				# if self.event3: 
					# actor.setPoint = actor.prograde + math.pi * 0.5
					



				# event4 is entering target SoI and adjusting orbit to not miss or crash
				# if actor.orbit is not None:
				# 	print(actor.orbit.tAn)
					# if (differenceBetweenAngles(actor.orbit.tAn,math.pi) < 0.1) and differenceBetweenAngles(actor.orbit.tAn,math.pi) > -0.1 :



				if(actor.orbiting is self.parameter2) and self.event2 and self.event3 and not self.event4:
					if actor.orbit is not None:

						if actor.orbit.isCrashing() or actor.orbit.getPeriapsis() < self.parameter2.radius or actor.orbit.getApoapsis() < self.parameter2.radius :
							print('crashing. need to raise orbit')
							actor.setPoint = actor.nadir + (math.pi * 0.5) +( math.pi * 0.5)
							angleDifference = actor.setPoint - actor.body.angle
							if (angleDifference < 0.3 and angleDifference > 0) or angleDifference > (2*math.pi - 0.3): # only fire engines if the ship is pointing vaguely in the right direction
								actor.keyStates['up'] = True
							else:
								actor.keyStates['up'] = False
						

						else:

							# if you're orbiting the right thing
							if actor.orbiting is self.parameter2:

								if  actor.orbit.getPeriapsis() < 2 * self.parameter2.radius:
									print('adjusting periapsis')
									actor.setPoint = actor.prograde +( math.pi * 0.5)
									angleDifference = differenceBetweenAngles( actor.setPoint, actor.body.angle)
									if angleDifference < 0.1 and angleDifference > -0.1: # only fire engines if the ship is pointing vaguely in the right direction
										actor.keyStates['up'] = True
									else:
										actor.keyStates['up'] = False
							


								else:

									#event5 is circularization around target
									print('arrived at target. adding circularization to queue.')

									actor.keyStates['up'] = False
									self.completed = True
									actor.maneuverQueue.remove(self)
									actor.maneuverQueue.append(Maneuver('circularize', 2 * self.parameter2.radius))

									print(actor.maneuverQueue)

							# you're not even orbiting the right thing, try to get there
							if actor.orbiting is not self.parameter2 or actor.orbiting is None:

								angleToTarget = math.atan2(self.parameter2.body.position[1] - actor.body.position[1], self.parameter2.body.position[0] - actor.body.position[0])
								actor.setPoint = angleToTarget +( math.pi * 0.5)
								angleDifference = differenceBetweenAngles( actor.setPoint, actor.body.angle)
								if angleDifference < 0.1 and angleDifference > -0.1: # only fire engines if the ship is pointing vaguely in the right direction
									actor.keyStates['up'] = True
								else:
									actor.keyStates['up'] = False
						




					else:
						print('decelerating out of hyperbolic orbit')
						velocityVector = math.atan2(actor.body.velocity[1], actor.body.velocity[0])
						decelVector = addRadians(velocityVector, math.pi)
						actor.setPoint = decelVector + math.pi * 0.5
						angleDifference = differenceBetweenAngles( actor.setPoint, actor.body.angle)
						print(actor.setPoint)
						print(angleDifference)
						if angleDifference < 0.1: # only fire engines if the ship is pointing vaguely in the right direction
							actor.keyStates['up'] = True
						else:
							actor.keyStates['up'] = False


				

				


			if self.maneuverType == 'circularize':
				# true to circ at periapsis. false to circ at apoapsis.
				if self.parameter1:
					actor.setPoint = actor.prograde + 0.5 * math.pi
				else:
					actor.setPoint = actor.retrograde + 0.5 * math.pi

				if actor.orbit is not None:
					if (actor.orbit.tAn > 0 and actor.orbit.tAn < 0.1 and self.parameter1) or (actor.orbit.tAn > math.pi -0.1 and actor.orbit.tAn < math.pi + 0.1 and not self.parameter1): # the ship is near the apoapsis
						self.event2 = True
					if self.event2:
						angleDifference = actor.setPoint - actor.body.angle
						if (angleDifference < 0.1 and angleDifference > 0) or angleDifference > (2*math.pi - 0.1): # only fire engines if the ship is pointing vaguely in the right direction
							actor.keyStates['up'] = True
						else:
							actor.keyStates['up'] = False
					if actor.orbit.getApoapsis() < actor.orbit.getPeriapsis() + (actor.orbit.getPeriapsis()* 0.01):
						actor.keyStates['up'] = False
						self.completed = True

			if self.maneuverType == 'match velocity':
				pass

			if self.maneuverType == 'rendezvous':
				# pass
				# the ship guides itself to a target which is orbiting the same attractor.
				# it is safer to always phase the target by going upwards, because then you never have to worry about hitting the planet.

				# figure out the time discrepancy between you and your target.

				if not self.event1:
					print('rendezvous')
					self.event1 = True
				elif not self.event2:
					pass


			if self.maneuverType == 'takeoff':
				#In takeoff, parameter1 is the orbit height to achieve, and parameter2 is the attractor you are taking off from.
					#event1 is placing the periapsis high enough, and event2 is reaching the periapsis and being ready to circularize.
				

				actorHeightFromAttractorCenter = mag(actor.body.position - self.parameter2.body.position)
				naturalHeight = ((actorHeightFromAttractorCenter - self.parameter2.radius) / self.parameter1) # this is a number between 0 and 1 which is signifies the actors depth into this atmosphere layer.
				
				# Go straight up at the start, and turn sideways once you get up a bit.
				if not self.event1:
					if naturalHeight < 0.5:
						actor.setPoint = actor.zenith + 0.5 * math.pi
						actor.keyStates['up'] = True
					else:
						actor.setPoint = (actor.zenith + 0.5 * math.pi) +( naturalHeight * 0.5 * math.pi )
						actor.keyStates['up'] = True

					if actor.orbit is not None:
						if actor.orbit.getPeriapsis() > self.parameter1 + self.parameter2.radius:
							print('the apoapsis is high enough')
							self.event1 = True
							actor.keyStates['up'] = False
					
				else:
					if actor.orbit is not None:

						# first, focus on raising the apoapsis
							actor.setPoint = actor.prograde + 0.5* math.pi #- 0.5* math.pi # why the heck is it off by so much. literally should just be prograde.
							if self.event2:
								if actor.orbit.tAn < 0.2 and actor.orbit.getPeriapsis() > (self.parameter1  + self.parameter2.radius) and actor.orbit.getApoapsis()  > (self.parameter1 + self.parameter2.radius):
									actor.keyStates['up'] = False
									self.completed = True
									print('maneuver completed')
							else:
								if( actor.orbit.tAn < math.pi + 0.05 and actor.orbit.tAn > math.pi - 0.05 )or self.event2:
									print('reached the periapsis, ready to circularize')
									self.event2 = True
									if abs(actor.setPoint) - abs(actor.body.angle) < 0.1: # only fire engines if the ship is pointing vaguely in the right direction
										actor.keyStates['up'] = True
									else:
										actor.keyStates['up'] = False

class Actor():
	def __init__(self, name, modulesList, position, velocity, angle, isPlayer=False):
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
			modulePoints = []
			for point in module.points:
				transformedPoint = [point[0] + module.offset[0], point[1] + module.offset[1]]
				transformedPoint = rotate_point(transformedPoint, module.angle, module.offset)
				modulePoints.append(transformedPoint)
			self.points += modulePoints
			for resource, quantity in list(module.initialStores.items()):
				if resource not in self.storagePool: self.storagePool[resource] = quantity
				else: self.storagePool[resource] += quantity
			for resource, quantity in list(module.initialStores.items()):
				if resource not in self.availableResources: self.availableResources[resource] = quantity
				else: self.availableResources[resource] += quantity
			for resource, quantity in list(module.stores.items()):
				if resource not in self.maximumStores: self.maximumStores[resource] = quantity
				else: self.maximumStores[resource] += quantity
		
		self.inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, self.inertia)
		self.body.position = position
		self.body.velocity = velocity
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = 0.9
		self.orbit = None 
		self.freefalling = True
		self.color = (200,50,50)
		self.keyStates = {
			'up': False,
			'down': False,
			'left': False,
			'right': False,
			'Fire': False,
			'face direction': None,
			'hyperdrive engage':False,
			'strafe forward':False,
			'strafe back':False,
			'strafe left':False,
			'strafe right':False
		}
		self.orbiting = None
		self.stepsToFreefall = 1 # After an object has collided or accelerated, it doesn't return to freefall straight away. It waits 1 or 0 turns. 
								# If you wait 0 turns, an orbit will be prepared for you, so you can view your trajectory while still burning your rockets.
								# if you wait 1 turn, no orbit will be prepared or drawn, which is good for objects that are sitting on the ground.
								# objects that are sitting on the ground will also not gravitate, to prevent jiggling. This is ended when something hits them or they accelerate.

		self.decompEnergy = 5000000
		self.desiredAngle = 0
		self.exemptFromGravity = False
		self.body.angle = angle
		self.orbitPoints = []

		self.setPoint = 0
		self.prograde = 0
		self.retrograde = 0
		self.nadir = 0
		self.zenith = 0

		self.maneuverQueue = []
		self.combatantType = 'defender' # defenders will shoot at you while still doing what they're doing. attackers will pursue you. missiles will pursue you with the intent to ram.
		self.autoPilotActive = False

		self.autoPilotGoals = [] # a list of goals which the ai will use sequences of maneuvers to accomplish.

		self.target = None # another actor that this one can lock with radar and scanners.
		self.selectedWeapon = None
		self.hyperdriveDestination = None
		self.jumping = False

	def leaveFreefall(self, stepsToFreefall=1):
		self.stepsToFreefall = stepsToFreefall
		self.freefalling = False
		self.orbit = None
		
	def doResources(self, timestepSize):
		# - tally the amount of stored resources. first, zero everything out
		for module in self.modules:
			for resource, quantity in list(module.resources.items()):
					self.availableResources[resource] = 0

		# add all the unstored resources being produced by active modules
		for module in self.modules:
			for resource, quantity in list(module.resources.items()):
				# 'instantaneous' resources are different to material quantities that are produced and stored. 
				if resource is 'torque' or resource is 'thrust': 
					adjustedQuantity = quantity
				else:
					adjustedQuantity = quantity * timestepSize
				if adjustedQuantity > 0 and module.enabled and module.active:
					if resource not in self.availableResources:
						self.availableResources[resource] = adjustedQuantity
					else:
						self.availableResources[resource] += adjustedQuantity

		# - turn modules on and off
		for module in self.modules:
			module.enabled = True
			for resource, quantity in list(module.quiescent.items()):
				if resource is 'torque' or resource is 'thrust': 
					adjustedQuantity = quantity
				else:
					adjustedQuantity = quantity * timestepSize
				if adjustedQuantity < 0:
					if resource in self.storagePool:
						availableAmount = self.storagePool[resource] + self.availableResources[resource]
					else:
						availableAmount = self.availableResources[resource]
					if availableAmount < abs(adjustedQuantity):
						module.enabled = False

		# - consume and produce resources
		for module in self.modules:
			if module.enabled and module.active:
				for resource, quantity in list(module.resources.items()):
					if resource is 'torque' or resource is 'thrust': 
						adjustedQuantity = quantity
					else:
						adjustedQuantity = quantity * timestepSize
					if adjustedQuantity > 0: # producing resource
						if resource in self.storagePool:
							remainingCapacity = self.maximumStores[resource] - self.storagePool[resource]
							if adjustedQuantity > remainingCapacity:
								self.storagePool[resource] += remainingCapacity
								self.availableResources[resource] += adjustedQuantity - remainingCapacity
							else:
								self.storagePool[resource] += adjustedQuantity

					else: # consuming resource
						self.availableResources[resource] += adjustedQuantity # adding a negative number is a subtraction
						if self.availableResources[resource] < 0:
							if resource in self.storagePool:
								self.storagePool[resource] += self.availableResources[resource]
								self.availableResources[resource] = 0
							else:
								module.enabled = False
								module.active = False

		for resource, quantity in list(self.storagePool.items()):
			if quantity > self.maximumStores[resource]: quantity = self.maximumStores[resource]			

	def doModuleEffects(self, keyStates, timestepSize):
		# doModuleEffects performs engine thrust and RCS torque. It reports whether or not the craft has been accelerated.

		ifThrustHasBeenApplied = False
		for module in self.modules:
			if module.enabled:
				for giveResource, giveQuantity in list(module.resources.items()):
					if giveResource == 'thrust':

							# thruster control is provided by engines, fuckin deal with it
							if self.keyStates['strafe forward']:
								# print('strafe forward')
								force = [0 , giveQuantity * timestepSize * 5  ]
								self.body.apply_impulse_at_local_point(force, (0,0))
								ifThrustHasBeenApplied = True
							if self.keyStates['strafe back']:
								force = [0 , -(giveQuantity * timestepSize * 5 )]
								self.body.apply_impulse_at_local_point(force, (0,0))
								ifThrustHasBeenApplied = True
							if self.keyStates['strafe left']:
								# print('strafe left')
								force = [-(giveQuantity * timestepSize * 5 ) , 0]
								self.body.apply_impulse_at_local_point(force, (0,0))
								ifThrustHasBeenApplied = True
							if self.keyStates['strafe right']:
								force = [giveQuantity * timestepSize * 5 , 0]
								self.body.apply_impulse_at_local_point(force, (0,0))
								ifThrustHasBeenApplied = True


							if module.active:
								force = [(giveQuantity * timestepSize * 500 * math.cos(addRadians(module.angle, math.pi * 0.5))), -giveQuantity * timestepSize * 500 * math.sin(addRadians(module.angle, math.pi * 0.5) )]
								self.body.apply_impulse_at_local_point(force, (0,0))
								ifThrustHasBeenApplied = True
					elif giveResource == 'torque':
							self.setPoint = self.setPoint % (2*math.pi)
							self.body.angle = self.body.angle % (2*math.pi)
							correctionDirection = self.setPoint - self.body.angle
							if abs(correctionDirection) > 0.001:
								module.active = True
								if correctionDirection > math.pi or correctionDirection < -math.pi:
									correctionDirection = -correctionDirection
								if correctionDirection == 0:
									sign = -1
								else:
									sign = -correctionDirection / abs(correctionDirection)
								torqueAmount = sign * giveQuantity * timestepSize * 100

								# apply two impulses, pushing in opposite directions, an equal distance from the center to create torque
								if timestepSize  * 3 * 100 < 50:
									self.body.apply_impulse_at_local_point([-torqueAmount,0], [0,-module.momentArm])
									self.body.apply_impulse_at_local_point([torqueAmount,0], [0,module.momentArm])
								# else:
									# self.body.angle = self.setPoint
									# self.body._set_angular_velocity(0)
							else:
								module.active = False

					elif giveResource == 'warp energy':
						if self.keyStates['hyperdrive engage'] and module.enabled:
							self.storagePool['warp energy'] += giveQuantity * timestepSize
							if self.storagePool['warp energy'] >= module.stores['warp energy']:
								self.jumping = True
								self.storagePool['warp energy'] = 0

							module.active = True
						else:
							module.active = False

		return ifThrustHasBeenApplied

	def flightComputer(self):
		# this function describes the AI flight behaviour and player autopilot

		# allow the player to hold attitude
		if self.keyStates['face direction'] is not None:
			if self.keyStates['face direction'] == 'retrograde': self.setPoint = self.retrograde + 0.5 * math.pi
			if self.keyStates['face direction'] == 'prograde': self.setPoint = self.prograde +  0.5 * math.pi
			if self.keyStates['face direction'] == 'nadir': self.setPoint = self.nadir +  0.5 * math.pi
			if self.keyStates['face direction'] == 'zenith': self.setPoint = self.zenith +  0.5 * math.pi
			if self.target is not None:
				if self.keyStates['face direction'] == 'target': self.setPoint = math.atan2(self.target.body.position[1] - self.body.position[1],(self.target.body.position[0] - self.body.position[0])) + 0.5 * math.pi
			if str.isnumeric(self.keyStates['face direction']): self.setpoint = float(self.keyStates['face direction'])

		# perform autopilot maneuvers for the player and for NPCs
		if self.autoPilotActive:
			if len(self.maneuverQueue) > 0:
				self.maneuverQueue[0].perform(self)


class FakeBody():
	def __init__(self, position):
		self.position = position

class Cloud():
	def __init__(self, points, color, outlineColor, position):
		self.points = points
		self.color = color
		self.outlineColor = outlineColor
		self.body = FakeBody(position)


def make_clouds(cloudlineRadius, attractorPosition, fluffyness, n_points, atmosphereDepth, planetRadius):


	n_points = n_points
	fluffyness = fluffyness
	# cloudlineRadius = 100000

	# returns an array of cloud polygons that are wrapped around a planet.
	innerLine = []
	outerline = []

	innerRandomWalk = 0
	outerRandomWalk = 0

	for i in range(0,n_points):
		innerRandomWalk += random.randint(-fluffyness,fluffyness)
		outerRandomWalk += random.randint(-fluffyness,fluffyness)

		# if innerRandomWalk+cloudlineRadius> planetRadius + (atmosphereDepth * 0.9):
		# 	innerRandomWalk -=random.randint(0,fluffyness)
		# elif innerRandomWalk> cloudlineRadius + atmosphereDepth * 0.5:
		# 	innerRandomWalk -=random.randint(0,0.5*fluffyness)

		# if innerRandomWalk<planetRadius:
		# 	innerRandomWalk +=random.randint(0,2*fluffyness)


		# if outerRandomWalk+cloudlineRadius> planetRadius + (atmosphereDepth * 0.9):
		# 	outerRandomWalk -=random.randint(0,fluffyness)
		# elif outerRandomWalk> cloudlineRadius + atmosphereDepth * 0.5:
		# 	outerRandomWalk -=random.randint(0,0.5*fluffyness)
		# if outerRandomWalk<planetRadius:
		# 	outerRandomWalk +=random.randint(0,2*fluffyness)

		innerLine.append(innerRandomWalk)
		outerline.append(outerRandomWalk)

	clouds = []
	workingCloud = []
	workingOnACloud = False
	workingCloudStartIndex = 0

	sliceSize = (2 * math.pi) / n_points

	for i in range(1,n_points-1):
		# print(i)
		if outerline[i-1] < innerLine[i-1] and outerline[i] > innerLine[i] and not workingOnACloud:
			# the beginning of a cloud
			# print('the beginngin of a cldoueg')
			workingOnACloud = True
			workingCloudStartIndex = i

		if workingOnACloud:
			workingCloud.append([ (outerline[i] + cloudlineRadius * math.cos(i * sliceSize)),( outerline[i] + cloudlineRadius * math.sin(i * sliceSize))])

		if outerline[i-1] > innerLine[i-1] and outerline[i] < innerLine[i]:
			# the end of a cloud
			# print('the genege of a cllosuef')
			workingOnACloud = False

			diffulence = i - workingCloudStartIndex

			for k in range(0, diffulence):
				j = i-k
				# print(j)
				workingCloud.append( [( innerLine[j] + cloudlineRadius * math.cos(j * sliceSize)),( innerLine[j] + cloudlineRadius * math.sin(j * sliceSize))])
			clouds.append( Cloud( workingCloud, [200,200,200,255], [250,250,250,255],attractorPosition ) )
			workingCloud = []
			workingCloudStartIndex = i
			i-=1

		if workingOnACloud and i >= workingCloudStartIndex + 50:
			# the end of a cloud
			# print('the genege of a cllosuef')
			# workingOnACloud = False

			diffulence = i - workingCloudStartIndex


			for k in range(0, diffulence):
				j = i-k
				# print(j)
				workingCloud.append( [( innerLine[j] + cloudlineRadius * math.cos(j * sliceSize)),( innerLine[j] + cloudlineRadius * math.sin(j * sliceSize))])
			clouds.append( Cloud( workingCloud, [200,200,200,255], [250,250,250,255],attractorPosition ) )
			workingCloud = []
			workingCloudStartIndex = i
			i-=1

	return clouds

			

	# print(clouds)
	# return Cloudset(clouds, [200,200,200,255], [250,250,250,255],attractorPosition )






class Attractor():
	def __init__(self,planetName,radius,density,friction,color,outlineColor,atmosphere, position, clouds=None, noise=0):
		self.planetName = planetName
		self.radius = radius
		self.density = density
		self.friction = friction
		self.color = color
		self.outlineColor = outlineColor
		self.atmosphere = atmosphere

		# create pymunk physical body and shape
		self.mass = mass_of_a_sphere(self.density, self.radius)
		size = self.radius
		self.points = make_circle(self.radius, 120, noise)
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.shape = pymunk.Poly(self.body, self.points)
		self.shape.friction = self.friction

		# create astropynamics orbit-able body
		self.APBody = APBody(self.planetName, self.mass * gravitationalConstant * 0.163, self.radius)

		self.clouds = clouds

class buildMenuItem():
	# a buildMenuItem is literally the tiles in the build menu you can click and drag to add modules to your ship.
	def __init__(self, module, boundingRectangle=((0,0),(1,1))):
		self.module = module
		self.quantity = 1
		self.boundingRectangle = boundingRectangle


class SolarSystem():
	# fills the world with a variety of preset planets and characters.
	def __init__(self, solarSystemName, contents, position, color, outlineColor, links, hyperspaceThreshold):
		self.contents = contents
		self.solarSystemName = solarSystemName
		self.position = position
		self.color = color
		self.outlineColor = outlineColor
		self.links = links
		self.hyperspaceThreshold = hyperspaceThreshold


class BackgroundStar():
	def __init__(self, position, color, size):
		self.position = position
		self.color = color
		self.size = size



# # shipyard
# dinghy = [Module('generator',[0,0]), Module('engine 10',[0,15]), Module('RCS',[0,-10]) ]
# lothar = [Module('generator',[0,0]), Module('engine 10',[-13,8], 0.6/math.pi), Module('engine 10',[13,8],-0.6/math.pi), Module('RCS',[-13,-10]), Module('RCS',[13,-10]) , Module('cannon 10',[0,-10]) ]
# boldang = [Module('spar 10',[0,-100], (0.5* math.pi)), Module('box 10',[0,0])]
# bigmolly = [Module('box 100',[0,0]), Module('spar 100',[1000,0], 0.5 * math.pi),Module('box 100',[-1000,0]),Module('box 100',[2000,0]), Module('box 100',[-2000,0]),  Module('box 100',[3000,0])]
# derelict_hyperunit = [Module('hyperdrive 10',[0,0])]
# ida_frigate = [Module('generator',[0,50]), Module('engine 10',[-13,58], 0.6/math.pi), Module('engine 10',[13,58],-0.6/math.pi), Module('RCS',[-13,40]), Module('RCS',[13,40]) , Module('box 10',[0,10]), Module('box 10',[0,-40]), Module('hyperdrive 10',[0,-75]), Module('RCS',[-13,-70]), Module('RCS',[13,-70]) ]


# ida_frigate_instance = Actor('NPC lothar', ida_frigate,(1, 121000), [17000,0], 0.6 * math.pi, True)