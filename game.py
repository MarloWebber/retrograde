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

#
# import galaxy.getSolarSystems
# from galaxy import *
# from retrograde import Module
# from retrograde import Actor
# from retrograde import Attractor
# from retrograde import Atmosphere
# from retrograde import ModuleEffect
# from retrograde import Illuminator
# from retrograde import SolarSystem
# from retrograde import gravitationalConstant
# from retrograde import BackgroundStar
from retrograde import *
from ships import *
from galaxy import *

absolute_resolution = (1600,900) 	# stays the same even when the hud changes the view size
resolution = (1600,900)			 	# current resolution for the view of the in-game world
resolution_half = (1600/2,900/2)	
xCompensation = 200					# how big the HUD menus on the left and right are.

topLimit = [0,0]
bottomLimit = [0,0]
window = pyglet.window.Window(width=1600, height=900, fullscreen=True)
# window = pyglet.window.Window(width=1600, height=900)
label = pyglet.text.Label('Abc', x=5, y=5)

white = [255]*4 

color_point = [0,0]


previousHUDlistQuantities = [None]*100
previousHUDstrings = [None]*100
previousHUDlabels = [None]*100

previousHUDQuantityLabels = [None]*100

def bv2rgb(bv):
  if bv < -0.4: bv = -0.4
  if bv > 2.0: bv = 2.0
  if bv >= -0.40 and bv < 0.00:
    t = (bv + 0.40) / (0.00 + 0.40)
    r = 0.61 + 0.11 * t + 0.1 * t * t
    g = 0.70 + 0.07 * t + 0.1 * t * t
    b = 1.0
  elif bv >= 0.00 and bv < 0.40:
    t = (bv - 0.00) / (0.40 - 0.00)
    r = 0.83 + (0.17 * t)
    g = 0.87 + (0.11 * t)
    b = 1.0
  elif bv >= 0.40 and bv < 1.60:
    t = (bv - 0.40) / (1.60 - 0.40)
    r = 1.0
    g = 0.98 - 0.16 * t
  else:
    t = (bv - 1.60) / (2.00 - 1.60)
    r = 1.0
    g = 0.82 - 0.5 * t * t
  if bv >= 0.40 and bv < 1.50:
    t = (bv - 0.40) / (1.50 - 0.40)
    b = 1.00 - 0.47 * t + 0.1 * t * t
  elif bv >= 1.50 and bv < 1.951:
    t = (bv - 1.50) / (1.94 - 1.50)
    b = 0.63 - 0.6 * t * t
  else:
    b = 0.0
  return (r, g, b)

def sign(x):
	return x * abs(x)

def round_to_n(x, n):
	return round(x, -int(math.floor(math.log10(x))) + (n - 1))

def getFarthestPointInPolygon(polygon):
	farthestPoint = [0,0]
	for point in polygon:
		distance = mag(numpy.array(point))
		if distance > mag(numpy.array(farthestPoint)):
			farthestPoint = point
	return farthestPoint

def boundPolygon(polygon): # returns a rectangle that is a bounding box around the polygon.
	mostX  = 0
	leastX = polygon[0][0]
	mostY  = 0
	leastY = polygon[0][1]
	for point in polygon:
		if point[0] < leastX: leastX = point[0]
		if point[1] < leastY: leastY = point[1]
		if point[0] > mostX: mostX = point[0]
		if point[1] > mostY: mostY = point[1]
	return [[leastX, leastY], [mostX, mostY]]

def pointInRect(point,rect):
    if (rect[0][0] < point[0] and point[0] < rect[1][0]):
        if (rect[0][1] < point[1] and point[1] < rect[1][1]):
            return True
    return False



# def rotate_point(point, angle_rad, center_point=(0, 0)):
#     """Rotates a point around center_point(origin by default)
#     Angle is in degrees.
#     Rotation is counter-clockwise
#     """
#     new_point = (point[0] - center_point[0], point[1] - center_point[1])  # Shift the point so that center_point becomes the origin
#     new_point = (new_point[0] * math.cos(angle_rad) - new_point[1] * math.sin(angle_rad),
#                  new_point[0] * math.sin(angle_rad) + new_point[1] * math.cos(angle_rad))
#     new_point = [new_point[0] + center_point[0], new_point[1] + center_point[1]] # Reverse the shifting we have done
#     return new_point

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

# def area_of_annulus(inner, outer):
# 	return math.pi * ((inner**2) - (outer**2))

# def make_circle(radius, points, noise):
# 	# returns a list of points defining a circle.
# 	circle = []
# 	angleStep = 2*math.pi / points
# 	for i in range(0,points):
# 		circle.append([(radius * math.cos(i * angleStep) + random.randint(0,noise)),(radius * math.sin(i * angleStep) ) + random.randint(0,noise) ])
# 	return circle

# def mass_of_a_sphere(density, radius):
# 	return density * 4/3 * math.pi * radius**3

# def mass_of_a_circle(density, radius):
# 	return density * 2 * math.pi * radius**2

def semiMinorAxis(a, e):
	return a * math.sqrt(1 - (math.pow(e,2))) # https://math.stackexchange.com/questions/1259945/calculating-semi-minor-axis-of-an-ellipse

def averageOfList(TheList):
    avg = sum(TheList) / len(TheList)
    return avg

def centroidOfPolygon(polygon):
	xValues = []
	yValues = []
	for point in polygon:
		xValues.append(point[0])
		yValues.append(point[1])

	return [averageOfList(xValues), averageOfList(yValues)]

def point_inside_rectangle(point, rect):
	xMin = 0
	xMax = 0
	yMin = 0
	yMax = 0
	for rPoint in rect:
		if rPoint[0] > xMax or xMax == 0: xMax = rPoint[0]
		if rPoint[1] > yMax or yMax == 0: yMax = rPoint[1]
		if rPoint[0] < xMin or xMin == 0: xMin = rPoint[0]
		if rPoint[1] < yMin or yMin == 0: yMin = rPoint[1]

	if point[0] > xMin and point[0] < xMax and point[1] > yMin and point[1] < yMax:
		return True
	else:
		return False

def transformPolygonForLines(polygon):
		n = 0
		points = []

		if len(polygon) == 0:
			return

		firstPoint = polygon[0]
		lastPoint = firstPoint
		points.extend([ int(firstPoint[0]) , int(firstPoint[1])])
		n +=1 
		
		for index, point in enumerate(polygon):
			points.extend([ int(point[0]) , int(point[1])])
			points.extend([ int(point[0]) , int(point[1])])
			lastPoint = point
			n +=2

		# make a segment joining the first one to the last one and then do a repeat to close it off.
		points.extend([ int(firstPoint[0]) , int(firstPoint[1])])
		points.extend([ int(lastPoint[0]) , int(lastPoint[1])])
		points.extend([ int(lastPoint[0]) , int(lastPoint[1])])
		n +=3
		return [n,points]

def transformForView( position ,viewpointObjectPosition, zoom, resolution):

	transformedPosition = position - viewpointObjectPosition 					# offset everything by the position of the viewpointObject, so the viewpoint is at 0,0
	transformedPosition = transformedPosition * zoom  							# shrink or expand everything around the 0,0 point
	transformedPosition[0] = int(transformedPosition[0] + resolution_half[0]) 	# add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
	transformedPosition[1] = int(-transformedPosition[1] + resolution_half[1]) 	# add half the height. and invert it so that it's the right way up in opengl.
	# t: transformed position
	# p: pos
	# v: viewpoint pos
	# z: zoom
	# o: offset
	# t = ((p - v) * z) + o
	return transformedPosition

def antiTransformForView( position ,viewpointObjectPosition, zoom, resolution): # the inverse of transformForView
	# ((t-o)/z)+v = p
	transformedPosition = [(position[0] - resolution_half[0]) / zoom, (-1 * position[1] + resolution_half[1]) / zoom]
	transformedPosition[0] += viewpointObjectPosition[0]
	transformedPosition[1] += viewpointObjectPosition[1]
	return transformedPosition

def isPointIlluminated(point, color, illuminators,viewpointObjectPosition, zoom, resolution):
	thePointIsIlluminated = False
	# n_lum = 0

	# for each illuminator, see if it is close enough to the point to light it up.
	for illuminator in illuminators:
		illuminator.isItLightingThisPoint = False
		illuminator.distanceVector = point - illuminator.transformedPosition #[point[0] - illuminator.transformedPosition[0],point[1] - illuminator.transformedPosition[1]]
		illuminator.distanceScalar = (illuminator.distanceVector[0] ** 2 + illuminator.distanceVector[1] ** 2) ** 0.5 # just like mag(), but doesn't need an numpy vector.
		if illuminator.distanceScalar < illuminator.radius :
			illuminator.isItLightingThisPoint = True
			thePointIsIlluminated = True
			# n_lum +=1 

	# exit out of the function quickly if it is not.
	if not thePointIsIlluminated:
		return color

	# if it is lighting, figure out how bright the light is.
	resultColor = copy.copy(color)
	for illuminator in illuminators:
		if illuminator.isItLightingThisPoint:
			scaledDistance = illuminator.distanceScalar/zoom
			if scaledDistance == 0:
				amountOfLight = illuminator.brightness
			else:
				amountOfLight = illuminator.brightness/scaledDistance
			
			# mix the light color into the point's original color
			resultColor[0] = int(resultColor[0] + amountOfLight * illuminator.color[0])
			if resultColor[0] > 255: resultColor[0] = 255

			resultColor[1] = int(resultColor[1] + amountOfLight * illuminator.color[1])
			if resultColor[1] > 255: resultColor[1] = 255

			resultColor[2] = int(resultColor[2] + amountOfLight * illuminator.color[2])
			if resultColor[2] > 255: resultColor[2] = 255
	return resultColor

def transformPolygonForLinesWithIlluminators(polygon, color, illuminators, viewpointObjectPosition, zoom, resolution):
		# this function prepares a polygon to be rendered by openGL. It accomodates lighting information from the game. It is used to render outlines.
		n = 0
		points = []
		colorstream = []
		firstPoint = polygon[0]
		lastPoint = firstPoint
		points.extend([ int(firstPoint[0]) , int(firstPoint[1])])
		colorstream.extend(isPointIlluminated(firstPoint, color, illuminators, viewpointObjectPosition, zoom, resolution))
		n +=1 
		
		for index, point in enumerate(polygon):
			points.extend([ int(point[0]) , int(point[1])])
			points.extend([ int(point[0]) , int(point[1])])
			colorstream.extend(isPointIlluminated(point, color, illuminators, viewpointObjectPosition, zoom, resolution))
			colorstream.extend(isPointIlluminated(point, color, illuminators, viewpointObjectPosition, zoom, resolution))
			lastPoint = point
			n +=2

		# make a segment joining the first one to the last one and then do a repeat to close it off.
		points.extend([ int(firstPoint[0]) , int(firstPoint[1])])
		colorstream.extend(isPointIlluminated(firstPoint, color, illuminators, viewpointObjectPosition, zoom, resolution))
		points.extend([ int(lastPoint[0]) , int(lastPoint[1])])
		colorstream.extend(isPointIlluminated(lastPoint, color, illuminators, viewpointObjectPosition, zoom, resolution))
		points.extend([ int(lastPoint[0]) , int(lastPoint[1])])
		colorstream.extend(isPointIlluminated(lastPoint, color, illuminators, viewpointObjectPosition, zoom, resolution))
		n +=3

		return [n,points, colorstream]

def transformPolygonForTriangleFan(polygon):
	# this function transforms a polygon into a fan-shaped array of triangles that can be rendered by openGL. It is used for rendering filled-in areas of convex objects.
	# the number of points returned by this function is always 1.5n + 5, where n is the number of points in polygon.
	centroid = centroidOfPolygon(polygon)
	centroid[0] = int(centroid[0])
	centroid[1] = int(centroid[1])
	transformedPoints = []
	length = len(polygon)
	n = 0
	# repeat start
	transformedPoints.append(polygon[0][0])
	transformedPoints.append(polygon[0][1])
	n+= 1

	# take every first and second point, and make a triangle between them and the centroid.
	closeTriangle = False
	for point in polygon:
		n+= 1
		transformedPoints.append(point[0])
		transformedPoints.append(point[1])
		if closeTriangle:
			transformedPoints.append(centroid[0])
			transformedPoints.append(centroid[1])
			n+= 1
			closeTriangle = False

		else:
			closeTriangle = True

	# you may need to run this again to close the last triangle if the polygon had an uneven number of points.
	if closeTriangle:
		transformedPoints.append(centroid[0])
		transformedPoints.append(centroid[1])
		n+= 1

	# finally, you need to create another triangle to fill the space between the first and last vertices.
	# this snippet adds that triangle plus the repeat end sequence required by pyglet. 
	transformedPoints.append(polygon[length-1][0])
	transformedPoints.append(polygon[length-1][1])
	transformedPoints.append(polygon[0][0])
	transformedPoints.append(polygon[0][1])
	transformedPoints.append(centroid[0])
	transformedPoints.append(centroid[1])
	transformedPoints.append(centroid[0])
	transformedPoints.append(centroid[1])
	n+= 4
	return [n,transformedPoints]

def unwrapAtmosphereForGradient(n, inside_verts, outside_verts, inside_color, outside_color): 
		# this function prepares an atmosphere from the game, which is an annulus composed of two concentric circles of points, to be rendered by openGL. It allows a a color gradient from the inner edge to the outer edge.
		bitstream = []
		colorstream = []
		nn = 0

		#repeat start
		bitstream.append(inside_verts[0][0])
		bitstream.append(inside_verts[0][1])
		colorstream.append(inside_color[0])
		colorstream.append(inside_color[1])
		colorstream.append(inside_color[2])
		colorstream.append(inside_color[3])
		nn += 1

		# now go around and for each pair of inner and outer points create two triangles and a color gradient.
		# the bitstream order is: from https://stackoverflow.com/questions/20394727/gl-triangle-strip-vs-gl-triangle-fan
		# inner n, outer n, inner n + 1, outer n + 1
		# the colorstream order is (index matching bitstream order above):
		# full color, no color, full color, no color
		for index in range(0,n):
			bitstream.append(inside_verts[index][0])
			bitstream.append(inside_verts[index][1])
			bitstream.append(outside_verts[index][0])
			bitstream.append(outside_verts[index][1])
			colorstream.append(inside_color[0])
			colorstream.append(inside_color[1])
			colorstream.append(inside_color[2])
			colorstream.append(inside_color[3])
			colorstream.append(outside_color[0])
			colorstream.append(outside_color[1])
			colorstream.append(outside_color[2])
			colorstream.append(outside_color[3])
			nn += 2

		# closing the annulus is done by creating a pair of triangles between the first and last segments of the inner and outer edges.
		# it may fail in situations with an uneven number of points, so far we have always used even numbers of points.
		bitstream.append(inside_verts[0][0])
		bitstream.append(inside_verts[0][1])
		colorstream.append(inside_color[0])
		colorstream.append(inside_color[1])
		colorstream.append(inside_color[2])
		colorstream.append(inside_color[3])
		nn += 1

		bitstream.append(outside_verts[0][0])
		bitstream.append(outside_verts[0][1])
		colorstream.append(outside_color[0])
		colorstream.append(outside_color[1])
		colorstream.append(outside_color[2])
		colorstream.append(outside_color[3])
		nn += 1

		bitstream.append(outside_verts[0][0])
		bitstream.append(outside_verts[0][1])
		colorstream.append(outside_color[0])
		colorstream.append(outside_color[1])
		colorstream.append(outside_color[2])
		colorstream.append(outside_color[3])
		nn += 1

		return [nn, bitstream, colorstream]

def renderAConvexPolygon(batch, polygon, viewpointObjectPosition, zoom, resolution, color, outlineColor=None, illuminators=None):
	# this function is used to route game objects to the appropriate rendering method, depending on their properties.
	if len(polygon) > 0:
		transformedPoints = transformPolygonForTriangleFan(polygon)
		batch.add(transformedPoints[0], pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',transformedPoints[1]), ('c4B',color*transformedPoints[0]))
		if outlineColor is not None:
			if illuminators is not None:
				transformedPoints = transformPolygonForLinesWithIlluminators(polygon, outlineColor, illuminators, viewpointObjectPosition, zoom, resolution)
				batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',transformedPoints[2]))
			else:
				transformedPoints = transformPolygonForLines(polygon)
				batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',outlineColor*transformedPoints[0]))


class SavedGame():
	def __init__(actors, currentSystem):
		self.actors = actors
		self.currentSystem = currentSystem



class World():
	def __init__(self):
		self.time = 0 							# the number of timesteps that have passed in-universe. used for physics and orbital calculations.
		self.space = pymunk.Space()				
		self.space.gravity = (0.0, 0.0)			# pymunk's own gravity is turned off, so we can use our own.
		# self.gravitationalConstant = gravitationalConstant		# sets the strength of gravity
		self.dragCoefficient = 0.005 			# Sets the strength of atmospheric drag
		self.actors = []
		self.attractors = []
		resolution = resolution_half
		self.ch = self.space.add_collision_handler(0, 0)
		self.ch.post_solve = self.handle_collision
		self.viewpointObject = None
		self.viewpointObjectIndex = 0
		self.player = None
		self.zoom = 1 # the actual applied zoom number.
		self.pan = [0,0]
		self.rotate = 0
		self.timestepSize = 0.2/60.0 #1.0/60.0
		self.showHUD = False
		self.hudbatch  = pyglet.graphics.Batch()
		# pyglet.gl.glLineWidth(2)# None # used to hold a precooked copy of the HUD so it doesn't need recomputed all the time.
		self.paused = True

		self.buildMenu = False				# used to toggle between the game and the build screen
		self.availableModules = []			# a list of potentially useable modules that the player has in 'inventory'
		self.modulesInUse = []  # a list of modules that the player has dragged onto the screen to make a ship
		self.buildDraggingModule = None
		self.mouseCursorPosition = []

		self.navCircleLines = []
		self.n_navcircle_lines = 32
		self.navcircleLinesLength = 10
		self.navcircleInnerRadius = 250

		self.projectiles = []
		self.illuminators = []

		self.topLimit = antiTransformForView( resolution  ,[0,0], self.zoom, resolution)
		self.bottomLimit = antiTransformForView( [0,0] ,[0,0], self.zoom, resolution)

		self.mapView = False
		self.galaxy = {} # a list of all the solar systems the player can visit.
		self.currentSystem = None
		self.hyperjumpDestinationIndex = 0

		self.playerTargetIndex = None # index of which actor in the list the player is targeting.
		self.targetfutureTime = 1
		self.scrollLockTargetFutureTime = None

	def gravityForce(self, actorPosition, attractorPosition, attractorMass):
		# this function returns the force of gravity between an actor and an attractor.
		distance = attractorPosition - actorPosition
		magnitude = mag(distance)
		gravity = gravitationalConstant * attractorMass
		appliedGravity = gravity/(magnitude**2)
		components = numpy.divide(distance, magnitude)
		force = components * appliedGravity * self.timestepSize * 4.75
		return force

	def gravitate(self, actor, force):
		# this function is used to apply a force to an actor. It is principally used to apply gravity to the actors when they are not freefalling.
		rotatedForce = Vec2d(force[0], force[1])
		rotatedForce = rotatedForce.rotated(-actor.body.angle)
		actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])

	def getModuleFromCursorPosition(self, cursorPositionRaw):
		cursorPosition = cursorPositionRaw
		smellyCursor = [0,0]
		smellyCursor[0] = cursorPositionRaw[0]
		smellyCursor[1] = -cursorPositionRaw[1] +  resolution[1]

		# when picking a module out of the list, you don't need to transform into game coordinates.
		for listItem in self.availableModuleListItems:
			if pointInRect(smellyCursor, listItem.boundingRectangle):
				self.availableModuleListItems.remove(listItem)
				return listItem.module

		# put bounding boxes around the modules in-game coordinates, then antitransform the mouse and see if it landed on any of them.
		for module in self.modulesInUse:
			transformedPoints = []
			for point in module.points:
				transformedPoint = [0,0]
				transformedPoint[0] = (point[0] + module.offset[0])
				transformedPoint[1] = (point[1] + module.offset[1])
				transformedPoints.append(transformedPoint)
			boundingBox = boundPolygon(transformedPoints)
			if pointInRect( self.antiTransformForBuild(cursorPositionRaw) , boundingBox):
				self.modulesInUse.remove(module)
				return module

	def add(self, thing):  
		self.space.add(thing.body, thing.shape)
		if thing.__class__ == Actor:
			self.actors.append(thing)
		elif thing.__class__ == Attractor:
			self.attractors.append(thing)

		self.player = self._getPlayer()
		# self.viewpointObject = self.player

	def destroyActor(self, actor):
		# 
		self.space.remove(actor.shape, actor.body)
		# if actor in self.actors:
		# 	print('removing actor: ' + str(actor.name))
		# 	self.actors.remove(actor)

	def destroyAttractor(self,attractor):
		# self.space.remove(attractor.shape, attractor.body)
		self.space.remove(attractor.shape, attractor.body)
		# self.space.remove(attractor.body)
		# if attractor in self.attractors:
		# 	print('removing attractor: ' + str(attractor.planetName))
		# 	self.attractors.remove(attractor)

	def decomposeActor(self, actor, modules):
		listLength = len(actor.modules)
		if listLength == 1:
			return # the actor is already fully decomposed, destroy it if you want
		else:
			# create a new actor, minus the module
			addAPlayerFragment = False
			if actor.isPlayer:
				addAPlayerFragment = True
			for index, module in enumerate(actor.modules):

				# create the module on it's own as a new actor
				fragmentPosition = [actor.body.position[0] + (module.offset[0] * math.cos(actor.body.angle)), actor.body.position[1] +  (module.offset[1] * math.sin(actor.body.angle))]
				module.offset = [0,0]

				if addAPlayerFragment:
					self.add(Actor(actor.name + ' fragment', [module], fragmentPosition, actor.body.velocity, True))
					self.player = self._getPlayer()
					addAPlayerFragment = False
				else:
					self.add(Actor(actor.name + ' fragment', [module], fragmentPosition, actor.body.velocity, False))

			# remove the actor from the space.
			self.destroyActor(actor)

	def physics(self):
		for actor in self.actors:
			actor.exemptFromGravity = False

		# iterate the physics engine. This is done first, so that changes can be reacted to (i.e. to exempt things from gravity if they collide with a planet, to prevent sinking).
		self.space.step(self.timestepSize)
		self.time += self.timestepSize

		for actor in self.actors:
			if self.timestepSize  * 3 * 100 > 50:
				actor.body.angle = actor.setPoint
				actor.body._set_angular_velocity(0)


			if math.isnan(mag(actor.body.position)):
				self.destroyActor(actor)
				continue

			# explode all the projectiles
			destroyed = False
			for module in actor.modules:
				if module.moduleType == "cannonshell 10":
					module.lifetime -= 1
					if module.lifetime < 0:
						self.destroyActor(actor)
						destroyed = True
						break
			if destroyed:
				continue

			actor.doResources(self.timestepSize)

			# figure out which attractor you are orbiting
			strongestForce = None
			strongestAttractor = None
			for attractor in self.attractors:
				force = self.gravityForce(actor.body.position, attractor.body.position, attractor.body.mass)
				if strongestAttractor is None or mag(force) > mag(strongestForce):
					strongestForce = force
					strongestAttractor = attractor

			# do atmospheric drag
			for attractor in self.attractors:
				if attractor.atmosphere is not None:
					actorHeightFromAttractorCenter = mag(actor.body.position - attractor.body.position)
					if actorHeightFromAttractorCenter < (attractor.radius + attractor.atmosphere.height):
						
						actor.leaveFreefall(0)

						naturalDepth = 1 - ((actorHeightFromAttractorCenter - attractor.radius) / attractor.atmosphere.height) # this is a number between 0 and 1 which is signifies the actors depth into this atmosphere layer.

						density = (attractor.atmosphere.topDensity + (naturalDepth * attractor.atmosphere.bottomDensity))**1.5 # not quite squared, but still exponential

						dragForceX = self.timestepSize * self.dragCoefficient * ((density * actor.body.velocity[0]**2) /2) * actor.body.mass # using mass as a placeholder because i don't have a drag frontal area calculation yet. but it still needs to apply to bigger things more.
						dragForceY = self.timestepSize * self.dragCoefficient * ((density * actor.body.velocity[1]**2) /2) * actor.body.mass # using mass as a placeholder because i don't have a drag frontal area calculation yet. but it still needs to apply to bigger things more.
						
						if actor.body.velocity[0] > 0:
							dragForceX = -abs(dragForceX)
						if actor.body.velocity[0] < 0:
							dragForceX = abs(dragForceX)

						if actor.body.velocity[1] > 0:
							dragForceY = -abs(dragForceY)
						if actor.body.velocity[1] < 0:
							dragForceY = abs(dragForceY)

						rotatedForce = Vec2d(dragForceX, dragForceY)
						rotatedForce = rotatedForce.rotated(-actor.body.angle)
						actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])
						
			# when you enter a new sphere of influence, regenerate the orbit information
			if strongestAttractor is not actor.orbiting or actor.orbiting is None:
				actor.leaveFreefall(0)
				actor.orbiting = strongestAttractor

			# shoot the guns
			for module in actor.modules:
				if module.moduleType == 'cannon 10':
					if actor.keyStates['Fire']:
						module.active= True
						actor.keyStates['Fire'] = False
					else:
						module.active = False
					
					if module.active:
						Nirn.shootABullet(module, Nirn.player)

			# figure out the angle to steer the ship.
			if actor.isPlayer:
				if actor.keyStates['left']:
					actor.setPoint -= 0.03
				if actor.keyStates['right']:
					actor.setPoint += 0.03

			angleDifference = actor.setPoint - actor.body.angle

			# all rotating actors experience a slight drag which slows their rotation. (it's more fun that way).
			if actor.body.angular_velocity != 0:
				correctionForce = actor.body.angular_velocity * 100 * actor.body.mass * self.timestepSize
				actor.body.apply_impulse_at_local_point([-correctionForce,0], [0,-10]) # 10 is the moment arm, in reality it should be equal to the actor's radius.
				actor.body.apply_impulse_at_local_point([correctionForce,0], [0,10])
				if abs(actor.body.angular_velocity) < 1/100000:
					actor.body.angular_velocity = 0

			# draw the module effects like gun flashes and engine flames
			for module in actor.modules:
				if module.enabled:		

					if module.moduleType == 'RCS':
						module.active = actor.keyStates['left'] or actor.keyStates['right']
					elif  module.moduleType == 'engine 10' or module.moduleType == 'engine 100':
						module.active = actor.keyStates['up']


					if module.active:
						if module.effect is not None and module.effect.illuminator is not None:
							position = actor.body.position + module.offset + module.effect.illuminator.offset
							position = rotate_point(position, actor.body.angle, actor.body.position)
							module.effect.illuminator.position = position
							self.illuminators.append(module.effect.illuminator)

				# turn the gun off until it has done a cooldown.
				if module.moduleType == 'cannon 10':
					module.active = False
	
			# figure out if the actor is freefalling by seeing if any engines or collisions have moved it.
			if actor.doModuleEffects(self.timestepSize):
				actor.leaveFreefall(0)

			# if it is freefalling, move it along the orbital track.
			if actor.orbiting is not None:
				actor.nadir = math.atan2( actor.orbiting.body.position[1] - actor.body.position[1], actor.orbiting.body.position[0] - actor.body.position[0] )
				actor.zenith = actor.nadir + math.pi
			
			if not actor.exemptFromGravity:
				if not actor.freefalling:
					# if it is not freefalling, add some gravity to it so that it moves naturally, and try to recalculate the orbit.
					if not actor.exemptFromGravity:
						if strongestForce is not None:
							self.gravitate(actor, strongestForce)

					if actor.stepsToFreefall > 0:
						actor.stepsToFreefall -= 1
					else:
						actor.freefalling = True
						actor.exemptFromGravity = False
				
					try:
						actor.orbit = Orbit.fromStateVector(numpy.array([actor.body.position[0] - actor.orbiting.body.position[0] ,actor.body.position[1] - actor.orbiting.body.position[1],1]), numpy.array([actor.body.velocity[0] ,actor.body.velocity[1] ,1]), actor.orbiting.APBody, Time('2000-01-01 00:00:00'), actor.name + " orbit around " + actor.orbiting.planetName)
						# generate the orbit points once only.
						if actor.orbit is not None:
							actor.orbitPoints = []
							n_points = 64
							sliceSize =  (2 * math.pi / n_points)
							onMirrorHalf = False
							for i in range(0,n_points):
								eccentricAnomaly =i *sliceSize
								trueAnomaly = 2 * math.atan2(  math.sqrt((1 + actor.orbit.e))* math.sin(eccentricAnomaly/2) , math.sqrt((1 - actor.orbit.e)) * math .cos(eccentricAnomaly/2)  )
								temp_vec3d = actor.orbit.cartesianCoordinates(trueAnomaly)	
								actor.orbitPoints.append((temp_vec3d[0] + actor.orbiting.body.position[0], temp_vec3d[1] + actor.orbiting.body.position[1]))
					except:
						actor.orbit = None

					if actor.orbit is not None:
						futureSteptAn = actor.orbit.tAnAtTime(self.timestepSize)
						futureStepCoordinates = actor.orbit.cartesianCoordinates(futureSteptAn)
						adjustedFutureStep = [futureStepCoordinates[0] + actor.orbiting.body.position[0] , futureStepCoordinates[1] + actor.orbiting.body.position[1]]
						actor.prograde = math.atan2( adjustedFutureStep[1] - actor.body.position[1], adjustedFutureStep[0] - actor.body.position[0] )
						actor.retrograde = actor.prograde + math.pi

				else:
					if actor.orbit is not None:
						actor.orbit.updTime(self.timestepSize)
						cartesian = actor.orbit.cartesianCoordinates(actor.orbit.tAn)
						actor.body.position =  [cartesian[0] + actor.orbiting.body.position[0] ,cartesian[1] + actor.orbiting.body.position[1]]

						# you also must update the actor's velocity, or else when it leaves the track it will have the same velocity it entered with, leading to weird jumps.
						trackSpeed = actor.orbit.getSpeed(actor.orbit.tAn)

						# there is almost definitely a way to figure this out from the ellipse's properties. You would need to find tangent to the ellipse. But I figured it out by computing the position one step into the future, and then finding the angle between positions.
						futureSteptAn = actor.orbit.tAnAtTime(self.timestepSize)
						futureStepCoordinates = actor.orbit.cartesianCoordinates(futureSteptAn)
						adjustedFutureStep = [futureStepCoordinates[0] + actor.orbiting.body.position[0] , futureStepCoordinates[1] + actor.orbiting.body.position[1]]
						actor.prograde = math.atan2( adjustedFutureStep[1] - actor.body.position[1], adjustedFutureStep[0] - actor.body.position[0] )
						actor.retrograde = actor.prograde + math.pi
						actor.body.velocity = [trackSpeed * math.cos(actor.prograde), trackSpeed * math.sin(actor.prograde)]

			# let the ai drive the ship. this comes after orbit calculation because it needs valid orbits
			actor.flightComputer()

			# perform hyperspace travel.
			if actor.jumping:
				actor.jumping = False
				self.hyperspaceJump(actor)

	def rotatePolygon(self, points, angle):
		return Rotate2D(points,(0,0),angle)

	def transformForBuild(self, position):
		# map a position in the game world, where 0,0 is in a corner and numbers are very large, onto pixels on the screen with 0,0 in the middle. Handles zooming and offsetting. This one is for what you see in the build menu.
		transformedPosition = [0,0] #* self.zoom  # shrink or expand everything around the 0,0 point
		transformedPosition[0] = -position[0] * self.zoom
		transformedPosition[1] = -position[1] * self.zoom
		transformedPosition[0] += 0.5 * resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
		transformedPosition[1] += 0.5 * resolution[1] # add half the height.
		return transformedPosition

	def antiTransformForBuild(self, position):
		# performs the inverse operation to transformForBuild, used to map the mouse cursor to coordinates in the game world.
		transformedPosition = [0,0] #* self.zoom  # shrink or expand everything around the 0,0 point
		transformedPosition[0] = (-position[0]) + 0.5 * resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
		transformedPosition[1] = (-position[1]) + 0.5 * resolution[1] # add half the height.
		transformedPosition[0] = transformedPosition[0] / self.zoom
		transformedPosition[1] = transformedPosition[1] / self.zoom
		return transformedPosition

	def drawModuleForList(self, main_batch, module, iconSize, position):
		transformedPoints = []
		for point in module.points:
			transformedPoint = [0,0]
			transformedPoint[0] = ((point[0] * iconSize) + position[0])
			transformedPoint[1] = -((point[1] * iconSize) + position[1]) + resolution[1]
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, resolution, module.color, module.outlineColor)

	def drawModuleForDrag(self,main_batch, module, position):
		rotatedPoints = rotate_polygon(module.points, module.angle)
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			rotatedPoint = self.transformForBuild(rotatedPoint)
			rotatedPoint[0] += position[0] - resolution_half[0]
			rotatedPoint[1] += position[1] - resolution_half[1]
			transformedPoint = rotatedPoint # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, resolution, module.color, module.outlineColor)

	def drawModuleForBuild(self, main_batch, module):
		rotatedPoints = rotate_polygon(module.points, module.angle)
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			rotatedPoint[0] += module.offset[0]
			rotatedPoint[1] += module.offset[1]
			transformedPoint = self.transformForBuild(rotatedPoint) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, resolution, module.color, module.outlineColor)

	def drawModuleEffects(self, main_batch, module, actor):
		rotatedPoints = module.effect.points
		rotatedPoints = rotate_polygon(rotatedPoints,module.angle)  # orient the polygon according to the body's current direction in space.
		rotatedPoints = rotate_polygon(rotatedPoints, actor.body.angle, [-(module.offset[0] + module.effect.offset[0]), -(module.offset[1] + module.effect.offset[1])])
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			try:
				transformedPoint = rotatedPoint + actor.body.position + module.offset + module.effect.offset
				transformedPoint = transformForView(transformedPoint,self.viewpointObject.body.position, self.zoom, resolution ) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
				transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
			except:
				pass
			
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, resolution,  module.effect.color, module.effect.outlineColor)#, self.illuminators)

	def drawColorIndicator(self, color,  position, size, main_batch):
		# draw a simple colored square that can be used for visual indication.

		points = [[-size,size],[size,size],[size,-size],[-size,-size]]
		newpos = transformForView(position, self.viewpointObject.body.position, self.zoom, resolution )
		nupoints = []
		for point in points:
			nupoints.append([int(point[0] + newpos[0]), int(point[1] + newpos[1])])
			
		renderAConvexPolygon(main_batch, nupoints, self.viewpointObject.body.position, self.zoom, resolution,  color)

	def drawModule(self, actor, module, main_batch): # draw the outline of the module.
		rotatedPoints = module.points
		rotatedPoints = rotate_polygon(rotatedPoints,module.angle)  # orient the polygon according to the body's current direction in space.
		rotatedPoints = rotate_polygon(rotatedPoints, actor.body.angle, [-module.offset[0], -module.offset[1]])
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			try:
				transformedPoint = transformForView(rotatedPoint + actor.body.position + module.offset,self.viewpointObject.body.position, self.zoom, resolution ) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
				transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
			except:
				pass
			
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, resolution,  module.color, module.outlineColor, self.illuminators)

		if module.effect is not None:
			if module.enabled and module.active:
				self.drawModuleEffects(main_batch, module, actor)
	
	def drawScreenFill(self, main_batch):
		# this function is used to fill the build menu with a white backdrop.
		color= [0,0,0,255]
		if self.buildMenu:
			color = [220,220,220,255]

		if self.showHUD:
			fillTriangles = [xCompensation,0, xCompensation,0, xCompensation,resolution[1], resolution[0]-xCompensation,resolution[1], resolution[0]-xCompensation,0, xCompensation,0, xCompensation,0 ]
			main_batch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',color*7))

		else:
			fillTriangles = [0,0, 0,0, 0,resolution[1], resolution[0],resolution[1], resolution[0],0, 0,0, 0,0 ]
			main_batch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',color*7))

	def drawHUDmask(self, main_batch):
		hudBackgroundColor = (25,25,25,255)

		fillTriangles = [100,0, 100,0, 100,resolution[1], 200,resolution[1], 200,0, 100,0, 100,0 ]
		main_batch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',hudBackgroundColor*7))

		xRightLimit = resolution[0] - 100
		fillTriangles = [resolution[0],0, resolution[0],0, resolution[0],resolution[1], xRightLimit,resolution[1], xRightLimit,0,  resolution[0],0,  resolution[0],0 ]
		main_batch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',hudBackgroundColor*7))


	def hudbatch_render(self):
		hudBackgroundColor = (25,25,25,255)

		# self.hudbatch.draw()
		self.hudbatch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)
		# self.hudbatch  = pyglet.graphics.Batch()
		fillTriangles = [0,0, 0,0, 0,resolution[1], 200,resolution[1], 200,0, 0,0, 0,0 ]
		self.hudbatch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',hudBackgroundColor*7))

		xRightLimit = resolution[0] - 200
		fillTriangles = [resolution[0],0, resolution[0],0, resolution[0],resolution[1], xRightLimit,resolution[1], xRightLimit,0,  resolution[0],0,  resolution[0],0 ]
		self.hudbatch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',hudBackgroundColor*7))

		self.drawHUDpermanentComponents(self.hudbatch)


	def drawActor(self, actor, main_batch):
		if actor.__class__ is Actor:
			for module in actor.modules:
					self.drawModule(actor, module, main_batch)

		if actor.__class__ is Atmosphere:

			transformedInnerPoints = []
			transformedOuterPoints = []

			for index, point in enumerate(actor.innerPoints): 
				transformedInnerPoints.append(transformForView(point,self.viewpointObject.body.position, self.zoom, resolution))
			for index, point in enumerate(actor.outerPoints): 
				transformedOuterPoints.append(transformForView(point,self.viewpointObject.body.position, self.zoom, resolution))

			rendering = unwrapAtmosphereForGradient(actor.n_points, transformedInnerPoints, transformedOuterPoints, actor.color, actor.outerColor)
			main_batch.add(rendering[0], pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',rendering[1]), ('c4B',rendering[2]))
		
		if actor.__class__ is Attractor or actor.__class__ is Cloud:
			transformedPoints = []
		
			# work backwards to figure out what coordinates in the game world correspond to being offscreen in the view.
			self.topLimit = antiTransformForView( resolution  ,self.viewpointObject.body.position, self.zoom, resolution)
			self.bottomLimit = antiTransformForView( [0,0] ,self.viewpointObject.body.position, self.zoom, resolution)

			prevPointInside = True
			pointInside = True
			firstPoint = actor.points[0]+ actor.body.position 
			prevPoint =	 firstPoint
			currentPoint = prevPoint
			
			# figure out the index at the end of the array.
			onTheVeryLastPoint = len(actor.points) - 1

			for index, point in enumerate(actor.points):
				# update the information for the current and previous points.
				prevPoint = currentPoint
				transformedPoint = point + actor.body.position
				currentPoint = transformedPoint
				prevPointInside = pointInside
				pointInside = False

				# figure out if the point is inside the view field.
				if transformedPoint[0] > self.bottomLimit[0] and transformedPoint[0] < self.topLimit[0]:
					pointInside = True

				if  transformedPoint[1] > self.bottomLimit[1] and transformedPoint[1] < self.topLimit[1]:
					pointInside = True

				# if one of the points is outside the view field, and one is inside, you want to draw it, so that the line leaving the screen looks correct.
				if (pointInside and not prevPointInside) or (prevPointInside and not pointInside):
					transformedPoint = transformForView(prevPoint,self.viewpointObject.body.position, self.zoom, resolution) 
					transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
					transformedPoint = transformForView(currentPoint,self.viewpointObject.body.position, self.zoom, resolution) 
					transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])

					# on the last point, you want to make a connection between the first and last point in the polygon. It should be treated the same as the other lines.
					if index == onTheVeryLastPoint:
						transformedPoint = transformForView(firstPoint,self.viewpointObject.body.position, self.zoom, resolution) 
						transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])

				# if both points are inside, you want to draw both. But you've already drawn lastPoint, so you only need to add the current point.
				elif pointInside and prevPointInside:
					transformedPoint = transformForView(currentPoint,self.viewpointObject.body.position, self.zoom, resolution) 
					transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])

					if index == onTheVeryLastPoint:
						transformedPoint = transformForView(firstPoint,self.viewpointObject.body.position, self.zoom, resolution) 
						transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])

			# if the polygon was completely offscreen, the length of the array will be 0, and you can skip it.
			if len(transformedPoints) > 0:
				renderAConvexPolygon(main_batch, transformedPoints,self.viewpointObject.body.position, self.zoom, resolution,  actor.color, actor.outlineColor)	

	def getActorFromBody(self, body):
		for actor in self.actors:
			if body == actor.body:
				return actor
		return None
	def getAttractorFromBody(self, body):
		for attractor in self.attractors:
			if body == attractor.body:
				return attractor
		return None

	def handle_collision(self, arbiter, space, data):
		shapes = arbiter._get_shapes()
		ke = arbiter._get_total_ke()

		bodyA = shapes[0]._get_body()
		actorA = self.getActorFromBody(bodyA)

		bodyB = shapes[1]._get_body()
		actorB = self.getActorFromBody(bodyB)

		if actorA is not None and actorB is not None:
			if (actorA.isPlayer and len(actorB.modules) == 1) and actorB.modules[0].moduleType != 'nano orbital':
				self.availableModules.append(actorB.modules[0])
				self.destroyActor(actorB)
				return

			if (actorB.isPlayer and len(actorA.modules) == 1) and actorA.modules[0].moduleType != 'nano orbital':
				self.availableModules.append(actorA.modules[0])
				self.destroyActor(actorA)
				return

		#handle collisions between a planet and a ship differently than collisions between two ships
		attractorCollision = False
		if actorA is None:
			attractorCollision = True
		else:
			# figure out the energy applied to each vessel, so you can tell how destructive it is.
			ke_A = ke / actorA.mass
			actorA.leaveFreefall(1)

			# if the energy is more than the actor can take, explode the actor.
			if ke_A > actorA.decompEnergy:
				self.decomposeActor(actorA, actorA.modules)

		if actorB is None:
			attractorCollision = True
		else:
			ke_B = ke / actorB.mass
			actorB.leaveFreefall(1)
			if ke_B > actorB.decompEnergy:
				self.decomposeActor(actorB, actorB.modules)

		if attractorCollision:
			if actorA is not None:
				actorA.exemptFromGravity = True

				#if the velocity is very small, kill it completely. Along with gravity exemption, this absolutely solves jitter and sinking.
				if mag(bodyA.velocity - bodyB.velocity) < 50:
					bodyA.velocity = bodyB.velocity

			if actorB is not None:
				actorB.exemptFromGravity = True

				if mag(bodyA.velocity - bodyB.velocity) < 50:
					bodyB.velocity = bodyA.velocity

	def _getPlayer(self):
		for actor in self.actors:
			if actor.isPlayer:
				return actor

	def drawAPOrbit(self, main_batch, actor, orbit, attractor, color):
		# This function visually renders an AstroPynamics orbit object.
		if actor.orbit is None:
			return
		if len(actor.orbitPoints) > 0:
			points = []
			for point in actor.orbitPoints:
				if not math.isnan(point[0]) and not math.isnan(point[1]):
					transformedPoint = transformForView(point,self.viewpointObject.body.position, self.zoom, resolution)
					points.append(transformedPoint)

			transformedPoints = transformPolygonForLines(points)
			if transformedPoints is not None:
				main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[100,100,100,255]*(transformedPoints[0])))

	def drawModuleListItem(self, main_batch, listItem, index):
		# draw one of the modules in the list in the build menu.
		buildListSpacing = 30
		itemSize = 2 * mag(numpy.array(getFarthestPointInPolygon(listItem.module.points)))
		iconSize = buildListSpacing / itemSize

		gnarlypoints = []
		for point in listItem.module.points:
			gnarlypoints.append([point[0] * iconSize + buildListSpacing,point[1] * iconSize + (index * buildListSpacing ) ])

		listItem.boundingRectangle = boundPolygon(gnarlypoints)

		self.drawModuleForList(main_batch, listItem.module, iconSize, [buildListSpacing, index * buildListSpacing] )

	def drawHUDListItemQuantity(self,quantity, index, listCorner, main_batch):
		HUDlistItemSpacing = 15
		listXPosition = 10 + 100
		fontSize = 12

		color = (200,200,200,255)

		if quantity is None:
			return index + 1

		if listCorner == 'bottom left':
			label = pyglet.text.Label(str(quantity) ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=listXPosition, y=HUDlistItemSpacing * index,
	                      color=color,
	                      align="left",
	                      batch=main_batch)
		elif listCorner == 'top right':
			label = pyglet.text.Label(str(quantity) ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=resolution[0] - 100, y=resolution[1]-(HUDlistItemSpacing * index),
	                      color=color,
	                      align="right",
	                      batch=main_batch)
		elif listCorner == 'top left':
			label = pyglet.text.Label(str(quantity) ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=listXPosition, y=resolution[1]-(HUDlistItemSpacing * index),
	                      color=color,
	                      align="left",
	                      batch=main_batch)
		elif listCorner == 'bottom right':
			label = pyglet.text.Label(str(quantity) ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=resolution[0] - 100, y=HUDlistItemSpacing * index,
	                      color=color,
	                      align="right",
	                      batch=main_batch)

		# previousHUDstrings[index] = string
		previousHUDlistQuantities[index] = quantity
		# previousHUDlabels[index] = label
		previousHUDQuantityLabels[index] = label
		# label.draw()
		# main_batch.add(label)

		return index + 1

	def drawHUDListItemLabel(self,string, index, listCorner, main_batch):
		HUDlistItemSpacing = 15
		listXPosition = 10 
		fontSize = 12

		color = (100,100,100,255)

		if len(string) == 0:
			return index + 1

		if listCorner == 'bottom left':
			label = pyglet.text.Label(string ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=listXPosition, y=HUDlistItemSpacing * index,
	                      color=color,
	                      align="left",
	                      batch=main_batch)
		elif listCorner == 'top right':
			label = pyglet.text.Label(string ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=resolution[0] - 200 + listXPosition, y=resolution[1]-(HUDlistItemSpacing * index),
	                      color=color,
	                      align="right",
	                      batch=main_batch)
		elif listCorner == 'top left':
			label = pyglet.text.Label(string ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=listXPosition, y=resolution[1]-(HUDlistItemSpacing * index),
	                      color=color,
	                      align="left",
	                      batch=main_batch)
		elif listCorner == 'bottom right':
			label = pyglet.text.Label(string ,
	                      font_name='Times New Roman',
	                      font_size=fontSize,
	                      x=resolution[0] - 200+ listXPosition, y=HUDlistItemSpacing * index,
	                      color=color,
	                      align="right",
	                      batch=main_batch)

		previousHUDstrings[index] = string
		previousHUDlabels[index] = label

		return index + 1

	def createHUDNavCircle(self): # this function predraws the nav circle. because it is just static lines, it does not need to be recalculated every frame.
		for n in range(0,self.n_navcircle_lines):
			angle = n * (2 * math.pi / self.n_navcircle_lines)
			start = ((self.navcircleInnerRadius * math.cos(angle)) + (resolution[0]*0.5) , (self.navcircleInnerRadius* math.sin(angle)) +(resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius + self.navcircleLinesLength) * math.cos(angle)+ (resolution[0]*0.5) , (self.navcircleInnerRadius + self.navcircleLinesLength) * math.sin(angle)+ (resolution[1]*0.5))
			self.navCircleLines.append([start, end])

	def drawHUDpermanentComponents(self,main_batch):
		if self.player is None:
			return
		i = 1
		hudList = {}
		for resource, quantity in list(self.viewpointObject.availableResources.items()):
			hudList[resource] = quantity
		for resource, quantity in list(self.viewpointObject.storagePool.items()):
			if resource in hudList:
				hudList[resource] += self.viewpointObject.storagePool[resource]
			else:
				hudList[resource] = self.viewpointObject.storagePool[resource]

		for availableResource, availableQuantity in list(hudList.items()):
			i = self.drawHUDListItemLabel(str(availableResource) + ': ', i, 'top left', self.hudbatch)


		i = self.drawHUDListItemLabel('', i, 'top left', self.hudbatch) # blank line as a separator
		i = 1

		printFreefalling = False
		if self.viewpointObject.freefalling or self.viewpointObject.stepsToFreefall == 0:
			printFreefalling = True

		i = self.drawHUDListItemLabel('freefalling: ', i, 'bottom left', self.hudbatch)
		i = self.drawHUDListItemLabel('landed: ', i, 'bottom left', self.hudbatch)
		if self.player.orbiting is not None:
			i = self.drawHUDListItemLabel('orbiting: ', i, 'bottom left', self.hudbatch)
		i = self.drawHUDListItemLabel('', i, 'bottom left', self.hudbatch) # blank line as a separator

		i = self.drawHUDListItemLabel('player: ', i, 'bottom left', self.hudbatch)
		i = self.drawHUDListItemLabel('time accel: ', i, 'bottom left', self.hudbatch)
		i = self.drawHUDListItemLabel('zoom: ', i, 'bottom left', self.hudbatch)
		i = self.drawHUDListItemLabel('paused: ', i, 'bottom left', self.hudbatch)

		i = 1

		# if self.player.target is not None:
		i = self.drawHUDListItemLabel('target: ', i, 'top right', self.hudbatch)
		i = self.drawHUDListItemLabel('weapon: ', i, 'top right', self.hudbatch)

		i = 1

		i = self.drawHUDListItemLabel('hyperdrive: ', i, 'bottom right', self.hudbatch)
		i = self.drawHUDListItemLabel('weapon: ', i, 'top right', self.hudbatch)

		
	def drawHUD(self, main_batch):
		if self.player is None:
			return
		# show the player what resources are available

		# layout = pyglet.text.layout.TextLayout(document, width, height, multiline=True)
		# if 

		# self.drawHUDmask(main_batch)

		i = 1
		hudList = {}
		for resource, quantity in list(self.viewpointObject.availableResources.items()):
			hudList[resource] = quantity
		for resource, quantity in list(self.viewpointObject.storagePool.items()):
			if resource in hudList:
				hudList[resource] += self.viewpointObject.storagePool[resource]
			else:
				hudList[resource] = self.viewpointObject.storagePool[resource]

		for availableResource, availableQuantity in list(hudList.items()):
			i = self.drawHUDListItemQuantity(round(availableQuantity,2), i, 'top left', main_batch)

		i = self.drawHUDListItemQuantity( None, i, 'top left', main_batch) # blank line as a separator
		i = 1

		printFreefalling = False
		if self.viewpointObject.freefalling or self.viewpointObject.stepsToFreefall == 0:
			printFreefalling = True

		i = self.drawHUDListItemQuantity( printFreefalling, i, 'bottom left', main_batch)
		i = self.drawHUDListItemQuantity( self.viewpointObject.exemptFromGravity, i, 'bottom left', main_batch)
		if self.player.orbiting is not None:
			i = self.drawHUDListItemQuantity( self.viewpointObject.orbiting.planetName, i, 'bottom left', main_batch)
		i = self.drawHUDListItemQuantity( None, i, 'bottom left', main_batch) # blank line as a separator

		i = self.drawHUDListItemQuantity( self.viewpointObject.isPlayer, i, 'bottom left', main_batch)
		i = self.drawHUDListItemQuantity( self.timestepSize * 3 * 100, i, 'bottom left', main_batch)
		i = self.drawHUDListItemQuantity(self.zoom, i, 'bottom left', main_batch)
		i = self.drawHUDListItemQuantity( self.paused, i, 'bottom left', main_batch)

		i = 1

		if self.player.target is not None:
			i = self.drawHUDListItemQuantity( self.player.target.name, i, 'top right', main_batch)
		i = self.drawHUDListItemQuantity( self.player.selectedWeapon, i, 'top right', main_batch)

		i = 1

		if self.player.hyperdriveDestination is not None:
			i = self.drawHUDListItemQuantity(self.player.hyperdriveDestination.solarSystemName, i, 'bottom right', main_batch)

	def drawInstructions(self, main_batch):
		# prints the help information to the screen.
		i = 1

		i = self.drawHUDListItemLabel('Game Commands', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('unpause / pause and view instructions: p', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('accelerate time: ,', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('decelerate time: .', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('zoom view in: =', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('zoom view out: -', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('look at next actor: ]', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('look at previous actor: [', i, 'top left', main_batch)


		i = self.drawHUDListItemLabel('show HUD: h', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)


		i = self.drawHUDListItemLabel('Piloting', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('accelerate: up arrow', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('setpoint left: left arrow', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('setpoint right: right arrow', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('prograde: w', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('retrograde: s ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('zenith: d', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('nadir: a', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('target: x', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('Astral Navigation', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('jump destination: \\', i, 'top left', main_batch) # backslash is an escape character, so it is doubled here in order to get a backslash character.
		i = self.drawHUDListItemLabel('engage hyperdrive: v', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('galaxy map: m', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('Combat', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('fire guns: space', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('select target: r', i, 'top left', main_batch)


		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('Build Menu', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('toggle build menu: b', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('save to file (ships/playerShip): ctrl + s', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('accept and fly ship: y', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('select part: mouse left click', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('rotate part left: q', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('rotate part right: e', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('mirror part: x', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('put part back in inventory: d', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel(' ', i, 'top left', main_batch)


		i = self.drawHUDListItemLabel('The object of the game is to explore space', i, 'top left', main_batch)

		i = self.drawHUDListItemLabel('Spaceships are made up of many different modules stuck together. Each module provides and consumes resources, and has useful abilities.', i, 'top left', main_batch)
		i = self.drawHUDListItemLabel('Vanquish enemies to smash them into parts. Collect single parts by flying into them. Add new parts to your ship in the build screen.', i, 'top left', main_batch)



	def shootABullet(self, gunModule, actor):
		if gunModule.enabled:
			if gunModule.moduleType == 'cannon 10':
				bulletPosition = [actor.body.position[0] + gunModule.offset[0] + gunModule.barrelHole[0], actor.body.position[1] + gunModule.offset[1] + gunModule.barrelHole[1]]
				bulletPosition = rotate_point(bulletPosition, gunModule.angle, actor.body.position + gunModule.offset )
				bulletPosition = rotate_point(bulletPosition, actor.body.angle, actor.body.position )
				bulletVelocity = [gunModule.muzzleVelocity * math.cos(gunModule.angle + actor.body.angle - 0.5*math.pi) , gunModule.muzzleVelocity * math.sin(gunModule.angle + actor.body.angle- 0.5*math.pi)]
				bullet = Actor('cannonshell 10', [Module('cannonshell 10', (0,0))], bulletPosition ,bulletVelocity,0 ,False)
				self.add(bullet)

	def loadShipIntoBuildMenu(self, actor):
		self.modulesInUse = []
		self.availableModuleListItems = []
		for module in actor.modules:
			self.modulesInUse.append(module)

		for module in self.availableModules:
			self.availableModuleListItems.append(buildMenuItem(module))
			
	def flyShipFromBuildMenu(self):
		playersNewShip = Actor(self.player.name, self.modulesInUse, self.player.body.position, self.player.body.velocity, self.player.body.angle, True )
		self.destroyActor(self.player)
		self.actors.remove(self.player)
		self.add(playersNewShip)

		self.player = self._getPlayer()
		self.viewpointObject = self.player

	def dropModuleIntoBuildArea(self, module, position):
		module.offset = self.antiTransformForBuild(position)
		self.modulesInUse.append(module)

	def buildMenuGraphics(self):
		window.clear()

		main_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		self.drawScreenFill(main_batch)

		# draw the modules the player has assembled 
		for module in self.modulesInUse:
			self.drawModuleForBuild(main_batch, module)

		# draw the inventory list up the left hand side
		i = 1
		for listItem in self.availableModuleListItems:
			self.drawModuleListItem(main_batch, listItem, i)
			i += 1

		# draw the module item the player is dragging, if applicable
		if self.buildDraggingModule is not None:
			self.drawModuleForDrag(main_batch, self.buildDraggingModule, self.mouseCursorPosition)
			pass

		main_batch.draw()


	def drawJumpLinksInMapView(self,solar_system, main_batch):

		for link in solar_system.links:

			aEnd = [int((( solar_system.position[0]) * self.zoom) + resolution_half[0]), int((( solar_system.position[1]) * self.zoom) + resolution_half[1])] # solar_system.position
			bEnd = None
				
			otherSolarSystem = self.galaxy[link]

			bEnd = ([int((( otherSolarSystem.position[0]) * self.zoom) + resolution_half[0]), int(((  otherSolarSystem.position[1]) * self.zoom) + resolution_half[1])]) #otherSolarSystem.position

			if bEnd is not None:
				transformedPoints = transformPolygonForLines([aEnd, bEnd])
				main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',[100,100,100,255]*transformedPoints[0]))

			if self.player.hyperdriveDestination is not None:
				aEnd = [int((( solar_system.position[0]) * self.zoom) + resolution_half[0]), int((( solar_system.position[1]) * self.zoom) + resolution_half[1])] 
				bEnd = ([int((( self.player.hyperdriveDestination.position[0]) * self.zoom) + resolution_half[0]), int(((  self.player.hyperdriveDestination.position[1]) * self.zoom) + resolution_half[1])]) #otherSolarSystem.position

				transformedPoints = transformPolygonForLines([aEnd, bEnd])
				main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',[100,100,100,255]*transformedPoints[0]))



	def drawSolarSystemInMapView(self, solar_system, main_batch):
		#self.drawColorIndicator([20,200,80,255], [rotatedPoint[0], rotatedPoint[1]], int(1*self.zoom), third_batch)
		size = int(5 )
		points = [[-size,size],[size,size],[size,-size],[-size,-size]]
		# newpos = transformForView(position, self.viewpointObject.body.position, self.zoom, resolution )
		nupoints = []
		for point in points:
			nupoints.append([int(((point[0] + solar_system.position[0]) * self.zoom) + resolution_half[0]), int(((point[1] + solar_system.position[1]) * self.zoom) + resolution_half[1])])
			
		renderAConvexPolygon(main_batch, nupoints, resolution_half, self.zoom, resolution,  solar_system.color, solar_system.outlineColor)

	def mapViewGraphics(self):
		window.clear()

		mapBatch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		for key, solar_system in self.galaxy.items():
			self.drawJumpLinksInMapView(solar_system, mapBatch)

		for key, solar_system in self.galaxy.items():
			self.drawSolarSystemInMapView(solar_system, mapBatch)

		mapBatch.draw()

	def graphics(self):
		# window.clear()
		
		if self.viewpointObject is None:
			self.viewpointObject = self.actors[0]

		first_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		self.drawScreenFill(first_batch)
		# self.hudbatch_render()
		# self.hudbatch.draw()

		self.drawBackgroundStars(first_batch)

	
		for attractor in self.attractors:
			if attractor.atmosphere != None:
				self.drawActor(attractor.atmosphere, first_batch)

		first_batch.draw()

		# self.hudbatch_render()


		second_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(1)

		
		# draw the actor's orbits
		if self.showHUD:
			for actor in self.actors:
				if actor.orbit is not None:
					self.drawAPOrbit(second_batch, actor, actor.orbit, actor.orbiting, (100,100,100))

		# # this part draws color dots on where you and your target will be some distance into the future. it's used as a pilot aid for rendezvous.
		if self.player is not None and self.player.target is not None and self.player.target.orbiting is not None:
			# self.futureTime 

			if self.scrollLockTargetFutureTime == '+':
				self.targetfutureTime += 0.5
			elif self.scrollLockTargetFutureTime == '-':
				self.targetfutureTime -= 0.5

			playerFuturetAn = self.player.orbit.tAnAtTime(self.targetfutureTime)
			playerFutureCoordinates = self.player.orbit.cartesianCoordinates(playerFuturetAn)	
			playerFutureCoordinates = [playerFutureCoordinates[0], playerFutureCoordinates[1]]
			playerFutureCoordinates[0] += self.player.orbiting.body.position[0]
			playerFutureCoordinates[1] += self.player.orbiting.body.position[1]

			targetFuturetAn = self.player.target.orbit.tAnAtTime(self.targetfutureTime)
			targetFutureCoordinates = self.player.target.orbit.cartesianCoordinates(targetFuturetAn)	
			targetFutureCoordinates = [targetFutureCoordinates[0], targetFutureCoordinates[1]]
			targetFutureCoordinates[0] += self.player.target.orbiting.body.position[0]
			targetFutureCoordinates[1] += self.player.target.orbiting.body.position[1]

			self.drawColorIndicator([20,200,80,255], playerFutureCoordinates, int(2), second_batch)
			self.drawColorIndicator([20,200,80,255], targetFutureCoordinates, int(2), second_batch)


		second_batch.draw()

		main_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		for illuminator in self.illuminators:
			illuminator.transformedPosition = transformForView( illuminator.position ,self.viewpointObject.body.position, self.zoom, resolution)

		for attractor in self.attractors:
			if attractor.clouds is not None and len(attractor.clouds) > 0:
				for cloud in attractor.clouds:
					self.drawActor(cloud, main_batch)

		main_batch.draw()

		main_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		for attractor in self.attractors:
			self.drawActor(attractor, main_batch)
			
		for actor in self.actors:
			self.drawActor(actor, main_batch)
			

		main_batch.draw()
		
		third_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		if self.showHUD:

			self.hudbatch.draw()

			# print the navcircle
			for line in self.navCircleLines:
				transformedPoints = transformPolygonForLines(line)
				third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[50,50,50,255]*(transformedPoints[0])))

			# blip for actual orientation
			blipLength = (self.navcircleInnerRadius)
			angle = self.viewpointObject.body.angle - 0.5 * math.pi
			start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius+self.navcircleLinesLength) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius+self.navcircleLinesLength) * math.sin(angle)+ (resolution[1]*0.5))
			transformedPoints = transformPolygonForLines([start,end])
			third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[200,40,0,255]*(transformedPoints[0])))

			# blip for setpoint
			blipLength = (self.navcircleInnerRadius)
			angle = self.viewpointObject.setPoint - 0.5 * math.pi
			start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius+self.navcircleLinesLength) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius+self.navcircleLinesLength) * math.sin(angle)+ (resolution[1]*0.5))
			transformedPoints = transformPolygonForLines([start,end])
			third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[250,200,0,255]*(transformedPoints[0])))

			if self.viewpointObject.freefalling and self.viewpointObject.orbit is not None:
				if not math.isnan( self.viewpointObject.prograde):
					# blip for prograde
					blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
					angle = self.viewpointObject.prograde
					start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
					end = ((self.navcircleInnerRadius) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius) * math.sin(angle)+ (resolution[1]*0.5))
					transformedPoints = transformPolygonForLines([start,end])
					third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[0,200,20,255]*(transformedPoints[0])))

					# blip for retrograde
					blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
					angle = self.viewpointObject.retrograde
					start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
					end = ((self.navcircleInnerRadius) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius) * math.sin(angle)+ (resolution[1]*0.5))
					transformedPoints = transformPolygonForLines([start,end])
					third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[0,100,10,255]*(transformedPoints[0])))

			# blip for nadir
			blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
			angle = self.viewpointObject.nadir
			start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius) * math.sin(angle)+ (resolution[1]*0.5))
			transformedPoints = transformPolygonForLines([start,end])
			third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[0,100,200,255]*(transformedPoints[0])))

			# blip for zenith
			blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
			angle = self.viewpointObject.zenith
			start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius) * math.cos(angle)+ (resolution[0]*0.5),- (self.navcircleInnerRadius) * math.sin(angle)+ (resolution[1]*0.5))
			transformedPoints = transformPolygonForLines([start,end])
			third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[0,50,100,255]*(transformedPoints[0])))

			# blip for target
			if self.player is not None:
				if self.player.target is not None:
					blipLength = (self.navcircleInnerRadius+self.navcircleLinesLength)
					angle = math.atan2(self.player.target.body.position[1] - self.player.body.position[1],(self.player.target.body.position[0] - self.player.body.position[0]))
					start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
					end = ((self.navcircleInnerRadius+self.navcircleLinesLength*2) * math.cos(angle)+ (resolution[0]*0.5),- ((self.navcircleInnerRadius+self.navcircleLinesLength*2)) * math.sin(angle)+ (resolution[1]*0.5))
					transformedPoints = transformPolygonForLines([start,end])
					third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[250,250,250,255]*(transformedPoints[0])))

			# blip for hyperdrive destination
			if self.player is not None:
				if self.player.hyperdriveDestination is not None:
					blipLength = (self.navcircleInnerRadius+self.navcircleLinesLength)
					angle = math.atan2(self.player.hyperdriveDestination.position[1] - self.currentSystem.position[1],(self.player.hyperdriveDestination.position[0] - self.currentSystem.position[0]))
					start = ((blipLength * math.cos(angle)) + (resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( resolution[1] * 0.5) )
					end = ((self.navcircleInnerRadius+self.navcircleLinesLength*2) * math.cos(angle)+ (resolution[0]*0.5),- ((self.navcircleInnerRadius+self.navcircleLinesLength*2)) * math.sin(angle)+ (resolution[1]*0.5))
					transformedPoints = transformPolygonForLines([start,end])
					third_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[160,0,255,255]*(transformedPoints[0])))

		if self.showHUD:
			self.drawHUD(third_batch)

		# this part draws color dots on all the modules so you can tell if they're active or enabled or off. But I think it looks ugly, so I disabled it.
		if self.showHUD and False:
			for module in self.player.modules:
				if module.enabled:
					rotatedPoint = rotate_point([self.player.body.position[0] +module.offset[0], self.player.body.position[1]+module.offset[1]], self.player.body.angle, self.player.body.position)
					if module.active:
						self.drawColorIndicator([20,200,80,255], [rotatedPoint[0], rotatedPoint[1]], int(1*self.zoom), third_batch)
					else:
						self.drawColorIndicator([255,150,50,255], [rotatedPoint[0], rotatedPoint[1]], int(1*self.zoom), third_batch)



		if self.paused:
			self.drawInstructions(third_batch)

		third_batch.draw()

	def generateBackgroundStars(self):
		self.backgroundStars = []
		self.backgroundStarColorstream = []
		n_stars = 350
		for i in range(n_stars):
			position = [random.randint(0,resolution[0]), random.randint(0,resolution[1])]
			color = [155,155,155,255]
			# color[0] += random.randint(0,100)
			# color[2] += random.randint(0,100)

			brightness = random.randint(0,10)

			# color[0] = int(color[0] * (brightness/10))
			# color[1] = int(color[1] * (brightness/10))
			# color[2] = int(color[2] * (brightness/10))

			bv_key = random.randint(0,173)
			bv_key -= 33
			bv_key = bv_key/100
			decodedcolor = bv2rgb(bv_key)

			# print(bv_key)
			# print(decodedcolor)


			color[0] = int(decodedcolor[0] * (brightness/10) * 255)
			color[1] = int(decodedcolor[1] * (brightness/10) * 255)
			color[2] = int(decodedcolor[2] * (brightness/10) * 255)


			# self.backgroundStars.append(BackgroundStar(position, color, random.randint(1,2)))
			self.backgroundStars.append(position[0])
			self.backgroundStars.append(position[1])
			self.backgroundStarColorstream.append(color[0])
			self.backgroundStarColorstream.append(color[1])
			self.backgroundStarColorstream.append(color[2])
			self.backgroundStarColorstream.append(255)

	def drawBackgroundStars(self, main_batch):
		# print(self.backgroundStars)
		main_batch.add( 350, pyglet.gl.GL_POINTS, None, ('v2i', self.backgroundStars), ('c4B',self.backgroundStarColorstream))
		
	def hyperspaceJump(self, actor) :

		# print(self.attractors)

		# teleports the player across space, by removing all the other actors and stuff and replacing them with new stuff.
		if actor.hyperdriveDestination is None:
			return
		if actor.isPlayer:
			previousSystem = self.currentSystem
			self.currentSystem = actor.hyperdriveDestination

			for otherActor in self.actors:
				# print(otherActor)
				if otherActor is not actor:
					self.destroyActor(otherActor)
			for attractor in self.attractors:
				# print(attractor)
				self.destroyAttractor(attractor)

			self.actors = []
			self.actors.append(actor)

			self.attractors = []

			for thing in self.currentSystem.contents:
				self.add(thing)

			playerInsertionAngle = math.atan2(self.currentSystem.position[1] - previousSystem.position[1], self.currentSystem.position[0] - previousSystem.position[1])
			actor.body.position = [self.currentSystem.hyperspaceThreshold * math.cos(playerInsertionAngle), self.currentSystem.hyperspaceThreshold * math.sin(playerInsertionAngle)]

			self.viewpointObject = self._getPlayer()
		else:
			self.destroyActor(actor)
			self.actors.remove(actor)


		self.generateBackgroundStars()

	def dock(self, actor):

		# find the station's docking ring
		for module in actor.target.modules:
			if module.moduleType == 'docking port':
				
				dockingPoint = rotate_point(-module.effect.offset,module.angle)  # orient the polygon according to the body's current direction in space.
				dockingPoint = rotate_point(rotatedPoints, actor.target.body.angle, [-(module.offset[0] + module.effect.offset[0]), -(module.offset[1] + module.effect.offset[1])])
				dockingPoint = dockingPoint + actor.target.body.position + module.offset + module.effect.offset
				
				# self.drawColorIndicator()

				# go through all the actors points and see if any of them are in the docking ring.
				# allActorPoints = []
				closeEnoughToDock = False
				for module in actor.modules:
					for point in module.points:
						# pass
						# allActorPoints.append(point)
						# for box in boxes:
						# 	if pointInRect(point, box):
						# 		closeEnoughToDock = True
						# rotatedPoints = module.points
						rotatedPoint = rotate_point(point,module.angle)  # orient the polygon according to the body's current direction in space.
						rotatedPoint = rotate_point(rotatedPoint, actor.body.angle, [-module.offset[0], -module.offset[1]])
						rotatedPoint = rotatedPoint + actor.body.position + module.offset

						if mag(rotatedPoint - dockingPoint) < 50:
							# succcessss
							pass


	def step(self):
		if self.buildMenu:
			self.buildMenuGraphics()
		elif self.mapView:
			self.mapViewGraphics()
		else:
			
			if not self.paused:
				# print(self.illuminators)
				self.illuminators = []
				self.physics()

			self.player = self._getPlayer()


			if self.viewpointObject is None:
				self.viewpointObject = self.player

			self.graphics()

	def addSolarSystem(self, solar_system):
		self.galaxy[solar_system.solarSystemName] = solar_system

	def setup(self):
		self.createHUDNavCircle()

		self.generateBackgroundStars()



		for module in shipyard('sport parts'):
			self.availableModules.append(copy.deepcopy(module))

		for solar_system in solarSystems:
			self.addSolarSystem(solar_system)
			
		self.currentSystem = self.galaxy['Sol III']
		print(self.currentSystem)

		for thing in self.currentSystem.contents:
			self.add(thing)

		loaddedbrige =  shipyard('starbridge')

		rockeyt =  shipyard('rocket_1')

		ida_frigate_instance = Actor('player ida_frigate', loaddedbrige,(-800000, -1), [1,50000], 0.6 * math.pi, True)

		fojgesogj = Actor('psefse', rockeyt,(-800000, -100), [1,50000], 0.6 * math.pi)

		ida_frigate_instance.maneuverQueue.append(Maneuver('transfer',earth, moon))

		self.add(ida_frigate_instance)
		self.add(fojgesogj)

		self.viewpointObject = self._getPlayer()

	def start(self):
		self.setup()

def saveShipFromBuildMenu():
	dill.dump(Nirn.modulesInUse,open('ships/playerShip', 'wb'))

def loadShipIntoBuildMenu():
	Nirn.modulesInUse = dill.load(open('ships/playerShip', 'rb'))

def save():
	newSave = SavedGame(Nirn.actors, Nirn.currentSystem, Nirn.availableModules, Nirn.player)
	dill.dump(newSave,open('save', 'wb'))


def load():
	
	newSave=dill.load(open('save', 'rb'))

	# exit the current system and load the one from the file
	Nirn.hyperspaceJump(newSave.currentSystem.solarSystemName, Nirn.player)
	Nirn.actors = newSave.actors
	Nirn.availableModules = newSave.availableModules
	Nirn.player = Nirn._getPlayer()
	Nirn.viewpointObject = Nirn.player

Nirn = World()

@window.event
def on_key_press(symbol, modifiers):
	if symbol == key.ESCAPE:
		exit()
	elif symbol == key.LEFT:
		if Nirn.player is not None:
			Nirn.player.keyStates['left'] = True
			Nirn.player.keyStates['face direction'] = None
	elif symbol == key.RIGHT:
		if Nirn.player is not None:
			Nirn.player.keyStates['right'] = True
			Nirn.player.keyStates['face direction'] = None
	elif symbol == key.UP:
		if Nirn.player is not None:
			Nirn.player.keyStates['up'] = True
	elif symbol == key.P:
		Nirn.paused = not Nirn.paused
	elif symbol == key.EQUAL:
		Nirn.zoom = round_to_n(Nirn.zoom + (Nirn.zoom * 0.5), 2)
	elif symbol == key.MINUS:
		Nirn.zoom = round_to_n(Nirn.zoom - (Nirn.zoom * 0.5), 2)
	elif symbol == key.H:
		Nirn.showHUD = not Nirn.showHUD

		# if Nirn.showHUD:
		# 	# resolution = (1400,900)
		# 	# resolution_half = (1400/2,900/2)
		Nirn.hudbatch_render()
		# else:
		# 	pass
			# resolution = (1600,900)
			# resolution_half = (1600/2,900/2)
	elif symbol == key.COMMA:
		Nirn.timestepSize = round_to_n(Nirn.timestepSize + (Nirn.timestepSize * 0.5), 2)
	elif symbol == key.PERIOD:
		Nirn.timestepSize = round_to_n(Nirn.timestepSize - (Nirn.timestepSize * 0.5), 2)
	elif symbol == key.B:
		if Nirn.buildMenu:
			Nirn.buildMenu = False
		else:
			Nirn.buildMenu = True
			Nirn.paused = True
			Nirn.loadShipIntoBuildMenu(Nirn.player)
	elif symbol == key.Y:
		if Nirn.buildMenu:
			Nirn.flyShipFromBuildMenu()
			Nirn.buildMenu = False
	elif symbol == key.D:
		if Nirn.buildMenu:
			if Nirn.buildDraggingModule is not None:
				# put buildDraggingModule back in to the player inventory
				Nirn.availableModuleListItems.append(buildMenuItem(Nirn.buildDraggingModule))
				Nirn.availableModules.append(Nirn.buildDraggingModule)
				Nirn.buildDraggingModule = None
		else:
			if Nirn.player is not None:
				Nirn.player.keyStates['face direction'] = 'zenith'
			
	elif symbol == key.E:
		if Nirn.buildMenu:
			if Nirn.buildDraggingModule is not None:
				Nirn.buildDraggingModule.angle += (1/32) * math.pi
	elif symbol == key.Q:
		if Nirn.buildMenu:
			if Nirn.buildDraggingModule is not None:
				Nirn.buildDraggingModule.angle -= (1/32) * math.pi
	elif symbol == key.BRACKETRIGHT:
		Nirn.viewpointObjectIndex += 1
		if Nirn.viewpointObjectIndex >= len(Nirn.actors):
			Nirn.viewpointObjectIndex = 0
		Nirn.viewpointObject = Nirn.actors[Nirn.viewpointObjectIndex]
	elif symbol == key.BRACKETLEFT:
		Nirn.viewpointObjectIndex -= 1
		if Nirn.viewpointObjectIndex < 0:
			Nirn.viewpointObjectIndex = len(Nirn.actors)-1
		Nirn.viewpointObject = Nirn.actors[Nirn.viewpointObjectIndex]
	elif symbol == key.SPACE:
		if Nirn.player is not None:
			Nirn.player.keyStates['Fire'] = True
	elif symbol == key.S:
		if modifiers & key.MOD_CTRL:
			if Nirn.buildMenu:
				saveShipFromBuildMenu()
		else:
			if Nirn.player is not None:
				Nirn.player.keyStates['face direction'] = 'retrograde'
	elif symbol == key.W:
		if Nirn.player is not None:
			Nirn.player.keyStates['face direction'] = 'prograde'
	elif symbol == key.A:
		if Nirn.player is not None:
			Nirn.player.keyStates['face direction'] = 'nadir'
	elif symbol == key.T:
		if Nirn.player is not None:
			Nirn.player.autoPilotActive = not Nirn.player.autoPilotActive
	elif symbol == key.L:
		if not Nirn.buildMenu and not Nirn.mapView:
			Nirn.player.keyStates['strafe right'] = True
		if modifiers & key.MOD_CTRL:
			if Nirn.buildMenu:
				loadShipIntoBuildMenu()
	elif symbol == key.V:
		if Nirn.player is not None:
			Nirn.player.keyStates['hyperdrive engage'] = True
	elif symbol == key.M:
		Nirn.mapView = not Nirn.mapView
	elif symbol == key.R:
		if Nirn.playerTargetIndex == None:
			Nirn.playerTargetIndex = 0
		else:
			Nirn.playerTargetIndex += 1
			if Nirn.playerTargetIndex > len(Nirn.actors) - 1:
				Nirn.playerTargetIndex = None

		if Nirn.playerTargetIndex is not None:
			if Nirn.player is not None:
				if Nirn.actors[Nirn.playerTargetIndex] is Nirn.player:
					Nirn.playerTargetIndex += 1
				Nirn.player.target = Nirn.actors[Nirn.playerTargetIndex]
		else:
			if Nirn.player is not None:
				Nirn.player.target = None

	elif symbol == key.BACKSLASH:
		if Nirn.player is not None:
			if Nirn.hyperjumpDestinationIndex == None:
				Nirn.hyperjumpDestinationIndex = 0
			else:
				Nirn.hyperjumpDestinationIndex += 1
				if Nirn.hyperjumpDestinationIndex > len(Nirn.currentSystem.links) - 1:
					Nirn.hyperjumpDestinationIndex = None

			if Nirn.hyperjumpDestinationIndex is not None:
				if Nirn.currentSystem.links[Nirn.hyperjumpDestinationIndex]  is Nirn.currentSystem:
					Nirn.hyperjumpDestinationIndex += 1
				Nirn.player.hyperdriveDestination =Nirn.galaxy[Nirn.currentSystem.links[Nirn.hyperjumpDestinationIndex]]
			else:
				Nirn.player.hyperdriveDestination = None
	elif symbol == key.X:
		if Nirn.buildMenu:
			Nirn.buildDraggingModule.points = mirror_polygon(Nirn.buildDraggingModule.points)
		else:
			Nirn.player.keyStates['face direction'] = 'target'

	elif symbol == key.I:
		Nirn.player.keyStates['strafe forward'] = True
	elif symbol == key.J:
		Nirn.player.keyStates['strafe left'] = True
	elif symbol == key.K:
		Nirn.player.keyStates['strafe back'] = True
	# elif symbol == key.L:
		# Nirn.player.keyStates['strafe right'] = True
	elif symbol == key.NUM_MULTIPLY:
		# Nirn.targetfutureTime += 0.5
		Nirn.scrollLockTargetFutureTime = '+'
	elif symbol == key.NUM_DIVIDE:
		# Nirn.targetfutureTime -= 0.5
		Nirn.scrollLockTargetFutureTime = '-'
	elif symbol == key.G:
		Nirn.dock()
		# dock the ship to a docking port that it is near

@window.event
def on_key_release(symbol, modifiers):
	if symbol == key.LEFT:
		if Nirn.player is not None:
			Nirn.player.keyStates['left'] = False
	elif symbol == key.RIGHT:
		if Nirn.player is not None:
			Nirn.player.keyStates['right'] = False
	elif symbol == key.UP:
		if Nirn.player is not None:
			Nirn.player.keyStates['up'] = False		
	elif symbol == key.V:
		if Nirn.player is not None:
			Nirn.player.keyStates['hyperdrive engage'] = False	
	elif symbol == key.I:
		if Nirn.player is not None:
			Nirn.player.keyStates['strafe forward'] = False
	elif symbol == key.J:
		if Nirn.player is not None:
			Nirn.player.keyStates['strafe left'] = False
	elif symbol == key.K:
		if Nirn.player is not None:
			Nirn.player.keyStates['strafe back'] = False
	elif symbol == key.L:
		if Nirn.player is not None:
			Nirn.player.keyStates['strafe right'] = False	
	elif symbol == key.NUM_MULTIPLY:
		Nirn.scrollLockTargetFutureTime = None
	elif symbol == key.NUM_DIVIDE:
		Nirn.scrollLockTargetFutureTime = None


def stepWithBatch(dt):
	pass

@window.event()
def on_mouse_press(x, y, button, modifiers):
	Nirn.mouseCursorPosition =[x,y]
	if Nirn.buildMenu:
		if mouse.LEFT is button:
			Nirn.buildDraggingModule = Nirn.getModuleFromCursorPosition(Nirn.mouseCursorPosition)

@window.event()
def on_mouse_release(x, y, button, modifiers):
	Nirn.mouseCursorPosition = [x,y]
	if Nirn.buildMenu:
		if mouse.LEFT is button:
			if Nirn.buildDraggingModule is not None:
				Nirn.dropModuleIntoBuildArea(Nirn.buildDraggingModule, Nirn.mouseCursorPosition)
				Nirn.buildDraggingModule = None

@window.event()
def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
    Nirn.mouseCursorPosition =[x,y]

@window.event()
def on_draw():
	Nirn.step()

def hello():
	Nirn.start()

	pyglet.clock.schedule_interval(stepWithBatch, 0.01)
	pyglet.app.run()

hello()