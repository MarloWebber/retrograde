import math, sys, random

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

import pyglet
from pyglet.gl import *
from pyglet.window import Window
from pyglet.window import key
from pyglet.window import mouse

import copy

import cProfile

resolution = (1280,780)
window = pyglet.window.Window(width=1280, height=780)
label = pyglet.text.Label('Abc', x=5, y=5)

white = [255]*4

def mag(x):
	return numpy.sqrt(x.dot(x))

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

	# print([[leastX, leastY], [mostX, mostY]])

	return [[leastX, leastY], [mostX, mostY]]

def pointInRect(point,rect):
    if (rect[0][0] < point[0] and point[0] < rect[1][0]):
        if (rect[0][1] < point[1] and point[1] < rect[1][1]):
            return True
    return False

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
	return density * 4/3 * math.pi * radius**3

def mass_of_a_circle(density, radius):
	return density * 2 * math.pi * radius**2

def semiMinorAxis(a, e):
	return a * math.sqrt(1 - (math.pow(e,2))) # https://math.stackexchange.com/questions/1259945/calculating-semi-minor-axis-of-an-ellipse

def averageOfList(TheList):
    avg = sum(TheList) / len(TheList)
    return avg

def centroidOfPolygon(polygon):
	xValues = []
	yValues = []
	# print(polygon)
	for point in polygon:
		xValues.append(point[0])
		yValues.append(point[1])

	return [averageOfList(xValues), averageOfList(yValues)]

# def point_inside_polygon(x, y, poly, include_edges=True):
#     '''
#     Test if point (x,y) is inside polygon poly.

#     poly is N-vertices polygon defined as 
#     [(x1,y1),...,(xN,yN)] or [(x1,y1),...,(xN,yN),(x1,y1)]
#     (function works fine in both cases)

#     Geometrical idea: point is inside polygon if horisontal beam
#     to the right from point crosses polygon even number of times. 
#     Works fine for non-convex polygons.
#     '''
#     n = len(poly)
#     inside = False

#     p1x, p1y = poly[0]
#     for i in range(1, n + 1):
#         p2x, p2y = poly[i % n]
#         if p1y == p2y:
#             if y == p1y:
#                 if min(p1x, p2x) <= x <= max(p1x, p2x):
#                     # point is on horisontal edge
#                     inside = include_edges
#                     break
#                 elif x < min(p1x, p2x):  # point is to the left from current edge
#                     inside = not inside
#         else:  # p1y!= p2y
#             if min(p1y, p2y) <= y <= max(p1y, p2y):
#                 xinters = (y - p1y) * (p2x - p1x) / float(p2y - p1y) + p1x

#                 if x == xinters:  # point is right on the edge
#                     inside = include_edges
#                     break

#                 if x < xinters:  # point is to the left from current edge
#                     inside = not inside

#         p1x, p1y = p2x, p2y

#     return inside

# def ccw(A,B,C):
#     return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x)

# # Return true if line segments AB and CD intersect
# def intersect(A,B,C,D):
#     return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

# def point_in_polygon(pt, poly, inf):
#     result = False
#     for i in range(len(poly)-1):
#         if intersect((poly[i][0], poly[i][1]), ( poly[i+1][0], poly[i+1][1]), (pt[0], pt[1]), (inf, pt[1])):
#             result = not result
#     if intersect((poly[-1][0], poly[-1][1]), (poly[0][0], poly[0][1]), (pt[0], pt[1]), (inf, pt[1])):
#         result = not result
#     return result

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
	# if seviewpointObject == None:
	# 	return position
	# else:
	transformedPosition = position - viewpointObjectPosition # offset everything by the position of the viewpointObject, so the viewpoint is at 0,0
	transformedPosition = transformedPosition * zoom  # shrink or expand everything around the 0,0 point
	transformedPosition[0] = int(transformedPosition[0] + 0.5 * resolution[0]) # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
	transformedPosition[1] = int(-transformedPosition[1] + 0.5 * resolution[1]) # add half the height. and invert it so that it's the right way up in opengl.
	return transformedPosition

def isPointIlluminated(point, color, illuminators,viewpointObjectPosition, zoom, resolution):
	resultColor = [color[0],color[1],color[2],255]
	for illuminator in illuminators:
		transformedPosition = transformForView( illuminator.position ,viewpointObjectPosition, zoom, resolution)
		distance =numpy.array([point[0] - transformedPosition[0],point[1] - transformedPosition[1]])
		

		if distance[0] < illuminator.radius and distance[1] < illuminator.radius :
			
			magnitude = mag(distance )
			magnitude = magnitude/zoom
			print(magnitude)
			if magnitude == 0:
				dispersion = 1
			else:
				dispersion = 1/magnitude
			# print(dispersion)
			resultColor[0] += int(dispersion * illuminator.color[0])
			if resultColor[0] > 255: resultColor[0] = 255
			resultColor[1] += int(dispersion * illuminator.color[1])
			if resultColor[1] > 255: resultColor[1] = 255
			resultColor[2] += int(dispersion * illuminator.color[2])
			if resultColor[2] > 255: resultColor[2] = 255
	return resultColor



def transformPolygonForLinesWithIlluminators(polygon, color, illuminators, viewpointObjectPosition, zoom, resolution):
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
		# shred one of the ingame annular atmospheres into a ribbon of triangles, and add a color gradient from the inner edge to the outer.

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
			nn += 2

			colorstream.append(inside_color[0])
			colorstream.append(inside_color[1])
			colorstream.append(inside_color[2])
			colorstream.append(inside_color[3])

			colorstream.append(outside_color[0])
			colorstream.append(outside_color[1])
			colorstream.append(outside_color[2])
			colorstream.append(outside_color[3])

		bitstream.append(outside_verts[index][0])
		bitstream.append(outside_verts[index][1])
		
		colorstream.append(outside_color[0])
		colorstream.append(outside_color[1])
		colorstream.append(outside_color[2])
		colorstream.append(outside_color[3])

		nn += 1

		return [nn, bitstream, colorstream]




def renderAConvexPolygon(batch, polygon, viewpointObjectPosition, zoom, resolution, color, outlineColor=None, illuminators=None):
	# nsfe = int(1.5*len(polygon) + 5)
	# print(nsfe)
	transformedPoints = transformPolygonForTriangleFan(polygon)
	# print(transformedPoints[0])
	batch.add(transformedPoints[0], pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',transformedPoints[1]), ('c4B',color*transformedPoints[0]))

	if outlineColor is not None:
		# transformedPoints = transformPolygonForLines(polygon)
		if illuminators is not None:
			transformedPoints = transformPolygonForLinesWithIlluminators(polygon, outlineColor, illuminators, viewpointObjectPosition, zoom, resolution)
			batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',transformedPoints[2]))
		else:
			transformedPoints = transformPolygonForLines(polygon)
			batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i',transformedPoints[1]), ('c4B',outlineColor*transformedPoints[0]))

class Atmosphere():
	def __init__(self,radius, planetPosition):
		self.height = 15000
		self.bottomDensity = 1
		self.topDensity = 0
		self.color = [70,145,220,255]
		self.outerColor = [0,0,0,255]
		self.outlineColor = [50,125,200,255]

		# each atmosphere layer is an annulus.
		# it is rendered as a triangle strip with a gradient between the inner and outer edges.
		self.n_points = 100
		self.innerPoints = make_circle(radius, 100)
		self.outerPoints = make_circle(radius+self.height, 100)

		# self.rendering = unwrapForGradient(n_points, self.innerPoints, self.outerPoints, self.color, self.outerColor)

		self.mass = ( (self.bottomDensity + self.topDensity) / 2 ) * area_of_annulus(radius+self.height, radius)
		# inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		# self.body = pymunk.Body(self.mass, inertia)
		# self.body.position = planetPosition

class ModuleEffect(): # a ModuleEffect is just a polygon that is visually displayed when a module is active.
	def __init__(self, offset=[0,0]):
		self.offset = offset
		self.radius = 3
		self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]

# class Bullet(): # a physical projectile that is too simple to be an actor.
# 	def __init__(self, position, velocity, bulletType)
# 		self.position = position
# 		self.velocity = velocity

# 		if bulletType == 'Cannon 10'
# 			self.

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
				'electricity': 0.1,
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
			self.radius = 5
			self.points = [[-self.radius, -self.radius], [-self.radius, self.radius], [self.radius,self.radius], [self.radius, -self.radius]]
			self.color = [75,10,10,255]
			self.outlineColor = [200,70,70,255]

		elif self.moduleType is 'engine':
			self.mass = 1
			self.resources = {
				'thrust': 100,
				'fuel': -0.1,
				'electricity':1
			}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size*2], [-size, size*2], [size,size*2], [size, -size*2]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]

			self.illuminatorOffset = [0,10]

			self.effect = ModuleEffect([0,0])

		elif self.moduleType is 'RCS':
			self.mass = 0.2
			self.resources = {
				'torque': 5,
				'electricity': -0.2
			}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size], [-size, size], [size,size], [size, -size]]
			self.color = [120,100,100,255]
			self.outlineColor = [170,150,150,255]

			self.momentArm = self.radius

		elif self.moduleType is 'spar 10':
			self.mass = 1
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size, -size*10], [-size, size*10], [size,size*10], [size, -size*10]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]

			self.momentArm = self.radius

		elif self.moduleType is 'box 10':
			self.mass = 10
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*10, -size*10], [-size*10, size*10], [size*10,size*10], [size*10, -size*10]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]

			self.momentArm = self.radius

		elif self.moduleType is 'box 100':
			self.mass = 100
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*100, -size*100], [-size*100, size*100], [size*100,size*100], [size*100, -size*100]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]

			self.momentArm = self.radius

		elif self.moduleType is 'spar 100':
			self.mass = 1
			self.resources = {}
			self.stores = {}
			self.initialStores = {}
			self.radius = 5
			size = self.radius
			self.points = [[-size*10, -size*100], [-size*10, size*100], [size*10,size*100], [size*10, -size*100]]
			self.color = [50,50,50,255]
			self.outlineColor = [100,100,100,255]

			self.momentArm = self.radius

		elif self.moduleType is 'cannonshell 10':
			self.mass = 1
			self.active = True
			self.resources = {}
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

			self.lifetime = 100 # the shell lasts for 1000 somethings and then explodes.


		elif self.moduleType is 'cannon 10':
			self.mass = 1
			self.active = False
			self.resources = {
				# 'cannonshell 10': 1
				'electricity':0.001
			}
			self.stores = {
				'cannonshell 10': 100
			}
			self.initialStores = {
				'cannonshell 10': 100
			}
			self.radius = 5
			self.points = [[-0.5*self.radius, -self.radius], [-0.5*self.radius, self.radius], [0.5*self.radius,self.radius], [0.5*self.radius, -self.radius]]
			self.color = [30,30,30,255]
			self.outlineColor = [100,100,30,255]

			self.barrelHole = [0,-self.radius + 1.5]
			self.muzzleVelocity = 500000
			self.cooldownTime = 100
			self.cooldownValue = 0


dinghy = [Module('generator',[0,0]), Module('engine',[0,8]), Module('RCS',[0,-10]) ]
lothar = [Module('generator',[0,0]), Module('engine',[-13,8], 0.6/math.pi), Module('engine',[13,8],-0.6/math.pi), Module('RCS',[-13,-10]), Module('RCS',[13,-10]) , Module('cannon 10',[0,-10]) ]
boldang = [Module('spar 10',[0,-100], (0.5* math.pi)), Module('box 10',[0,0])]
bigmolly = [Module('box 100',[0,0]), Module('spar 100',[1000,0], 0.5 * math.pi),Module('box 100',[-1000,0]),Module('box 100',[2000,0]), Module('box 100',[-2000,0]),  Module('box 100',[3000,0])]

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
		
		inertia = pymunk.moment_for_poly(self.mass, self.points, (0,0))
		self.body = pymunk.Body(self.mass, inertia)
		self.body.position = position
		self.body.velocity = velocity
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
		self.stepsToFreefall = 1
		self.decompEnergy = 5000000
		self.desiredAngle = 0
		self.exemptFromGravity = False
		self.body.angle = angle

	def leaveFreefall(self, stepsToFreefall=1):
		self.stepsToFreefall = stepsToFreefall
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



			if module.enabled:
				if module.active:
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
		self.planetName = planetType # just set the planetName to something easy for now. # the name of the individual instance of this planet type.
		
		if planetType == 'earth':
			self.radius = 320000
			self.density = 1
			self.friction = 0.9
			self.color = [180,170,145,255]
			self.outlineColor = [200,190,155,255]
			self.atmosphere = Atmosphere(self.radius, position)

		elif planetType == 'moon':
			self.radius = 80000
			self.density = 1
			self.friction = 0.9
			self.color = [45,45,45,255]
			self.outlineColor = [145,145,145,255]
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
		self.APBody = APBody(self.planetName, self.mass * gravitationalConstant * 0.163, self.radius)

class buildMenuItem():
	# a buildMenuItem is literally the tiles in the build menu you can click and drag to add modules to your ship.
	def __init__(self, module, boundingRectangle=((0,0),(1,1))):
		self.module = module
		self.quantity = 1
		self.boundingRectangle = boundingRectangle

class Illuminator():
	def __init__(self, position):
		self.radius = 1000
		self.color = [255,255,0]
		self.intensity = 1
		self.position = position #(100, -320030)

class World():
	def __init__(self):
		# pygame.init()
		# self.clock = pygame.time.Clock() # the pygame clock is NOT the same as the simulation clock.
		self.time = 0 # the number of timesteps that have passed in-universe. used for physics and orbital calculations.
		self.space = pymunk.Space()
		self.space.gravity = (0.0, 0.0)
		self.gravitationalConstant = 0.03
		self.dragCoefficient = 0.005 # atmospheric drag
		self.actors = []
		self.attractors = []
		self.resolution = resolution

		# self.screen = pygame.display.set_mode(self.resolution)
		# self.draw_options = pymunk.pygame_util.DrawOptions(self.screen)
		self.ch = self.space.add_collision_handler(0, 0)
		# self.ch.data["surface"] = self.screen
		self.ch.post_solve = self.handle_collision

		# self.window = pyglet.window.Window()

		self.viewpointObject = None
		self.viewpointObjectIndex = 0
		self.player = None
		self.zoom = 1 # the actual applied zoom number.
		self.pan = [0,0]
		self.rotate = 0
		self.timestepSize = 0.2/60.0 #1.0/60.0
		# pygame.key.set_repeat(50,50) # holding a key down repeats the instruction. https://www.pygame.org/docs/ref/key.html
		# self.font = pygame.font.SysFont('dejavusans', 12)
		self.showHUD = False
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

	def getModuleFromCursorPosition(self, cursorPositionRaw):
		cursorPosition = cursorPositionRaw
		# cursorPosition[1] = -cursorPosition[1] + 0.5 * resolution[0]
		smellyCursor = [0,0]
		smellyCursor[0] = cursorPositionRaw[0]
		smellyCursor[1] = -cursorPositionRaw[1] +  resolution[1]
		for listItem in self.availableModuleListItems:
			# if listItem.boundingRectangle.collidepoint(cursorPosition):

			# xMost = 0
			# for point in listItem.boundingRectangle:
			# 	if point[0] > xMost:
			# 		xMost = point[0]

			# if point_in_polygon(cursorPosition, listItem.boundingRectangle,xMost ):

			# print(smellyCursor)
			# print(self.transformForBuild(cursorPosition))
			# print(listItem.boundingRectangle)
			if pointInRect(smellyCursor, listItem.boundingRectangle):
				# print('click a list item')
				self.availableModuleListItems.remove(listItem)
				return listItem.module

		for module in self.modulesInUse:
			transformedPoints = []
			for point in module.points:
				transformedPoint = [0,0]
				transformedPoint[0] = (point[0] + module.offset[0])
				transformedPoint[1] = (point[1] + module.offset[1])
				transformedPoints.append(transformedPoint)

			boundingBox = boundPolygon(transformedPoints)
			print(self.antiTransformForBuild(cursorPositionRaw))
			print(boundingBox)
			if pointInRect( self.antiTransformForBuild(cursorPositionRaw) , boundingBox):
				# print('click a module in use')
				self.modulesInUse.remove(module)
				return module

	def add(self, thing):  
		self.space.add(thing.body, thing.shape)
		if thing.__class__ == Actor:
			self.actors.append(thing)
		elif thing.__class__ == Attractor:
			self.attractors.append(thing)

		self.player = self._getPlayer()
		self.viewpointObject = self.player

	def destroyActor(self, actor):
		self.space.remove(actor.body)
		self.space.remove(actor.shape)
		if actor in self.actors:
				self.actors.remove(actor)

	def decomposeActor(self, actor, modules):
		listLength = len(actor.modules)
		if listLength == 1:
			return # the actor is already fully decomposed, destroy it if you want
		else:
		
			# create a new actor, minus the module
			for index, module in enumerate(actor.modules):

				# create the module on it's own as a new actor
				fragmentPosition = [actor.body.position[0] + (module.offset[0] * math.cos(actor.body.angle)), actor.body.position[1] +  (module.offset[1] * math.sin(actor.body.angle))]

				if actor.isPlayer and index == listLength-1:
					self.add(Actor(actor.name + ' fragment', [module], fragmentPosition, actor.body.velocity, True))
					self.player = self._getPlayer()
					self.viewpointObject = self.player
				else:
					self.add(Actor(actor.name + ' fragment', [module], fragmentPosition, actor.body.velocity, False))

			# remove the actor from the space.
			self.destroyActor(actor)


			# if not actor.modules:
			# 	pass
			# else:
			# 	self.add( Actor(actor.name, actor.modules, actor.body.position, actor.body.velocity, actor.isPlayer ) ) # add the remaining parts of the original actor back into the space

	def physics(self):
		for actor in self.actors:
			actor.exemptFromGravity = False

		# iterate the physics engine. This is done first, so that changes can be reacted to (i.e. to exempt things from gravity if they collide with a planet, to prevent sinking).
		self.space.step(self.timestepSize)
		self.time += self.timestepSize

		for actor in self.actors:

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


			actor.doResources()

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


						#drag = coefficient * (density * velocity**2 / 2) * reference area
						density = (attractor.atmosphere.topDensity + (naturalDepth * attractor.atmosphere.bottomDensity))**1.5 # not quite squared, but still exponential

						# print (density)
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

						# actor.body.velocity[0] += dragForceX
						# actor.body.velocity[1] += dragForceY

						rotatedForce = Vec2d(dragForceX, dragForceY)
						rotatedForce = rotatedForce.rotated(-actor.body.angle)
						actor.body.apply_impulse_at_local_point(rotatedForce, [0,0])
						# actor.body.apply_impulse_at_local_point([dragForceX,dragForceY], [0,0])
						
			# when you enter a new sphere of influence, regenerate the orbit information
			if strongestAttractor is not actor.orbiting or actor.orbiting is None:
				actor.leaveFreefall(0)
				actor.orbiting = strongestAttractor
				
			# figure out if the actor is freefalling by seeing if any engines or collisions have moved it.
			if actor.doModuleEffects(actor.keyStates, self.timestepSize):
				actor.leaveFreefall(0)
				for module in actor.modules:
					if module.enabled:
						if module.active:
							for giveResource, giveQuantity in list(module.resources.items()): #module.produces.items():
								if giveResource == 'thrust':
									if actor.keyStates['up']:
										position = actor.body.position + module.offset + module.illuminatorOffset
										position = rotate_point(position, actor.body.angle, actor.body.position)
										self.illuminators.append(Illuminator(position))

			# if it is freefalling, move it along the orbital track.
			if actor.freefalling and actor.orbit is not None:
				actor.orbit.updTime(self.timestepSize)
				cartesian = actor.orbit.cartesianCoordinates(actor.orbit.tAn)
				actor.body.position =  [cartesian[0] + actor.orbiting.body.position[0] ,cartesian[1] + actor.orbiting.body.position[1]]

				# you also must update the actor's velocity, or else when it leaves the track it will have the same velocity it entered with, leading to weird jumps.
				trackSpeed = actor.orbit.getSpeed(actor.orbit.tAn)

				# there is almost definitely a way to figure this out from the ellipse's properties. You would need to find tangent to the ellipse. But I figured it out by computing the position one step into the future, and then finding the angle between positions.
				futureSteptAn = actor.orbit.tAnAtTime(self.timestepSize)
				futureStepCoordinates = actor.orbit.cartesianCoordinates(futureSteptAn)
				adjustedFutureStep = [futureStepCoordinates[0] + actor.orbiting.body.position[0] , futureStepCoordinates[1] + actor.orbiting.body.position[1]]
				angle = math.atan2( adjustedFutureStep[1] - actor.body.position[1], adjustedFutureStep[0] - actor.body.position[0] )

				actor.body.velocity = [trackSpeed * math.cos(angle), trackSpeed * math.sin(angle)]

			else:
				# if it is not freefalling, add some gravity to it so that it moves naturally, and try to recalculate the orbit.
				if not actor.exemptFromGravity:
					self.gravitate(actor, strongestForce)
				if actor.stepsToFreefall > 0:
					actor.stepsToFreefall -= 1
				else:
					actor.freefalling = True
					actor.exemptFromGravity = False
					try:
						actor.orbit = Orbit.fromStateVector(numpy.array([actor.body.position[0] - actor.orbiting.body.position[0] ,actor.body.position[1] - actor.orbiting.body.position[1],1]), numpy.array([actor.body.velocity[0] ,actor.body.velocity[1] ,1]), actor.orbiting.APBody, Time('2000-01-01 00:00:00'), actor.name + " orbit around " + actor.orbiting.planetName)
					except:
						actor.orbit = None

	def rotatePolygon(self, points, angle):
		return Rotate2D(points,(0,0),angle)

	def transformForBuild(self, position):
		# map a position in the game world, where 0,0 is in a corner and numbers are very large, onto pixels on the screen with 0,0 in the middle. Handles zooming and offsetting. This one is for what you see in the build menu.
		transformedPosition = [0,0] #* self.zoom  # shrink or expand everything around the 0,0 point
		transformedPosition[0] = -position[0] * self.zoom
		transformedPosition[1] = -position[1] * self.zoom
		transformedPosition[0] += 0.5 * self.resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
		transformedPosition[1] += 0.5 * self.resolution[1] # add half the height.
		return transformedPosition

	def antiTransformForBuild(self, position):
		# performs the inverse operation to transformForBuild, used to map the mouse cursor to coordinates in the game world.
		transformedPosition = [0,0] #* self.zoom  # shrink or expand everything around the 0,0 point
		transformedPosition[0] = (-position[0]) + 0.5 * self.resolution[0] # add half the width of the screen, to get to the middle. 0,0 is naturally on the corner.
		transformedPosition[1] = (-position[1]) + 0.5 * self.resolution[1] # add half the height.
		transformedPosition[0] = transformedPosition[0] / self.zoom
		transformedPosition[1] = transformedPosition[1] / self.zoom
		# transformedPosition[0] = -transformedPosition[0] + 0.5 * self.resolution[0]
		# transformedPosition[1] = -transformedPosition[1] + 0.5 * self.resolution[1]
		# transformedPosition = rotate_point(transformedPosition, math.pi, [0.5 * self.resolution[0], 0.5 * self.resolution])




		
		return transformedPosition



	def drawCircle(self,color, position, radius):
		# pygame.draw.circle(self.screen, color, [int(position[0]), int(position[1])], int((radius * self.zoom)))
		pass

	def drawModuleForList(self, main_batch, module, iconSize, position):
		transformedPoints = []
		for point in module.points:
			transformedPoint = [0,0]
			transformedPoint[0] = ((point[0] * iconSize) + position[0])
			transformedPoint[1] = -((point[1] * iconSize) + position[1]) + self.resolution[1]
			# transformedPoints.append(transformedPoint) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		
		# print(transformedPoints)
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, self.resolution, module.color, module.outlineColor)

		# try:
		# 	# return pygame.draw.lines(self.screen, module.color, True, transformedPoints) # return the bounding rectangle
		# 	pass
		# except:
		# 	print('drawModuleForList error')

	def drawModuleForDrag(self,main_batch, module, position):
		# rotatedPoints = rotate_polygon(module.points, module.angle)  # orient the polygon according to the body's current direction in space.
		# # rotatedPoints = rotate_polygon(rotatedPoints,module.angle, -module.offset)  # orient the polygon according to the body's current direction in space.
		# transformedPoints = []
		# for rotatedPoint in rotatedPoints:
		# 	rotatedPoint[0] = (rotatedPoint[0] * self.zoom ) + position[0]
		# 	rotatedPoint[1] = (rotatedPoint[1] * self.zoom ) + position[1]
		# 	transformedPoints.append(rotatedPoint) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
		# try:
		# 	# pygame.draw.lines(self.screen, module.color, True, transformedPoints)
		# 	pass
		# except:
		# 	print('drawModuleForBuild error')
		# 	# print(transformedPoints)


		rotatedPoints = rotate_polygon(module.points, module.angle)
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			rotatedPoint[0] = (rotatedPoint[0] * self.zoom ) + position[0]
			rotatedPoint[1] = (rotatedPoint[1] * self.zoom ) + position[1]
			transformedPoint = rotatedPoint # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		
		# print(transformedPoints)
		renderAConvexPolygon(main_batch, transformedPoints,self.viewpointObject.body.position, self.zoom, self.resolution,  module.color, module.outlineColor)


	def drawModuleForBuild(self, main_batch, module):
		# rotatedPoints = rotate_polygon(module.points, module.angle)  # orient the polygon according to the body's current direction in space.
		# rotatedPoints = rotate_polygon(rotatedPoints,module.angle, -module.offset)  # orient the polygon according to the body's current direction in space.
		# transformedPoints = []
		# for rotatedPoint in rotatedPoints:
		# 	rotatedPoint[0] += module.offset[0]
		# 	rotatedPoint[1] += module.offset[1]
		# 	transformedPoints.append(self.transformForBuild(rotatedPoint)) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
		


		rotatedPoints = rotate_polygon(module.points, module.angle)
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			rotatedPoint[0] += module.offset[0]
			rotatedPoint[1] += module.offset[1]
			transformedPoint = self.transformForBuild(rotatedPoint) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		
		# print(transformedPoints)
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, self.resolution, module.color, module.outlineColor)


	def drawModuleEffects(self, module, actor):
		pass
		# put a circle in the middle if it is enabled, and a smaller red circle in the middle of that, if it is activated.
		# if module.enabled:
		# 	activeCircle = self.transformForView(module.offset + actor.body.position)
		# 	activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
		# 	self.drawCircle(module.color, activeCircle, 2)
		# 	if module.active:
		# 		self.drawCircle((255,0,0), activeCircle, 1)

		# # rocket engines have a line coming out of them.
		# if module.enabled and module.active:
		# 	for giveResource, giveQuantity in list(module.resources.items()): 
		# 		if giveResource == 'thrust':
		# 			if actor.keyStates['up']:
		# 				forceAngle = actor.body.angle + module.angle
		# 				force = [(giveQuantity * math.cos(addRadians(forceAngle, math.pi * 0.5))), giveQuantity * math.sin(addRadians(forceAngle, math.pi * 0.5) )]
		# 				activeCircle = self.transformForView(module.offset + actor.body.position)
		# 				activeCircle = rotate_point(activeCircle, actor.body.angle, self.transformForView(actor.body.position))
		# 				ananas = (int(activeCircle[0] + force[0] * self.zoom), int(activeCircle[1]+force[1] * self.zoom ) )
						# print ananas
						# pygame.draw.lines(self.screen, (255,255,200), True, [activeCircle,ananas])

	def drawModule(self, actor, module, main_batch): # draw the outline of the module.
		rotatedPoints = module.points
		rotatedPoints = rotate_polygon(rotatedPoints, actor.body.angle, [-module.offset[0], -module.offset[1]])
		rotatedPoints = rotate_polygon(rotatedPoints,module.angle)  # orient the polygon according to the body's current direction in space.
		transformedPoints = []

		for index, rotatedPoint in enumerate(rotatedPoints):
			# transformForView(self, position, viewpointObjectPosition, zoom, resolution):
			transformedPoint = transformForView(rotatedPoint + actor.body.position + module.offset,self.viewpointObject.body.position, self.zoom, self.resolution ) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
			transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
		
		renderAConvexPolygon(main_batch, transformedPoints, self.viewpointObject.body.position, self.zoom, self.resolution,  module.color, module.outlineColor, self.illuminators)

		self.drawModuleEffects(module, actor)

	def drawScreenFill(self, main_batch):

		fillTriangles = [0,0, 0,0, 0,resolution[1], resolution[0],resolution[1], resolution[0],0, 0,0, 0,0 ]

		main_batch.add(7, pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',fillTriangles), ('c4B',[255,255,255,255]*7))

	def drawActor(self, actor, main_batch):
		if actor.__class__ is Actor:
			for module in actor.modules:
					self.drawModule(actor, module, main_batch)

		if actor.__class__ is Atmosphere:

			transformedInnerPoints = []
			transformedOuterPoints = []

			#transformForView(self, position, viewpointObjectPosition, zoom, resolution):
			for index, point in enumerate(actor.innerPoints): 
				transformedInnerPoints.append(transformForView(point,self.viewpointObject.body.position, self.zoom, self.resolution))
			for index, point in enumerate(actor.outerPoints): 
				transformedOuterPoints.append(transformForView(point,self.viewpointObject.body.position, self.zoom, self.resolution))

			rendering = unwrapAtmosphereForGradient(actor.n_points, transformedInnerPoints, transformedOuterPoints, actor.color, actor.outerColor)

			main_batch.add(rendering[0], pyglet.gl.GL_TRIANGLE_STRIP, None, ('v2i',rendering[1]), ('c4B',rendering[2]))
		
		if actor.__class__ is Attractor:
			transformedPoints = []
			for point in actor.points:
				transformedPoint = transformForView(point + actor.body.position,self.viewpointObject.body.position, self.zoom, self.resolution) # transformForView does operations like zooming and mapping 0 to the center of the screen. 
				transformedPoints.append([int(transformedPoint[0]), int(transformedPoint[1])])
			
			renderAConvexPolygon(main_batch, transformedPoints,self.viewpointObject.body.position, self.zoom, self.resolution,  actor.color, actor.outlineColor)	

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
			# print('attractorCollision')
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

	def drawAPOrbit(self, main_batch, orbit, attractor, color):
		points = []
		n_points = 100
		sliceSize =  (2 * math.pi / n_points)
		for i in range(0,n_points):
			temp_vec3d = orbit.cartesianCoordinates(i *sliceSize)
			point = (temp_vec3d[0] + attractor.body.position[0], temp_vec3d[1] + attractor.body.position[1])
			point = transformForView(point,self.viewpointObject.body.position, self.zoom, self.resolution)
			points.append(point)

		transformedPoints = transformPolygonForLines(points)

		main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[100,100,100,255]*(transformedPoints[0])))

	def drawModuleListItem(self, main_batch, listItem, index):
		# draw one of the modules in the list in the build menu.
		buildListSpacing = 30

		itemSize = 2 * mag(numpy.array(getFarthestPointInPolygon(listItem.module.points)))

		iconSize = buildListSpacing / itemSize

		# listItem.boundingRectangle = [[buildListSpacing, index * buildListSpacing], [buildListSpacing, index+1 * buildListSpacing], [2 * buildListSpacing, index+1 * buildListSpacing], [2 * buildListSpacing, index * buildListSpacing]]

		gnarlypoints = []
		for point in listItem.module.points:
			gnarlypoints.append([point[0] * iconSize + buildListSpacing,point[1] * iconSize + (index * buildListSpacing ) ])


		listItem.boundingRectangle = boundPolygon(gnarlypoints)

		#self, main_batch, module, iconSize, position
		self.drawModuleForList(main_batch, listItem.module, iconSize, [buildListSpacing, index * buildListSpacing] )


	def drawHUDListItem(self,string, quantity, index):
		HUDlistItemSpacing = 15
		listXPosition = 30

		if quantity is None:
			# textsurface = self.font.render(string, False, (255, 255, 255))
			pass
		else:
			# textsurface = self.font.render(string + str(quantity), False, (255, 255, 255))
			pass

		# self.screen.blit(textsurface,(listXPosition,index * HUDlistItemSpacing))
		return index + 1

	def createHUDNavCircle(self):
		# this function predraws the nav circle. because it is just static lines, it does not need to be recalculated every frame.
		for n in range(0,self.n_navcircle_lines):
			angle = n * (2 * math.pi / self.n_navcircle_lines)
			start = ((self.navcircleInnerRadius * math.cos(angle)) + (self.resolution[0]*0.5) , (self.navcircleInnerRadius* math.sin(angle)) +( self.resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius + self.navcircleLinesLength) * math.cos(angle)+ (self.resolution[0]*0.5), (self.navcircleInnerRadius + self.navcircleLinesLength) * math.sin(angle)+ (self.resolution[1]*0.5))
			self.navCircleLines.append([start, end])
			# pygame.draw.lines(self.screen, (100,100,100), True, (start,end))


		

	def drawHUD(self, main_batch):
		if self.player is None:
			return

		# show the player what resources are available
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
			i = self.drawHUDListItem(str(availableResource) + ': ', availableQuantity, i)
		i = self.drawHUDListItem('', None, i) # blank line as a separator

		i = self.drawHUDListItem('freefalling: ', self.viewpointObject.freefalling, i)
		i = self.drawHUDListItem('landed: ', self.viewpointObject.exemptFromGravity, i)
		if self.player.orbiting is not None:
			i = self.drawHUDListItem('attractor: ', self.viewpointObject.orbiting.planetName, i)
			if self.player.orbit is not None:
				i = self.drawHUDListItem('orbit valid', None, i)
			else:
				i = self.drawHUDListItem('no orbit', None, i)
		i = self.drawHUDListItem('', None, i) # blank line as a separator

		i = self.drawHUDListItem('player: ', self.viewpointObject.isPlayer, i)
		i = self.drawHUDListItem('warp: ', self.timestepSize * 3 * 100, i)
		i = self.drawHUDListItem('zoom: ', self.zoom, i)
		i = self.drawHUDListItem('paused: ', self.paused, i)
		
		# print the navcircle
		for line in self.navCircleLines:

			transformedPoints = transformPolygonForLines(line)
			main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[50,50,50,255]*(transformedPoints[0])))


		# blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
		# angle = self.viewpointObject.desiredAngle
		# start = ((blipLength * math.cos(angle)) + (self.resolution[0]*0.5) , (blipLength* math.sin(angle)) +( self.resolution[1] * 0.5) )
		# end = ((self.navcircleInnerRadius) * math.cos(angle)+ (self.resolution[0]*0.5), (self.navcircleInnerRadius) * math.sin(angle)+ (self.resolution[1]*0.5))
		# # pygame.draw.lines(self.screen, (200,0,10), True, (start,end))
		# transformedPoints = transformPolygonForLines([start,end])
		# main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[100,0,0,255]*(transformedPoints[0])))

		
		# # draw the actor's orbits
		# for actor in self.actors:
		# 	if actor.orbit is not None:
		# 		self.drawAPOrbit(main_batch, actor.orbit, actor.orbiting, (100,100,100))


		

	def shootABullet(self, gunModule, actor):
		if gunModule.moduleType is 'cannon 10':
			#(self, name, modulesList, position, velocity, angle, isPlayer=False):
			bulletPosition = [actor.body.position[0] + gunModule.offset[0] + gunModule.barrelHole[0], actor.body.position[1] + gunModule.offset[1] + gunModule.barrelHole[1]]
			bulletPosition = rotate_point(bulletPosition, gunModule.angle, actor.body.position + gunModule.offset )
			bulletPosition = rotate_point(bulletPosition, actor.body.angle, actor.body.position )
			bulletVelocity = [gunModule.muzzleVelocity * math.cos(gunModule.angle + actor.body.angle - 0.5*math.pi) , gunModule.muzzleVelocity * math.sin(gunModule.angle + actor.body.angle- 0.5*math.pi)]
			bullet = Actor('cannonshell 10', [Module('cannonshell 10', (0,0))], bulletPosition ,bulletVelocity,0 ,False)
			self.add(bullet)
			self.illuminators.append(Illuminator(bulletPosition))


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
		self.add(playersNewShip)

	def dropModuleIntoBuildArea(self, module, position):
		module.offset = self.antiTransformForBuild(position)
		print(self.antiTransformForBuild(position))
		self.modulesInUse.append(module)

	def buildMenuGraphics(self, main_batch):
		window.clear()

		main_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		# self.screen.fill((200,200,200))
		self.drawScreenFill(main_batch)
		# main_batch.draw()

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

	def graphics(self):

		window.clear()

		main_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(2)

		for attractor in self.attractors:
			if attractor.atmosphere != None:
				self.drawActor(attractor.atmosphere, main_batch)
			self.drawActor(attractor, main_batch)
		for actor in self.actors:
			self.drawActor(actor, main_batch)


		if self.showHUD:
			blipLength = (self.navcircleInnerRadius-self.navcircleLinesLength)
			angle = self.viewpointObject.body.angle - 0.5 * math.pi
			start = ((blipLength * math.cos(angle)) + (self.resolution[0]*0.5) , -(blipLength* math.sin(angle)) +( self.resolution[1] * 0.5) )
			end = ((self.navcircleInnerRadius) * math.cos(angle)+ (self.resolution[0]*0.5),- (self.navcircleInnerRadius) * math.sin(angle)+ (self.resolution[1]*0.5))
			# pygame.draw.lines(self.screen, (200,0,10), True, (start,end))
			transformedPoints = transformPolygonForLines([start,end])
			main_batch.add(transformedPoints[0], pyglet.gl.GL_LINES, None, ('v2i', transformedPoints[1]), ('c4B',[200,0,0,255]*(transformedPoints[0])))

		
		main_batch.draw()
		second_batch = pyglet.graphics.Batch()
		pyglet.gl.glLineWidth(1)

		if self.showHUD:
			self.drawHUD(second_batch)
			# draw the actor's orbits
			for actor in self.actors:
				if actor.orbit is not None:
					self.drawAPOrbit(second_batch, actor.orbit, actor.orbiting, (100,100,100))


		second_batch.draw()


			

	def step(self):


		if not self.buildMenu:
			self.player = self._getPlayer()
			if not self.paused:
				self.physics()
			self.graphics()

			self.illuminators = []
		else:
			self.buildMenuGraphics()


		

	def setup(self):

		self.createHUDNavCircle()

		planet_erf = Attractor('earth', [1,1], self.gravitationalConstant)
		planet_moon = Attractor('moon',[1000000,-1000000], self.gravitationalConstant)
		self.add(planet_erf)
		self.add(planet_moon)
		dinghy_instance = Actor('NPC dinghy', dinghy,(1000000, -1080100), [20000,0], 0)
		lothar_instance = Actor('NPC lothar', lothar,(10000, -345050), [45000,0], 0.6 * math.pi)
		lothar_instance2 = Actor('player Lothar', lothar,(100, -320030), [0,0], 0, True)
		boldang_instance = Actor('NPC boldang', boldang,(-100, -320050), [0,0],0)
		bigmolly_instance = Actor('NPC molly', bigmolly,(100, -350050), [45000,0],0)
		self.add(dinghy_instance)
		self.add(lothar_instance)
		self.add(lothar_instance2)
		self.add(bigmolly_instance)
		for module in lothar:
			self.availableModules.append(copy.deepcopy(module))

	def start(self):
		self.setup()		

Nirn = World()

@window.event
def on_key_press(symbol, modifiers):

	if symbol == key.ESCAPE:
		exit()

	elif symbol == key.LEFT:
		Nirn.player.keyStates['left'] = True
	elif symbol == key.RIGHT:
		Nirn.player.keyStates['right'] = True
	elif symbol == key.UP:
		Nirn.player.keyStates['up'] = True
	elif symbol == key.P:
		Nirn.paused = not Nirn.paused
	elif symbol == key.EQUAL:
		Nirn.zoom += Nirn.zoom * 0.5
	elif symbol == key.MINUS:
		Nirn.zoom -= Nirn.zoom * 0.5
	elif symbol == key.H:
		Nirn.showHUD = not Nirn.showHUD
	elif symbol == key.COMMA:
		Nirn.timestepSize += Nirn.timestepSize * 0.5
	elif symbol == key.PERIOD:
		Nirn.timestepSize -= Nirn.timestepSize * 0.5
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
		for module in Nirn.player.modules:
			if module.moduleType == 'cannon 10':
				Nirn.shootABullet(module, Nirn.player)

	# elif event.type == pygame.MOUSEBUTTONUP:
	# 	if self.buildMenu:
	# 		if event.button == 1:
	# 			if self.buildDraggingModule is not None:
	# 				self.dropModuleIntoBuildArea(self.buildDraggingModule, pygame.mouse.get_pos())
	# 				self.buildDraggingModule = None
				

@window.event
def on_key_release(symbol, modifiers):
    if symbol == key.LEFT:
    	Nirn.player.keyStates['left'] = False
    elif symbol == key.RIGHT:
    	Nirn.player.keyStates['right'] = False
    elif symbol == key.UP:
    	Nirn.player.keyStates['up'] = False

  # def inputs(self):
		# pass
	# 	for event in pygame.event.get():
	# 		if event.type == KEYDOWN and event.key == K_ESCAPE:
	# 			self.running = False
	# 		elif event.type == KEYDOWN and event.key == K_RIGHTBRACKET:
	# 			self.viewpointObjectIndex += 1
	# 			if self.viewpointObjectIndex >= len(self.actors):
	# 				self.viewpointObjectIndex = 0
	# 			self.viewpointObject = self.actors[self.viewpointObjectIndex]
	# 		elif event.type == KEYDOWN and event.key == K_LEFTBRACKET:
	# 			self.viewpointObjectIndex -= 1
	# 			if self.viewpointObjectIndex < 0:
	# 				self.viewpointObjectIndex = len(self.actors) - 1
	# 			self.viewpointObject = self.actors[self.viewpointObjectIndex]
	# 		elif event.type == KEYDOWN and event.key == K_EQUALS:
	# 			self.zoom += self.zoom * 0.5
	# 		elif event.type == KEYDOWN and event.key == K_MINUS:
	# 			self.zoom -= self.zoom * 0.5
	# 		elif event.type == KEYDOWN and event.key == K_COMMA:
	# 			self.timestepSize += self.timestepSize * 0.5
	# 		elif event.type == KEYDOWN and event.key == K_PERIOD:
	# 			self.timestepSize -= self.timestepSize * 0.5
	# 		elif event.type == KEYDOWN and event.key == K_LEFT:
	# 			self.player.keyStates['left'] = True
	# 		elif event.type == KEYUP and event.key == K_LEFT:
	# 			self.player.keyStates['left'] = False
	# 		elif event.type == KEYDOWN and event.key == K_RIGHT:
	# 			self.player.keyStates['right'] = True
	# 		elif event.type == KEYUP and event.key == K_RIGHT:
	# 			self.player.keyStates['right'] = False
	# 		elif event.type == KEYDOWN and event.key == K_UP:
	# 			self.player.keyStates['up'] = True
	# 		elif event.type == KEYUP and event.key == K_UP:
	# 			self.player.keyStates['up'] = False
	# 		elif event.type == KEYDOWN and event.key == K_DOWN:
	# 			self.player.keyStates['down'] = True
	# 		elif event.type == KEYUP and event.key == K_DOWN:
	# 			self.player.keyStates['down'] = False
	# 		elif event.type == KEYDOWN and event.key == K_h:
	# 			self.showHUD = not self.showHUD
	# 		elif event.type == KEYDOWN and event.key == K_p:
	# 			self.paused = not self.paused
	# 		elif event.type == KEYDOWN and event.key == K_b:
	# 			if self.buildMenu:
	# 				self.buildMenu = False
	# 			else:
	# 				self.buildMenu = True
	# 				self.paused = True

	# 				self.loadShipIntoBuildMenu(self.player)
	# 		elif event.type == KEYDOWN and event.key == K_y:
	# 			if self.buildMenu:
	# 				self.flyShipFromBuildMenu()
	# 		elif event.type == pygame.MOUSEBUTTONDOWN:
	# 			# event.button can equal several integer values:# 1 - left click# 2 - middle click# 3 - right click# 4 - scroll up# 5 - scroll down
	# 			if self.buildMenu:
	# 				if event.button == 1:
	# 					self.buildDraggingModule = self.getModuleFromCursorPosition(pygame.mouse.get_pos())
	# 		elif event.type == pygame.MOUSEBUTTONUP:
	# 			if self.buildMenu:
	# 				if event.button == 1:
	# 					if self.buildDraggingModule is not None:
	# 						self.dropModuleIntoBuildArea(self.buildDraggingModule, pygame.mouse.get_pos())
	# 						self.buildDraggingModule = None
						


def stepWithBatch(dt):
	pass


@window.event()
def on_mouse_press(x, y, button, modifiers):
    # pass
	# print('on_mouse_release')
	Nirn.mouseCursorPosition =[x,y]# Nirn.antiTransformForBuild([x,y])
    	# elif event.type == pygame.MOUSEBUTTONDOWN:
	# 	# event.button can equal several integer values:# 1 - left click# 2 - middle click# 3 - right click# 4 - scroll up# 5 - scroll down
	if Nirn.buildMenu:
		# if event.button == 1:
		if mouse.LEFT is button:
			# print('LEFT')
			Nirn.buildDraggingModule = Nirn.getModuleFromCursorPosition(Nirn.mouseCursorPosition)

@window.event()
def on_mouse_release(x, y, button, modifiers):
	Nirn.mouseCursorPosition = [x,y]#Nirn.antiTransformForBuild([x,y])
	# print('on_mouse_release')
	if Nirn.buildMenu:
		if mouse.LEFT is button:
			# print('drop')
			if Nirn.buildDraggingModule is not None:
				Nirn.dropModuleIntoBuildArea(Nirn.buildDraggingModule, Nirn.mouseCursorPosition)
				Nirn.buildDraggingModule = None

@window.event()
def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
    Nirn.mouseCursorPosition =[x,y]# Nirn.antiTransformForBuild([x,y])

@window.event()
def on_draw():

	Nirn.step()

	# main_batch.draw()


def hello():

	Nirn.start()


	# pyglet.gl.glLineWidth(2)

	pyglet.clock.schedule_interval(stepWithBatch, 0.01)
	pyglet.app.run()

cProfile.run('hello()')