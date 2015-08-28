# Thomas Macrina
# 150710
# Python implementation to take correspondences and produce color plots

import os
import csv
import math
from PIL import Image, ImageDraw
import numpy as np

# Pulled from Variable Scope
# http://variable-scope.com/posts/hexagon-tilings-with-python
class HexagonGenerator(object):
	"""Returns a hexagon generator for hexagons of the specified size."""
	def __init__(self, edge_length, row_factor):
		# self.short_radius = edge_length * math.sqrt(3) / 2
		self.row_factor = row_factor
		self.edge_length = edge_length

	@property
	def col_width(self):
		# return self.short_radius
		return self.edge_length

	@property
	def row_height(self):
		return self.edge_length * self.row_factor

	def __call__(self, row, col):
		x = col * self.col_width
		y = row * self.row_height
		for angle in range(0, 360, 60):
			x += math.sin(math.radians(angle)) * self.col_width
			y += math.cos(math.radians(angle)) * self.row_height
			yield x
			yield y

# Adapted from mpicbg.ij.util.Util.java
# https://github.com/axtimwalde/mpicbg/blob/master/mpicbg/src/main/java/mpicbg/ij/util/Util.java
def color_vector(d):
	dx = d[0]
	dy = d[1]
	r, g, b = 0, 0, 0
	a = math.sqrt(dx*dx + dy*dy)
	if a != 0.0:
		o = (math.atan2(dx/a, dy/a) + math.pi) / math.pi * 3
		if o < 3:
			r = min(1.0, max(0, 2.0 - o)) * a
		else:
			r = min(1.0, max(0, o - 4.0)) * a
		
		o += 2
		if o >= 6:
			o -= 6
		if o < 3:
			g = min(1.0, max(0, 2.0 - o)) * a
		else:
			g = min(1.0, max(0, o - 4.0)) * a

		o += 2
		if o >= 6:
			o -= 6
		if o < 3:
			b = min(1.0, max(0, 2.0 - o)) * a
		else:
			b = min(1.0, max(0, o - 4.0)) * a

	return int(r*255), int(g*255), int(b*255)

def create_color_plot_from_meshes(initial_nodes, final_nodes):
	x_min = np.min(initial_nodes[:,0])
	y_min = np.min(initial_nodes[:,1])

	x_unq = np.unique(initial_nodes[:,0])
	y_unq = np.unique(initial_nodes[:,1])

	x_dist = x_unq[1] - x_unq[0]
	y_dist = y_unq[1] - y_unq[0]

	# x_dist = node_dist / 2.0
	# y_dist = node_dist * math.sqrt(3) / 2

	# This will cause gaps where rounding doesn't perfectly match integers
	# Could try snapping points to predetermined grid
	# Good enough for now, though
	x_idx = np.round((initial_nodes[:,0] - x_min)/x_dist).astype(int)
	y_idx = np.round((initial_nodes[:,1] - y_min)/y_dist).astype(int)

	d = final_nodes - initial_nodes
	max_dist = np.max(np.linalg.norm(d, axis=1))
	d /= max_dist

	colors = np.apply_along_axis(color_vector, 1, d)

	hex_edge = 10
	hexagon_generator = HexagonGenerator(hex_edge, y_dist/x_dist)
	width = int((np.max(x_idx) + 1.5) * hexagon_generator.col_width)
	height = int((np.max(y_idx) + 1.5) * hexagon_generator.row_height)
	size = width, height+18
	img = Image.new("RGB", size, "white")
	draw = ImageDraw.Draw(img)
	for (i, j, color) in zip(x_idx, y_idx, colors):
		hexagon = hexagon_generator(j, i)
		draw.polygon(list(hexagon), outline=tuple(color), fill=tuple(color))
	draw.text((1, height+1), "max dist: " + str(max_dist), (0,0,0))
	return img

def write_image_from_points(initial_nodes, final_nodes, filename):
	im = create_color_plot_from_meshes(initial_nodes, final_nodes)
	im.save(filename)

def demo():
	a = [i for i in range(1,10)*9]
	b = np.sort(a)
	initial_nodes = np.array([a, b]).T
	final_nodes = np.array([a, b]).T
	final_nodes[3,1] = 5
	final_nodes[4,1] = 6
	final_nodes[5,1] = 7
	write_image_from_points(initial_nodes, final_nodes, 'test.png')