import numpy as np
import pandas as pd
import points

class MatchInspect():
	"""Object to control display of correspondence point pairs with filters so
	they can be manually inspected and editted.

	df: pandas dataframe with pre & post points along with current filters
	points_controller: points.Controller object for state server interaction
	"""

	def __init__(self, df, points_controller):
		self.df = df
		self.pc = points_controller
		self.filter = False

	def display(self):
		sub_df = self.df[self.df['filter'] .== self.filter]
		pts = list(zip(map(int,sub_df.pre_x), 
                        map(int,sub_df.locs_2), 
                        map(int,sub_df.locs_3))), list(map(int,sub_df.index))
		self.pc.set()

	def sync(self):
