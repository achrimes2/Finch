def resample(galaxy):
	# -*- coding: utf-8 -*-
	"""
	Created on Tue Dec 12 12:34:39 2017
	
	@author: Ashch
	"""
	# Code for resampling an image to vary it by the error on each pixel.

	import numpy as np

	Inum = 1
	Snum = 2
		
	image = galaxy[Inum].data
	sigma = galaxy[Snum].data

	image = np.random.normal(image,sigma)
    
	temp = galaxy
	temp[1].data = image  
	
	return temp


