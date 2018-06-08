	import numpy as np
	import matplotlib.pyplot as plt
	from sklearn.cluster import KMeans
	from sklearn.cluster import MiniBatchKMeans
	from scipy.ndimage.interpolation import shift
	import pickle
	import time


	def __init__(self, data, i3, wvl, create_spectral_map,
				create_km_map, create_k-means_maps,
				time_import, verbose = true):
  		'''
		initializes variables tests
		''' 

		if verbose:
			print("True")
    		self.create_spectral_map = create_spectral_map
    		self.create_data_cube = create_data_cube
    		self.k-means_maps = k-means_maps
    		self.i3 = i3


	def read_data(self, i3, wvl, rbfilename='output_ray_l2d90x40r.ncdf', 
				npzfilename = "/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz"):
		'''
		reads data, loads data, saves the loading of the data
		'''

		self.data = np.load(npzfilename)
		self.i3 = data["arr_0"]
		self.wvl = data["arr_1"]
		pick_in = open(kmeansfilename, rbfilename) 
		self.km = pickle.load(pick_in)


	def wavelength_distinction(self, delt)
		'''
		Find where the value of delwvl is larger than the previous index by delta, appends that value into limits array
		Loops through the limits array, and sees if the value for one index in limits has a greater difference than 1 than the previous index, then the profile at inte is added to the array new_inte
		'''

		delta = delt
		delwvl = wvl-shift(wvl, 1)
		new_delwvl = shift(wvl, -1) - wvl
		limits = [v for v in range(1,len(delwvl)) if np.abs(delwvl[v]-new_delwvl[v]) > delta]
		limits.append(len(wvl))

		new_inte = {}
	
		count = 0
		for i in range(0, np.size(limits)-1):
			wavelengthmin = limits[i]
			wavelengthmax = limits[i+1]-1
			if(limits[i+1]-limits[i]>1):
				new_inte[count] = inte[:,:,wavelengthmin:wavelengthmax] # var
				count+=1

	def interp():
		'''
		interpolation of the axis 
		plots wvl against new_inte in an uniform axis(wvlAx)
		'''

		count = 0
		new_inte1 = {}
		wvlAx = {}
		for i in range(0, np.size(limits)-1):
			if(limits[i+1]-limits[i] >1):
				print(wvl[limits[i]],limits[i+1])
				mindelwvl = np.min(delwvl[limits[i]:limits[i+1]-1])
				n_points = np.min(((np.max(wvl)-np.min(wvl))/mindelwvl, 3000))
				max_value = np.max(wvl[limits[i]:limits[i+1]-1]) 
				min_value = np.min(wvl[limits[i]:limits[i+1]-1])
				wvlAx[count] = np.linspace(min_value, max_value, num = n_points)
				print(n_points, mindelwvl, np.shape(wvlAx[count]))
				inte1 = np.zeros((inte.shape[0], inte.shape[1], np.shape(wvlAx[count])[0]))
				print(inte.shape[0], inte.shape[1], np.shape(inte1))
				print('wvl',np.shape(wvlAx[count]),np.min(wvlAx[count]),np.max(wvlAx[count]), np.min(wvl[limits[i]:limits[i+1]-1]),np.max(wvl[limits[i]:limits[i+1]-1]))
				for j in range(0, len(new_inte[count][:,0,0])):
					print('j=',j)
					for k in range(0, len(new_inte[count][0,:,0])):
						inte1[j,k,:] = np.interp(wvlAx[count], wvl[limits[i]:limits[i+1]-1], new_inte[count][j, k, :])

				new_inte1[count] = inte1 
				print('new_inte1',count,np.shape(new_inte1[count]))
				if count == 2: 
					mgII_cube = new_inte1[2]
					np.save('new_inte1_02.npy', mgII_cube)

				count +=1

	def time_import(self, tm, inertia, t0):
		'''
		uses the MiniBatchKMeans function to fit the i3_2D data into clusters
   		computes the inertia of the MiniBatchKMeans
    	inputs: tm, inertia, t0, outputs: 
    	'''

		tm = np.zeros(30)
		inertia = np.zeros(30)
		t0 = time.time()
		for i in range(0,30):
    		print(i)
    		mini_km = MiniBatchKMeans(n_clusters = (i+1)*10, n_init = 10).fit(i3_2D[:,:])
    		t_mini_batch = time.time() - t0
    		tm[i] = t_mini_batch/((i+1)*10)
    		inertia[i] = mini_km.inertia_
		print(tm)
		plt.subplot(2,1,1)
		plt.plot(tm)
		plt.subplot(2,1,2)
		plt.plot(inertia)
		plt.show()

	def mini_batch_fit(self, t0):
		'''
		uses the MiniBatchKMeans function to fit the i3_2D data into clusters
		inputs: t0
		'''

		t0 = time.time()
		mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i3[:,1000:2000])
		t_mini_batch = time.time() - t0
		print("time = " t_mini_batch)
		print("inertia = " mini_km.inertia_)
		print("init = " mini_km.n_init)


	def create_km_map(self):
		'''
		creates the image of the km_map_datacube  
   		shows the locations of the k-means clusters for each labels
   		outputs: prints the image of the km_map datacube
   		'''

		plt.figure(figsize = (5,5))
		axes_style = {'linewidth':2}
		plt.subplots_adjust(bottom=0.15, top=.9, left=0.15)
		plt.imshow(km_map_datacube[:,60:,2], extent=[0,157,0,3000], aspect='auto')
		plt.xlabel('Time [Seconds]')
		plt.ylabel('Y [Pixels]')
		plt.savefig('map_2.eps')
		plt.show()

		plt.figure(figsize = (5,5))
		axes_style = {'linewidth':2}
		plt.subplots_adjust(bottom=0.15, top=.9, left=0.15)
		plt.imshow(km_map_datacube[:,60:,19], extent=[0,157,0,3000], aspect='auto')
		plt.xlabel('Time [Seconds]')
		plt.ylabel('Y [Pixels]')
		plt.savefig('map_19.eps')
		plt.show()

		plt.figure(figsize = (5,5))
		axes_style = {'linewidth':2}
		plt.subplots_adjust(bottom=0.15, top=.9, left=0.15)
		plt.imshow(km_map_datacube[:,60:,22], extent=[0,157,0,3000], aspect='auto')
		plt.xlabel('Time [Seconds]')
		plt.ylabel('Y [Pixels]')
		plt.savefig('map_22.eps')
		plt.show()



	def create_k_means_maps(self, wvl):
		'''
		creates the k_means maps 
		plots the spectral profile for the different k-means labels
		wavelength on the x axis and intensity on the y axis
		inputs: wvl, outputs: prints the image of the k-means labels
		'''

		plt.figure(figsize = (30,30))
		for i in range(0, 30):
        	plt.subplot(5,6,i+1)
        	plt.xlabel('Wavelength - 2796.2 [$\AA$]', fontsize=20)
        	plt.ylabel('Intensity/I_{c} [$erg s^{-1} cm^{-2} Hz^{-1} ster^{-1}$]', fontsize=20)
       		plt.plot(wvl*10.-2795.37, km.cluster_centers_[i,0:2000]*1e3/3.454e-6)
       	plt.show()

 	def create_spectral_map(self, i3, wvl, ax):
		'''
		reshapes into 2D array
   		interpolation
   		adjusts wvl axis
   		creates a spectral profile map showing the location of the clusters in all the labels as a whole
   		inputs: i3, wvl, ax, outputs: prints the image of the spectral map
   		'''

		wvl_new = wvl*10.-2795.37
		plt.subplots_adjust(bottom=0.2, top=.9, left=0.15)
		axes_style = {'linewidth':2}
		plt.imshow(i3[0,:,:], aspect = 'auto', extent=(np.min(wvl_new),np.max(wvl_new), 1570, 0))
		plt.title('Spectral Profile Map for Mg II k & h', fontsize = 15)
		plt.xlabel('Wavelength - 2796.2 [$\AA$]', fontsize = 15)
		axes_style = {'linewidth':2}
		plt.ylabel('Time (s)', fontsize = 15)
		ax = plt.gca()
		ax.invert_yaxis()
		plt.savefig('spectral_map.eps')
		plt.show()
		plt.show()
