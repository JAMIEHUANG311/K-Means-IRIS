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
    	initializes variables
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

	def interp():
		'''
		interpolation of the axis 
		'''

		count = 0
		new_inte1 = {}
		wvlAx = {}
		for k in range(0, len(wvlAx[count][0,:,0])):
			wvlAx[j,k,:] = np.interp(wvlAx[count], wvl[limits[i]:limits[i+1]-1], new_inte[count][j, k, :])


	def wavelength_distinctions(wvl):
		'''
		separates between wavelength range
		'''

		km_map = np.zeros(dim_i3[0]*dim_i3[1]).reshape(dim_i3[0],dim_i3[1])

		for k in range(30):
   			w=np.where(km.labels_ == k)[0]
    		print(w.shape)
    		x = (w/dim_i3[1])
    		y = w%dim_i3[1]
    		print(x.shape, y.shape)
    		print(x.max(), y.max())
    		print(dim_i3, w.shape[0], km_map.shape)
    	km_map[x.astype(int),y.astype(int)] = k


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

	def fit(self, t0):
		'''
		uses the MiniBatchKMeans function to fit the i3_2D data into clusters
		inputs: t0
		'''

		t0 = time.time()
		mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i3[:,1000:2000])
		t_mini_batch = time.time() - t0


	def create_km_map(self):
		'''
		creates the image of the km_map_datacube  
   		shows the locations of the k-means clusters for each labels
   		outputs: prints the image of the km_map datacube
   		'''

		plt.figure(figsize = (5,5))
		axes_style = {'linewidth':2}
		plt.subplots_adjust(bottom=0.15, top=.9, left=0.15)
		km_map_all = np.sum(km_map_datacube[:,60:,:], axis=2)
		plt.imshow(km_map_all, extent=[0,157,0,3000], aspect='auto', cmap='rainbow')
		plt.xlabel('Time [Seconds]')
		plt.ylabel('Y [Pixels]')
		plt.colorbar()
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
