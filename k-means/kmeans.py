'''
This allows to handle something
'''
# import pickle
import time
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from scipy.ndimage.interpolation import shift
import helita.io.ncdf as nd

class KMeans():
    '''
    This is the class
    '''

    def __init__(self, verbose=True):
        '''
        initializes variables tests
        '''

        if verbose:
            print("True")

    #def read_lim_interp()
    #    self.read_data_npz....
    #    self.wavelength_distinction...
    #    self.individual_spectral_data...
    #    self.interp....

    def read_data_npz(self, npzfilename="/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz"):
        '''
        reads data, loads data, saves the loading of the data
        '''

        self.data = np.load(npzfilename)
        self.i_3 = self.data["arr_0"]
        self.wvl = self.data["arr_1"]
        print("True")

    def read_data_ncdf(self, rbfilename='output_ray_l2d90x40r.ncdf'):
        self.wvl1 = nd.getvar(rbfilename,'wavelength')
        self.inte = nd.getvar(rbfilename, 'intensity', memmap=True)
        print("True")

    def read_data_pck(self, kmeansfilename='k-means.pck'):
        pick_in = open(kmeansfilename, 'rb')
        self.k_m = pickle.load(pick_in)

    def wavelength_distinction(self, delta=5):
        '''
        This allows to cut the wavelength range in different spectral lines.
        Find where the value of delwvl is larger than the previous index
        by delta, appends that value into limits array
        Loops through the limits array, and sees if the value for one index
        in limits has a greater difference than 1 than the previous index.
        '''

        delwvl = self.wvl1-shift(self.wvl1, 1)
        new_delwvl = shift(self.wvl1, -1) - self.wvl1
        self.limits = [v for v in range(1, len(delwvl)) if np.abs(delwvl[v]-new_delwvl[v]) > delta]
        self.limits.append(len(self.wvl1))

    def individual_spectral_data(self, delta=5):
        '''
        The profile at inte is added to the array new_inte
        '''

        self.new_inte = {}

        if not hasattr(self, 'limits'):
            self.wavelength_distinction(delta)

        count = 0
        for ind in range(0, np.size(self.limits)-1):
            if self.limits[ind+1]-self.limits[ind] > 1:
                self.new_inte[count] = self.inte[:, :, self.limits[ind]:self.limits[ind+1]-1] # var
                count = count + 1

    def interp(self):
        '''
        interpolation of the axis
        plots wvl against new_inte in an uniform axis(wvlax)
        '''

        delwvl = self.wvl1-shift(self.wvl1, 1)
        count = 0
        self.new_inte1 = {}
        self.wvlax = {}
        for ind in range(0, np.size(self.limits)-1):
            if self.limits[ind+1]-self.limits[ind] > 1:
                print(self.wvl[self.limits[ind]], self.limits[ind+1])
                mindelwvl = np.min(delwvl[self.limits[ind]:self.limits[ind+1]-1])
                n_points = np.min(((np.max(self.wvl)-np.min(self.wvl))/mindelwvl, 3000))
                
                if (n_points>5)
                    max_value = np.max(self.wvl[self.limits[ind]:self.limits[ind+1]-1])
                    min_value = np.min(self.wvl[self.limits[ind]:self.limits[ind+1]-1])
                    self.wvlax[count] = np.linspace(min_value, max_value, num=n_points)
                    print(n_points, mindelwvl, np.shape(self.wvlax[count]))
                    inte1 = np.zeros((self.inte.shape[0],
                                      self.inte.shape[1], np.shape(self.wvlax[count])[0]))
                    print(self.inte.shape[0], self.inte.shape[1], np.shape(inte1))
                    print('wvl', np.shape(self.wvlax[count]),
                          np.min(self.wvlax[count]), np.max(self.wvlax[count]),
                          np.min(self.wvl[self.limits[ind]:self.limits[ind+1]-1]),
                          np.max(self.wvl[self.limits[ind]:self.limits[ind+1]-1]))

                    for ind2 in range(0, len(self.new_inte[count][:, 0, 0])):
                        print('ind2=', ind2)
                        for ind3 in range(0, len(self.new_inte[count][0, :, 0])):
                            ind3[ind2, ind3, :] = np.interp(self.wvlax[count],
                                                            self.wvl[self.limits[ind]:
                                                                self.limits[ind3+1]-1],
                                                            self.new_inte[count][ind2, ind3, :])

                self.new_inte1[count] = inte1
                print('new_inte1', count, np.shape(self.new_inte1[count]))
                if count == 2:
                    mg2_cube = self.new_inte1[2]
                    np.save('new_inte1_02.npy', mg2_cube)

                count += 1

    def time_import(self):
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters
            computes the inertia of the MiniBatchKMeans
        inputs: t_m, inertia, t0, outputs:
        '''

        i3_2d = self.i_3.reshape((self.i_3.shape[0]*self.i_3.shape[1], self.i_3.shape[2]))

        self.t_m = np.zeros(30)
        self.inertia = np.zeros(30)
        self.t_zero = time.time()
        for i in range(0, 30):
            print(i)
            mini_km = MiniBatchKMeans(n_clusters=(i+1)*10, n_init=10).fit(i3_2d[:, :])
            t_mini_batch = time.time() - t_zero
            t_m[i] = t_mini_batch/((i+1)*10)
            inertia[i] = mini_km.inertia_
        print(t_m)
        plt.subplot(2, 1, 1)
        plt.plot(t_m)
        plt.subplot(2, 1, 2)
        plt.plot(inertia)
        plt.show()

    def mini_batch_fit(self):
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters
        inputs: t0
        '''

        t_zero = time.time()
        mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i_3[:, 1000:2000])
        t_mini_batch = time.time() - t_zero
        print("time = ", t_mini_batch)
        print("inertia = ", mini_km.inertia_)
        print("init = ", mini_km.n_init)


    def create_km_map(self):
        '''
        creates the image of the km_map_datacube
        shows the locations of the k-means clusters for each labels
        outputs: prints the image of the km_map datacube
        '''

        dim_i_3 = np.shape(self.i_3)

        plt.figure(figsize=(33, 33))
        km_map_datacube = np.zeros((dim_i_3[0], dim_i_3[1], 30))

        for i in range(0, 30):
            plt.subplot(5, 6, i+1)
            plt.xlabel('X [DNs]')
            plt.ylabel('Time [Seconds]')
            www = np.where(self.k_m.labels_ == i)[0]
            xxx = www/dim_i_3[1]
            yyy = www%dim_i_3[1]
            km_map_datacube[xxx.astype(int), yyy.astype(int), i] = i
            plt.imshow(km_map_datacube[:, :, i], extent=[0, 157, 0, 3000], aspect='auto')
            plt.show()


    def create_profiles(self):
        '''
        creates the k_means maps
        plots the spectral profile for the different k-means labels
        wavelength on the x axis and intensity on the y axis
        inputs: wvl, outputs: prints the image of the k-means labels
        '''

        wvl = self.data["arr_1"]
        plt.figure(figsize=(30, 30))
        for i in range(0, 30):
            plt.subplot(5, 6, i+1)
            plt.xlabel('Wavelength - 2796.2 [$\AA$]', fontsize=20)
            plt.ylabel('Intensity/I_{c} [$erg s^{-1} cm^{-2} Hz^{-1} ster^{-1}$]', fontsize=20)
            plt.plot(wvl*10.-2795.37, self.k_m.cluster_centers_[i, 0:2000]*1e3/3.454e-6)
        plt.show()

    def create_spectral_map(self):
        '''
        reshapes into 2D array
        interpolation
        adjusts wvl axis
        creates a spectral profile map
        showing the
        location of the clusters in all the labels as a whole
        inputs: i3, wvl, ax, outputs: prints the image of the spectral map
        '''
        wvl = self.data["arr_1"]
        wvl_new = wvl*10.-2795.37
        plt.subplots_adjust(bottom=0.2, top=.9, left=0.15)
        plt.imshow(self.i_3[0, :, :], aspect='auto',
                   extent=(np.min(wvl_new), np.max(wvl_new), 1570, 0))
        plt.title('Spectral Profile Map for Mg II k & h', fontsize=15)
        plt.xlabel('Wavelength - 2796.2 [$\AA$]', fontsize=15)
        plt.ylabel('Time (s)', fontsize=15)
        self.axi = plt.gca()
        self.axi.invert_yaxis()
        plt.savefig('spectral_map.eps')
        plt.show()
        plt.show()
