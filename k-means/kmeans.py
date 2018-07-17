'''
This allows to handle something
JMS Explain things here.
'''
import pickle
# import time
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.cluster import KMeans
from sklearn.cluster import MiniBatchKMeans
from scipy.ndimage.interpolation import shift
import helita.io.ncdf as nd

class kmeans():
    '''
    JMS No need of this
    This is the class
    '''

    def __init__(self, verbose=True):
        '''
        JMS No need of this
        initializes variables tests
        '''

        if verbose:
            self.verbose = True

    def read_data_npz(self, npzfilename="new_inte1_02.npy.npz"):
        '''
        reads data, loads data, saves the loading of the data
        ???
        JMS Lets first check if we have something that saves
        these information in npz format (np.savez)
        '''

        self.data = np.load(npzfilename)
        self.i_3 = self.data["arr_0"]
        self.wvl = self.data["arr_1"]
        if self.verbose:
            print("warning")


    def read_data_ncdf(self, ncdffilename='output_ray_l2d90x40r.ncdf'):
        '''
        Reading from the original RH code ncdf files
        '''

        self.wvl_orig = nd.getvar(ncdffilename, 'wavelength')
        self.int_orig = nd.getvar(ncdffilename, 'intensity', memmap=True)


    def read_data_pck(self, filename='k-means.pck'):

        def save_all_pkl(self, filename='k-means.pck'):
            data = self.__dict__
            filehandler = open(filename+'.opy', 'wb')
            pickle.dump(data, filehandler, protocol=4)
        def load_all_pkl(self, filename='k-means.pck'):
            filehandler = open(filename, 'rb')
            data = pickle.load(filehandler)
            for idata in data.keys():
                setattr(self, idata, data[idata])

        # pick_in = open(kmeansfilename, 'rb')
        # self.k_m = pickle.load(pick_in)


    def spect_lines_limits(self, wvl_delta=5):
        '''
        This allows to cut the wavelength range in different spectral lines.
        Find where the value of delwvl_up is larger than the previous index
        by wvl_delta, appends that value into wvl_lmts array
        Loops through the wvl_lmts array, and sees if the value for one index
        in wvl_lmts has a greater difference than 1 than the previous index.
        JMS follow standards of documentation.
        http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
        '''

        delwvl_up = self.wvl_orig-shift(self.wvl_orig, 1)
        delwvl_dn = shift(self.wvl_orig, -1) - self.wvl_orig
        self.wvl_lmts = [v for v in range(1, len(delwvl_up)) 
                         if np.abs(delwvl_up[v]-delwvl_dn[v]) > wvl_delta]
        self.wvl_lmts.append(len(self.wvl_orig))

    def print_line_list(self, wvl_delta=5):
        '''
        The profile at inte is added to the array int_indv
        '''

        if not hasattr(self, 'wvl_lmts'):
            self.spect_lines_limits(wvl_delta)

        count = 0
        for ind in range(0, np.size(self.wvl_lmts)-1):
            if self.wvl_lmts[ind+1]-self.wvl_lmts[ind] > 1:
                print('%i,wvl=%0.1f' %(count, self.wvl_orig[self.wvl_lmts[ind]]))
                count = count + 1


    def individual_spectral_data(self, wvl_delta=5): # break_lines
        '''
        The profile at inte is added to the array int_indv
        '''
        if not hasattr(self, 'int_indv'): 
            self.int_indv = {}
        if not hasattr(self, 'wvl_indv'):
            self.wvl_indv = {}
            print('True')
        if not hasattr(self, 'wvl_lmts'):
            self.spect_lines_limits(wvl_delta)

        count = 0
        for ind in range(0, np.size(self.wvl_lmts)-1):
            if self.wvl_lmts[ind+1]-self.wvl_lmts[ind] > 1:
                self.int_indv[count] = self.int_orig[:, :, self.wvl_lmts[ind]:self.wvl_lmts[ind+1]-1] # var
                self.wvl_indv[count] = self.wvl_orig[self.wvl_lmts[ind]:self.wvl_lmts[ind+1]-1]
                count = count + 1

    def linear_spect(self, id_lines=0, savefile=None, nwvl_min=500):
        '''
        interpolation of the axis
        plots wvl against int_indv in an uniform axis(wvl_linear)
        '''

        self.int_linear = {}

        if not hasattr(self, 'wvl_linear'):
            self.wvl_linear = {}

        delwvl_up = np.gradient(self.wvl_indv[id_lines])

        if not hasattr(self, 'int_linear'):
            self.int_linear = {}

        mindelwvl = np.min(delwvl_up)
        n_points = np.min([(np.max(self.wvl_indv[id_lines])-np.min(self.wvl_indv[id_lines]))/mindelwvl, nwvl_min])

        if n_points > 5:
            max_value = np.max(self.wvl_indv[id_lines])
            min_value = np.min(self.wvl_indv[id_lines])
            self.wvl_linear[id_lines] = np.linspace(min_value, max_value, num=n_points)
            inte1 = np.zeros((self.int_orig.shape[0],
                              self.int_orig.shape[1], np.shape(self.wvl_linear[id_lines])[0]))

            for ind2 in range(0, len(self.int_indv[id_lines][:, 0, 0])):
                for ind3 in range(0, len(self.int_indv[id_lines][0, :, 0])):
                    inte1[ind2, ind3, :] = np.interp(self.wvl_linear[id_lines],
                                                     self.wvl_indv[id_lines],
                                                     self.int_indv[id_lines][ind2, ind3, :])

        self.int_linear[id_lines] = inte1

        if self.verbose:
            print('int_linear', id_lines, np.shape(self.int_linear[id_lines]))
        if savefile != None:
            mg2_cube = self.int_linear[id_lines]
            np.save(savefile + '.npy', mg2_cube)


    def read_lim_interp(self, ncdffilename='/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/output_ray_l2d90x40r.ncdf', ind=0, wvl_delta=5):
        self.read_data_ncdf(ncdffilename=ncdffilename)
        self.individual_spectral_data(wvl_delta=wvl_delta)
        self.linear_spect(ind=ind)



    def km_clusters(self, max_cluster_niter=5):
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters,
        computes the inertia of the MiniBatchKMeans
        outputs: graph of t_m and inertia
        '''

        i3_2d = self.i_3.reshape(self.i_3.shape[0]*self.i_3.shape[1], self.i_3.shape[2]) ## int_linear[line_id]
        self.inertia = np.zeros(max_cluster_niter) #km_inertia
        for i in range(0, max_cluster_niter):
            self.mini_km = MiniBatchKMeans(n_clusters=(i+1)*10, n_init=10).fit(i3_2d)
            self.inertia[i] = self.mini_km.inertia_
            ### save in an array (i+1)*10
            if self.verbose:
                print('description iteration=', i)
        plt.plot(self.t_m, self.inertia)
        plt.show()
        print(self.t_m)

    def mini_batch_fit(self, ind=0):
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters
        outputs: prints time and inertia
        '''

        int_tmp = np.reshape(self.int_linear[0], (self.int_linear[0].shape[0]*self.int_linear[0].shape[1], self.int_linear[0].shape[2]))
        self.mini_km = MiniBatchKMeans(n_clusters=30).fit(int_tmp) # km_means_list
        if self.verbose:
            print("inertia = ", self.inertia)

    def create_km_map(self, show=True):
        '''
        creates the image of the km_map_datacube
        shows the locations of the k-means clusters for each labels
        outputs: prints the image of the km_map datacube
        '''

        int_linear_shape = np.shape(self.i_3)

        km_map_datacube = np.zeros((int_linear_shape[0], int_linear_shape[1], 30))
        # km_map_datacube = np.zeros((int_linear_shape[0], int_linear_shape[1]))

        for i in range(0, 30): # 30 must be n_cluster from mini_batch_fit
            plt.subplot(5, 6, i+1)
            plt.xlabel('X [DNs]', fontsize=8)
            plt.ylabel('Time [Seconds]', fontsize=8)
            www = np.where(self.k_m.labels_ == i)[0]
            xxx = www/int_linear_shape[1]
            yyy = www%int_linear_shape[1]
            km_map_datacube[xxx.astype(int), yyy.astype(int), i] = i
            #km_map_datacube[xxx.astype(int), yyy.astype(int)] = i
        '''
        plt.figure(figsize=(15, 15))
        plt.subplot(5, 6, i+1)
        plt.xlabel('X [DNs]', fontsize=8)
        plt.ylabel('Time [Seconds]', fontsize=8)
        plt.imshow(km_map_datacube[:, :], extent=[0, 157, 0, 3000], aspect='auto')
        '''
        plt.show()


    def create_profiles(self):
        '''
        creates the k_means maps
        plots the spectral profile for the different k-means labels
        wavelength on the x axis and intensity on the y axis
        outputs: prints the image of the k-means labels
        '''
        # plt.subplots_adjust(bottom=0.15, top=.9, left=0.15)

        wvl = self.data["arr_1"]
        plt.figure(figsize=(20, 20))
        for i in range(0, 30):
            plt.subplot(5, 6, i+1)
            plt.xlabel('Wavelength - 2796.2 [$\AA$]', fontsize=10)
            plt.ylabel('Intensity/I_{c} [$erg s^{-1} cm^{-2} Hz^{-1} ster^{-1}$]', fontsize=10)
            plt.plot(wvl*10.-2795.37, self.k_m.cluster_centers_[i, 0:2000]*1e3/3.454e-6)
        plt.show()

    def create_spectral_map(self):
        '''
        interpolation, creates a spectral profile map showing the
        location of the clusters in all the labels as a whole
        outputs: prints the image of the spectral map
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
