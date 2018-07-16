'''
This allows to handle something
JMS Explain things here.
'''
import pickle
import time
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
        print("True") #JMS remove this print statement.
        '''
        if self.verbose:
            print(some meaningfull message), message, warning.
        '''


    def read_data_ncdf(self, rbfilename='output_ray_l2d90x40r.ncdf'): ## ncdfilename
        '''
        Reading from the original RH code ncdf files
        '''

        self.wvl1 = nd.getvar(rbfilename, 'wavelength') ## wvl1 -> wvl_orig
        self.inte = nd.getvar(rbfilename, 'intensity', memmap=True) ## int_orig


    def read_data_pck(self, kmeansfilename='k-means.pck'):
        '''
        def save_all_pkl(self,filename=''):
            data = self.__dict__
            filehandler = open(filename+'.opy', 'wb')
            pickle.dump(data, filehandler,protocol=4)
        def load_all_pkl(self,filename):
            filehandler = open(filename, 'rb')
            data=pickle.load(filehandler)
            for idata in data.keys():
                setattr(self,idata,data[idata])
        '''
        pick_in = open(kmeansfilename, 'rb')
        self.k_m = pickle.load(pick_in)


    def wavelength_distinction(self, delta=5): # spect_lines_limits
        '''
        This allows to cut the wavelength range in different spectral lines.
        Find where the value of delwvl is larger than the previous index
        by delta, appends that value into limits array
        Loops through the limits array, and sees if the value for one index
        in limits has a greater difference than 1 than the previous index.
        JMS follow standards of documentation.
        http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
        '''

        delwvl = self.wvl1-shift(self.wvl1, 1) ### delwvl_up
        new_delwvl = shift(self.wvl1, -1) - self.wvl1  ### delwvl_dn
        self.limits = [v for v in range(1, len(delwvl)) if np.abs(delwvl[v]-new_delwvl[v]) > delta] ### wvl_lmts
        self.limits.append(len(self.wvl1))

    def print_list_wvl(self, delta=5): # print_line_list
        '''
        The profile at inte is added to the array new_inte
        '''

        self.new_inte = {}  #Remove this line.

        if not hasattr(self, 'limits'):
            self.wavelength_distinction(delta)

        count = 0
        for ind in range(0, np.size(self.limits)-1):
            if self.limits[ind+1]-self.limits[ind] > 1:
                print('%i,wvl=%0.1f' %(count, self.wvl1[self.limits[ind]]))
                count = count + 1


    def individual_spectral_data(self, delta=5): # break_lines, delta->wvl_delta
        '''
        The profile at inte is added to the array new_inte
        '''
        if not hasattr(self, 'new_inte'):
            self.new_inte = {} ### int_indv
        if not hasattr(self, 'new_wvl'):
            self.new_wvl = {} ### wvl_indv
            print('True')

        self.new_inte = {} #remove
        self.new_wvl = {} #remove

        if not hasattr(self, 'limits'):
            self.wavelength_distinction(delta)

        count = 0
        for ind in range(0, np.size(self.limits)-1):
            if self.limits[ind+1]-self.limits[ind] > 1:
                self.new_inte[count] = self.inte[:, :, self.limits[ind]:self.limits[ind+1]-1] # var
                self.new_wvl[count] = self.wvl1[self.limits[ind]:self.limits[ind+1]-1]
                count = count + 1

    def km_interp(self, ind=0, savefile=None, nwvl_min=500): #linear_spect, ind->id_lines,
        '''
        interpolation of the axis
        plots wvl against new_inte in an uniform axis(wvlax)
        '''

        self.new_inte1 = {}

        if not hasattr(self, 'wvlax'): # wvl_linear
            self.wvlax = {}

        delwvl = np.gradient(self.new_wvl[ind])

        if not hasattr(self, 'new_inte1'): # int_linear
            self.new_inte1 = {}

        mindelwvl = np.min(delwvl)
        n_points = np.min([(np.max(self.new_wvl[ind])-np.min(self.new_wvl[ind]))/mindelwvl, nwvl_min])


        inte1 = {} # remove

        if n_points > 5:
            max_value = np.max(self.new_wvl[ind])
            min_value = np.min(self.new_wvl[ind])
            self.wvlax[ind] = np.linspace(min_value, max_value, num=n_points)
            inte1 = np.zeros((self.inte.shape[0],
                              self.inte.shape[1], np.shape(self.wvlax[ind])[0]))

            for ind2 in range(0, len(self.new_inte[ind][:, 0, 0])):
                for ind3 in range(0, len(self.new_inte[ind][0, :, 0])):
                    inte1[ind2, ind3, :] = np.interp(self.wvlax[ind],
                                            self.new_wvl[ind],
                                            self.new_inte[ind][ind2, ind3, :])

        self.new_inte1[ind] = inte1
        '''
        use verbose for print messages, make a more meaningfull messeage, i.e.,
        no new_inte1  but something else.
        '''
        print('new_inte1', ind, np.shape(self.new_inte1[ind]))
        if savefile != None:
            mg2_cube = self.new_inte1[ind]
            np.save(savefile + '.npy', mg2_cube)


    def read_lim_interp(self, rbfilename='/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/output_ray_l2d90x40r.ncdf', ind=0, delta=5):
        self.read_data_ncdf(rbfilename=rbfilename)
        self.individual_spectral_data(delta=delta)
        self.km_interp(ind=ind)



    def time_import(self, maxnum=5): # km_clusters_???? max_cluster_niter
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters,
        computes the inertia of the MiniBatchKMeans
        outputs: graph of t_m and inertia
        '''

        i3_2d = self.i_3.reshape(self.i_3.shape[0]*self.i_3.shape[1], self.i_3.shape[2]) ## int_linear[line_id]
        self.t_m = np.zeros(maxnum) # remove
        self.inertia = np.zeros(maxnum) #km_inertia
        self.t_zero = time.time() #remove
        self.t_mini_batch = time.time() - self.t_zero #  remove
        for i in range(0, maxnum):
            self.mini_km = MiniBatchKMeans(n_clusters=(i+1)*10, n_init=10).fit(i3_2d)
            self.t_m[i] = self.t_mini_batch/((i+1)*10) # remove
            self.inertia[i] = self.mini_km.inertia_
            ### save in an array (i+1)*10
            ### use verbose if statement.
            print('description iteration=',i)
        '''
        plt.plot(the array that contains the infor of n_clusters,self.inertia)
        plt.show()
        '''
        print(self.t_m)
        return self.mini_km #Remove

    def mini_batch_fit(self, ind=0):
        '''
        uses the MiniBatchKMeans function to fit the i3_2d data into clusters
        outputs: prints time and inertia
        '''

        self.t_zero = time.time() #remove

        self.a = np.reshape(self.new_inte1[0], (self.new_inte1[0].shape[0]*self.new_inte1[0].shape[1], self.new_inte1[0].shape[2])) # a -> int_tmp (no self)
        self.mini_km = MiniBatchKMeans(n_clusters=30).fit(self.a) # km_means_list
        print("time = ", self.t_mini_batch) #remove
        print("inertia = ", self.inertia) #use verbose

    def create_km_map(self,show=True):
        '''
        creates the image of the km_map_datacube
        shows the locations of the k-means clusters for each labels
        outputs: prints the image of the km_map datacube
        '''

        dim_i_3 = np.shape(self.i_3) # int_linear_shape

        km_map_datacube = np.zeros((dim_i_3[0], dim_i_3[1], 30))
        # km_map_datacube = np.zeros((dim_i_3[0], dim_i_3[1]))

        for i in range(0, 30): # 30 must be n_cluster from mini_batch_fit
            plt.subplot(5, 6, i+1)
            plt.xlabel('X [DNs]', fontsize=8)
            plt.ylabel('Time [Seconds]', fontsize=8)
            www = np.where(self.k_m.labels_ == i)[0]
            xxx = www/dim_i_3[1]
            yyy = www%dim_i_3[1]
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
