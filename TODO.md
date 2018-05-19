##NOTES: 
  If you full fill a task, add [x], i..e,

   [x] bla

## List files:

   [x] It looks like you did a folder with the same file except that
   the latest one (in the folder) has [*] before each function. This
   is not correct. [] marks are for TODO.dodo document. Let me know
   if this is clear.

## TEST:

   [] Create a TESTS folder.

     [] Add scripts that shows the tests.

## Style:

   [x] Dont leave empty lines between the  description of the code and 
   and function name. 

   [x] Make sure that you use the right (same) indent for each function. 

   [x] Make sure that you follow pep8 requirements, e.g., some lines are way 
     too long. 

   [x] rbfilename is not a filename, it is a format. So remove 
   this input parameter and line 28 should be: 

        self.pick_in = open(kmeansfilename, "rb") 

   [x] self.pick_in is not need to store in self. so: 
      
      pick_in = open(kmeansfilename, rbfilename) 

      self.km = pickle.load(pick_in) 

   [x] put an input parameter (verbose) for print statements, i.e., 
      	  __init__ ...  verbose=True
     and each time that you call a print: 
     	 if verbose: 
	    print(...)

## Github:

    [x] Invite Alberto to this github

    	[] Alberto, did you accepted the invitation?

## Documentation:

   [x] Complete README.md with: 

      [x python libraries used

      [x] helita, how to get it and setup

   [x] Do the intro more generic: 
      	 This code allows to read RH (in netCDF format) profiles. 
	 It allows to selects specific lines, interpolate in an uniform 
	 spectra and perform k-mean analysis. 

   [] Consider use github wiki for a more complete documentation
     (at least for the bifrost.md). This can be found here
     (https://github.com/JAMIEHUANG311/K-Means-IRIS/wiki). This will allow
     various sections.

## Code:

   [] Add other functions that you did. All the functions that does the interpolation,
   reads the synthetic data, separates between wavelength range etc. See the code 
   that you have in helita_old/helita/sim and in cutting... 

   [x] Avoid harcoding, e.g., "/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz"
      Instead, introduce a function that reads using an input parameter. You could consider to
      add a default if you like that that one is ok to be hardcode for now:
      __init__(self, data, i3, wvl, ...., npzfilename=""/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz":
      and then self.data = np.load(npzfilename) in line 16.

    [x] Similarly with k-means.pck

    [x] All the functions within a class must have (self) as a input parameter.

    [x] The variables between different functions can use the variables from others
    if they have been storaged in the library, i.e., self.i3 could be use in
    function fit:
      mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i3[:,1000:2000])
    compare with your version.
    Try to find which variables are useful to keep in self for other functions.
   
