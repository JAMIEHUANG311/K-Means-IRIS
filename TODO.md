NOTES: If you full fill a task, add [x], i..e,
 [x] bla

List files:
   [] It looks like you did a folder with the same file except that
   the latest one (in the folder) has [*] before each function. This
   is not correct. [] marks are for TODO.dodo document. Let me know
   if this is clear.

TEST:
   [] Create a TESTS folder.
     [] Add scripts that shows the tests.

Style:
   [] The description for each function must be in betweeen """ """ instead
   of #. See, for instance, helita/sim/bifrost.py.

Github:
    [] Invite Alberto to this github


Documentation:
   [] Fill README.md with information related to this library.
   [] Consider use github wiki for a more complete documentation
     (at least for the bifrost.md). This can be found here
     (https://github.com/JAMIEHUANG311/K-Means-IRIS/wiki). This will allow
     various sections.

code:
   [] Add other functions that you did. All the functions that does the interpolation,
   reads the synthetic data, separates between wavelength range etc.
   [] Try to separate them in class or ordered them by functionality, i.e.,
      1st functions to read, 2nd functions treat the data, 3rd K-min and similar
      4th sanity functions, such as plotting etc.
   [] Avoid harcoding, e.g., "/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz"
      Instead, introduce a function that reads using an input parameter. You could consider to
      add a default if you like that that one is ok to be hardcode for now:
      __init__(self, data, i3, wvl, ...., npzfilename=""/net/opal/Volumes/Amnesia/mpi3drun/2Druns/genohm/rain/new_inte1_02.npy.npz":
      and then self.data = np.load(npzfilename) in line 16.
    [] Similarly with k-means.pck
    [] All the functions within a class must have (self) as a input parameter.
    [] The variables between different functions can use the variables from others
    if they have been storaged in the library, i.e., self.i3 could be use in
    function fit:
      mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i3[:,1000:2000])
    compare with your version.
    Try to find which variables are useful to keep in self for other functions.
    [] Do not hesitate to ask me if the bullets above are unclear.
    