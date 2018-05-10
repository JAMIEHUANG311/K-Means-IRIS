## NOTES: 

If you full fill a task, add [x], i..e,

 [x] bla

## List files:

   [x] It looks like you did a folder with the same file except that
   the latest one (in the folder) has [*] before each function. This
   is not correct. [] marks are for TODO.dodo document. Let me know
   if this is clear.

   [] Still there are two files (one inside k-means and other 
   in the main directory), delete the oldest one. No need for now of a folder neither. 
   let me know if you need help on reorganize the files in github. 

## TEST:

   [] Create a TESTS folder.

     [] Add scripts that shows the tests.

     [] The tests must use the k-means library, i.e., import k-means

## Style:

   [] I modified one example of your (first one in k-means.py), note the indent
   of the description and how '''''' is added in separate lines. Use this style. 
   for each function.  #. See, for instance, helita/sim/bifrost.py.

## Github:

## Documentation:

   [] README.md: There are options that allows you to put bold, big letter, etc. 
       See one of the examples in helita for instance. 

   [] Consider use github wiki for a more complete documentation
     (at least for the bifrost.md). This can be found here
     (https://github.com/JAMIEHUANG311/K-Means-IRIS/wiki). This will allow
     various sections.

## code:

   [] Add other functions that you did. All the functions that does the interpolation,
   reads the synthetic data, separates between wavelength range etc.

   [x] Try to separate them in class or ordered them by functionality, i.e.,
      1st functions to read, 2nd functions treat the data, 3rd K-min and similar
      4th sanity functions, such as plotting etc.

   [] There are typos in line 16. Only one " in each side (see the code and comments in ### JMS)

   [] time_import is missing self 

   [x] The variables between different functions can use the variables from others
    if they have been storaged in the library, i.e., self.i3 could be use in
    function fit:
      mini_km = MiniBatchKMeans(n_clusters=30).fit(self.i3[:,1000:2000])
    compare with your version.
    Try to find which variables are useful to keep in self for other functions.
    
