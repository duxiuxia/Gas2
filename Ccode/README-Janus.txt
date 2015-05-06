On OS X, Install brew from http://brew.sh

$ brew tap homebrew/science
$ brew install netcdf
$ brew install openmpi


To compile Janus

$ make -j8
 --> 8 is the number of threads to use compiling (match CPU count)


If netcdf is installed in non-standard location:

$ export CPPFLAGS="-I/path/to/netcdf/include"
$ export LDFLAGS="-L/path/to/netcdf/lib"
$ make -j8