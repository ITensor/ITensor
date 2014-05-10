# Installation Instructions

## Building the Libraries

1. Download and unpack the [boost C++ library](http://www.boost.org). You do not have
   to compile it, since ITensor only uses the header files. 
   Boost is a free, high quality, open source C++ library used 
   extensively by ITensor and potentially helpful for your own application development.

2. Edit the user configurable options. Start by making a copy of the sample
   options file: 

   `cp options.mk.sample options.mk`

   Then edit options.mk in  your favorite text editor.

   You must make sure that BOOST_DIR points to the location of
   your unpacked boost libraries folder (the default is a folder called
   'boost' located one directory above this file).

   If you are using a system other than a Mac, edit `PLATFORM`,
   `BLAS_LAPACK_INCLUDEFLAGS` and `BLAS_LAPACK_LIBFLAGS` to reflect the
   type and location of your BLAS/LAPACK libraries. The list of currently
   available platforms is: macos, mkl, acml (see matrix/lapack_wrap.h). 
   The `PLATFORM` variable selects which function signature definitions will be used to 
   wrap vendor-specific BLAS/LAPACK fortran calls into C.

3. Finally, at the top level of the library (same directory as this file),
   type "make".
   If all goes well, the built library files should appear in the LIBDIR
   folder specified in options.mk.

Note: sometimes the library has issues compiling if the make "-j" flag is used 
(this flag enables parallel compilation on multi-core machines). Try 
disabling it (e.g. explicitly type "make -j 1") if you have compilation 
errors.

## Building the sample and sandbox apps

We have provided sample applications under the "sample" directory. Also, we
have provided a "sandbox" folder with some apps for you to modify to
experiment with the library's features.

To build these apps, simply 'cd' into each folder and type 'make'. Or, to 
build an individual app type 'make <appname>'.

## Linking your own applications to the libraries

We strongly recommend placing your own client/driver code *outside* the ITensor library source directory.
The location you choose is up to you. For sample Makefiles, the one in the sample can be a useful template.

The `options.mk` file at the top level of the ITensor source directory defines a number of Makefile 
variables that you may find useful in writing your own Makefiles. To include these variables,
at the top of your Makefile put the lines

    LIBRARY_DIR=/path/to/itensor
    include $(LIBRARY_DIR)/this_dir.mk
    include $(LIBRARY_DIR)/options.mk

where of course `/path/to/itensor/` should be replaced with the actual location of the ITensor source
directory on your machine. 

Including the `options.mk` file in this way defines certain useful variables such as 

* `ITENSOR_INCLUDEFLAGS`: compiler flags (of the form `-I/folder/name`) specifying paths to
  ITensor header files, Boost header files, and BLAS/LAPACK header files.

* `ITENSOR_LIBDIR`: the path to the lib/ subdirectory of the ITensor source directory

* `ITENSOR_LIBFLAGS`: flags that specify the names of the statically linked ITensor 
  library files, for example <br/> `-litensor -lmatrix -lutilities`.

* `OPTIMIZATIONS`: user-defined compiler optimization flags, such as `-O2`. It can be helpful for these to 
  match the ones used to build ITensor.

* `DEBUGFLAGS`: user-defined compiler debug-mode flags.
