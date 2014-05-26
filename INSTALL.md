# Installation Instructions

## Building the Libraries

1. Determine if your compiler supports C++11, the latest version of the
   C++ standard (this is true for Clang v3.0; G++ v4.7 or so; Intel C++ v13.0).

   If your compiler does not support C++11, download and unpack 
   the boost C++ library [http://www.boost.org](http://www.boost.org). 
   You do not have to compile boost, since ITensor only uses the header 
   files. Boost is a free, high quality, open source C++ library.

   If your compiler does support C++11, set the appropriate options in the 
   `options.mk` file described in the next step.

2. Edit the user configurable options. Start by making a copy of the sample
   Makefile options file: 

   `cp options.mk.sample options.mk`

   Then edit options.mk in your favorite text editor.

   2.1. Set USE_CPP11=yes if your compiler supports C++11.

   2.2. Set which compiler to use (the CCCOM variable). 
        Make sure to include the flag -std=c++11 or similar 
        to enable C++11 if available.
   
   2.3. If you are not using C++11, make sure the variable BOOST_DIR 
        points to the location of your unpacked boost library folder.
        This folder should contain another one called "boost" inside.

   2.4. If you are using a system other than a Mac, edit `PLATFORM`,
        `BLAS_LAPACK_INCLUDEFLAGS` and `BLAS_LAPACK_LIBFLAGS` to reflect the
        type and location of your BLAS/LAPACK libraries. The list of currently
        available platforms is: macos, mkl, acml, lapack
        (for details see matrix/lapack_wrap.h). The `PLATFORM` variable 
        selects which function signature definitions will be used to wrap 
        vendor-specific BLAS/LAPACK fortran calls into C.

3. Finally, at the top level of the library (same directory as this file),
   type "make".
   If all goes well, the built library files should appear in the LIBDIR
   folder specified in options.mk.

Note: sometimes the library has issues compiling if the make "-j" flag is used 
(this flag enables parallel compilation on multi-core machines). Try 
disabling it (e.g. explicitly type "make -j 1") if you have compilation 
errors.


## Building the sample and sandbox apps

We have provided sample applications under the "sample" directory. If you 
would like to experiment with these, consider making a copy of this folder 
to keep the original sample codes as a reference (and experiment on the copy).

To build the sample apps, simply 'cd' into the "sample" folder and type 'make'.
To build an individual app type 'make <appname>'.


## Linking your own applications to the libraries

We strongly recommend placing your own client/driver code *outside* the 
ITensor library source directory. The location you choose is up to you. 

To get started quickly on your own driver code, we have put a folder
called "project_template" under the "tutorial" folder. Copy the project_template
folder to your personal software folder then follow the instructions in the
Makefile to customize it.


## Helpful Makefile variables in options.mk

The `options.mk` file at the top level of the ITensor source directory 
defines a number of Makefile variables that you may find useful in writing 
your own Makefiles. To include these variables, at the top of your Makefile 
put the lines

    LIBRARY_DIR=/path/to/itensor
    include $(LIBRARY_DIR)/this_dir.mk
    include $(LIBRARY_DIR)/options.mk

where of course `/path/to/itensor/` should be replaced with the actual 
location of the ITensor source directory on your machine. 

Including the `options.mk` file in this way defines certain useful 
variables such as 

* `ITENSOR_INCLUDEFLAGS`: compiler flags (of the form `-I/folder/name`) specifying paths to
  ITensor header files, Boost header files, and BLAS/LAPACK header files.

* `ITENSOR_LIBDIR`: the path to the lib/ subdirectory of the ITensor source directory

* `ITENSOR_LIBFLAGS`: flags that specify the names of the statically linked ITensor 
  library files, for example <br/> `-litensor -lmatrix -lutilities`.

* `OPTIMIZATIONS`: user-defined compiler optimization flags, such as `-O2`. It can be helpful for these to 
  match the ones used to build ITensor.

* `DEBUGFLAGS`: user-defined compiler debug-mode flags.
