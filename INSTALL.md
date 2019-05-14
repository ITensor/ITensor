# Installation Instructions

## Building the Libraries

1. Edit the user configurable options. Start by making a copy 
   of the sample make options file: 

   `cp options.mk.sample options.mk`

   Then begin editing options.mk in your favorite text editor
   and follow the remaining instructions.

2. Set which compiler to use (the `CCCOM` variable). 
   Make sure to include the flag -std=c++17 or similar 
   to enable C++17.
   
3. If you are using a system other than a Mac, edit `PLATFORM`,
   `BLAS_LAPACK_INCLUDEFLAGS` and `BLAS_LAPACK_LIBFLAGS` to reflect the
   type and location of your BLAS/LAPACK libraries. The list of currently
   available platforms is: mkl, lapack, macos, acml
   (for details see matrix/lapack_wrap.h). The `PLATFORM` variable 
   selects which function signature definitions will be used to wrap 
   vendor-specific BLAS/LAPACK fortran calls into C.

4. Finally, at the top level of the library (same directory as this file),
   run the commmand "make" on the command line.
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
called `project_template` under the `tutorial` folder. Copy the `project_template`
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
