<img width=500 src="doc/mapfl_logo.png" alt="MapFL" />  
  
# MapFL: Map Field Lines  

[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
MapFL is a Fortran tool to trace field lines through 3D fields, including generating diagnostics such as open field maps, expansion factors, and Q-maps.

--------------------------------  
   
## HOW TO BUILD MAPFL
  
The included `build.sh` script will take a configuration file and generate a Makefile and build the code.  
The folder `conf` contains example configuration files for various compilers and systems.  
We recommend copying the configuration file closest to your setup and then modifying it to confomr to your compiler and system (such as `HDF5` library paths/flags, compiler flags, etc.).  
  
Given a configure script `conf/my_custom_build.conf`, the build script is invoked as:  
```
> ./build.sh ./conf/my_custom_build.conf
```

--------------------------------  

## HOW TO RUN MAPFL
  
### Setting Input Options  
  
`MapFL` uses a namelist in an input text file.  
The standard name for the input text file is `mapfl.in` but it can be named anything the user wants.
  
A full working input file with all the default parameter options is provided in the file:  
  
`doc/mapfl.in`  
   
A detailed description of each parameter is also given in that file, and (in addition to this README) is the current main documentation of the code.  
  
We have also provided example input file and dataset in the `example/` folder.  

### Launching the Code ###
    
To run `MapFL`, set the desired run parameters into a file called  `mapfl.in`.  

Set the system environment variable `OMP_NUM_THREADS` to the desired number of CPU threads to run on  
(e.g. `export OMP_NUM_THREADS=12`).  

Finally, copy or link the `mapfl` executable into the same directory as the input file and run the command:  
  
`./mapfl`  

--------------------------------




