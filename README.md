# Rotational Tensor Optimization

Reads *.csv (comma separated value) files with one 2nd rank R6x6 tensor per row.

Assuming that each tensor corresponds to a specific control volume (i.e., domain) of human cancellous bone, the program transforms the tensor in its "most orthotropic" orientation and stores the results in output *.csv files.

Tensors with more than 100% or less than 0.1% of the specified young modulus are discarded.

## Meta Template
Located in: 
```
./datasets/M-ROTO.meta.template
```
For use with previously used data sets:
```
cat ./datasets/M-ROTO.meta.template >> Your_Meta_File.meta
```

## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with Bash or Zsh
* GNU Compiler Collection (GCC), especially with gfortran
* An installation of Open-MPI

The program must be compiled with:
* Global integer kind=64Bit, signed
* Meta-format integer kind=64Bit, signed
* MPI integer kind=32Bit

The installation of Open MPI is simplified with the install script of the repository "Overview" of the biomechanics-hlrs-gebert organization @GitHub.

### Optional: Gnu debugging
* [gdb](https://www.gnu.org/software/gdb/)
* [tmpi](https://github.com/Azrael3000/tmpi)
* [tmux](https://github.com/tmux/tmux/wiki)

## Build
It's tested and recommended to build and run the program as follows. For developing the program on a laptop, "Julius" is the appropriate system.
### Set up the Environment
```vim ./central_src/auxiliaries/system_environments/<system>.sh```
```source ./environment.sh <system>``` 

* Set an architecture/a system
  * Give the absolute base path of your mpi-installation
  * Alternatively give the proper module names of your compute cluster

### Run make:
Build the program:    ```make```
Create documentation: ```make docs```

### Uninstall:
```make clean && rm -r <your program directory>```

## Usage
For example for testing on julius:
```
mpirun ./bin/roto_v1.0.0_x86_64 -np 4 ./datasets/SC00-0_tc_Dev_ctif_G3S11Sig10.meta
```
### Datasets
... are transfered via file exchange and are not pushed into the repository. 

#### \*.stte.\*, Stiffness Tensor Data
The 2nd rank R6x6 stiffness tensors (i.e., stiffness matrices) are store in files with the filename nomenclature:
```
path/basename.stte.<state>
```
With following states:
* covo - Original position of the stiffness tensors in respect to the control volume
* mono - Tensors optimized by a monotropic objective
* orth - Tensors optimized by a orthotropic objective
* ani1 - Tensors optimized by an anisotropic objective (type1)
* ani2 - Tensors optimized by an anisotropic objective (type2)

### External Sources
Plain text headers are parsed via a [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
### Arbitrary
Use this program at your own risk.

