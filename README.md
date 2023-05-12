# EdynCFD
Compressible fluid simulation codes for combustion 

## Run example codes
build and run a sample sod problem
```
$ git clone git@github.com:Keisuke043/EdynCFD.git
$ cd WORK/Examples/CH4_GRI3.0/SOD
$ mkdir BUILD && cd BUILD
$ cmake ../../../..
$ make
$ cd ..
$ ./run.sh
```

## Install libraries
### OSX
Install gcc and hdf5 with package manager Homebrew
```
$ brew install gcc
$ brew install cmake
$ brew install hdf5-mpi
```

### Intel Mac
Link library path in CMakeLists.txt
```CMakeLists.txt
include_directories(/usr/local/Cellar/hdf5-mpi/1.12.1/include)
link_directories(/usr/local/Cellar/hdf5-mpi/1.12.1/lib)
```

### M1 Mac
```CMakeLists.txt
include_directories(/opt/homebrew/Cellar/hdf5-mpi/1.12.1/include)
link_directories(/opt/homebrew/Cellar/hdf5-mpi/1.12.1/lib)
```

### LINUX
### CentOS Linux release 8.5.2111
Install gcc, hdf5 and openmpi with package manager Spack
#### Install Spack
```
$ cd ~
$ git clone https://github.com/spack/spack.git
$ echo export SPACK_ROOT="/home/---user---/spack" >> ~/.bashrc
$ echo export PATH=$SPACK_ROOT/bin:$PATH >> ~/.bashrc
$ echo . /home/---user---/spack/share/spack/setup-env.sh >> ~/.bashrc
$ source ~/.bashrc
```

#### Install libraries with Spack
```
$ spack install gcc@10.2.0
$ spack load gcc@10.2.0
$ spack install hdf5%gcc@10.2.0+fortran+cxx+szip+hl
$ spack install openmpi@4.0.5%gcc@10.2.0
```

#### Governing equations
<img src="https://latex.codecogs.com/gif.latex?\int_a^bf(x)dx" />
<img src="https://latex.codecogs.com/gif.latex?\inline&space;F_s&space;\frac{n}{N}\&space;(n&space;\in&space;\mathbb{N})" />
```
$$
\frac{\partial Q}{\partial t}+\frac{\partial E}{\partial x}=0 \\
Q^{n+1}=Q^n-\frac{\Delta t}{\Delta x}(\tilde{E}_{j+1/2}^n-\tilde{E}_{j-1/2}^n) \\
$$
```







