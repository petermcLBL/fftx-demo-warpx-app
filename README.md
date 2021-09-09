WarpX-like Interface to FFTX: A Demo
====================================

This application calls a built FFTX library using an interface like that in WarpX.

### Building Demo FFTX External Application

To use and build this demo application you must have installed FFTX and spiral-software. Ensure your environment sets **FFTX_HOME** and **SPIRAL_HOME** to point
to the locations of FFTX and Spiral.

### Installing Pre-requisites

Clone **spiral-software** to a location on your computer.  E.g., do:
```
cd ~/work
git clone https://www.github.com/spiral-software/spiral-software
```
This location is known as *SPIRAL HOME* and you must set an environment variable
**SPIRAL_HOME** to point to this location later.

You must also install two spiral packages, do the following:
```
cd ~/work/spiral-software/namespaces/packages
git clone https://www.github.com/spiral-software/spiral-package-fftx fftx
git clone https://www.github.com/spiral-software/spiral-package-simt simt
```
**NOTE:** The spiral packages must be installed under the directory
**$SPIRAL_HOME/namespaces/packages** and must be placed in folders with the
prefix *spiral-package* removed. 

Follow the build instructions for **spiral-software** (see the **README**
[**here**](https://github.com/spiral-software/spiral-software/blob/master/README.md) ).

You should be using the `develop` branch of all of these repos:
```
pushd $SPIRAL_HOME
git checkout develop
pushd namespace/packages/fftx
git checkout develop
pushd ../simt
git checkout develop
popd
popd
```

### Installing FFTX

Clone **FFTX** to a location on your computer.  E.g., do:
```
cd ~/work
git clone https://www.github.com/spiral-software/fftx
```

Then **FFTX_HOME** should be set to `~/work/fftx` . Use the `develop` branch:
```
pushd $FFTX_HOME
git checkout develop
```

In order to build the FFTX _library_, you first need to generate files with `./create_lib_code.sh` .  It is set up to work with CUDA, but you can have it work with HIP by simply changing references to "cuda" in that script to "hip".
```
cd $FFTX_HOME/examples/library
./create_lib_code.sh
```

Now follow the build instructions for **FFTX** (see the **README**
[**here**](https://github.com/spiral-software/FFTX/blob/master/README.md) ).

### Install and Build the Demo Application

Ensure you have valid settings for **FFTX_HOME** and **SPIRAL_HOME**.

This application also needs **AMREX_SRC** and **AMREX_BUILD** to be set to AMReX source directories.  For example:
```
setenv AMREX_SRC $HOME/warpx_directory/WarpX/build/_deps/fetchedamrex-src/Src/Base
setenv AMREX_BUILD $HOME/warpx_directory/WarpX/build/_deps/fetchedamrex-build
```

Clone the demo application.  E.g., do:
```
cd ~/work
git clone https://www.github.com/spiral-software/fftx-demo-warpx-app
cd fftx-demo-warpx-app
cmake -S . -B build
cmake --build build --target install
```
The demo application will be installed at `~/work/fftx-demo-warpx-app/build/bin/testwarpx` .
