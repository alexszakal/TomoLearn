This is a reminder on how to install CUDA and which versions are compatible with my ancient hw. 

# GeForce GTX650

GK107, Kepler architecture, CUDA compute capability: 3.0

Installation on Ubuntu 20.04.4 LTS
## Install the correct driver
- This setup needs nvidia-driver-470, install it first!
- If X is not loading after reboot for some reason: `sudo apt-get purge nvidia-*`

## Install CUDA Toolkit 10.0
- Install the Toolkit with the runfile method. Start the .sh file with the `--override` switch to ignore the error because of incompatible gcc version
- add /usr/local/cuda-10.0/bin to path (in `.bashrc` script and/or `/etc/environment`)
- include `/usr/local/cuda-10.0/lib64` in `/etc/ld.so.conf` and run `sudo ldconfig`
- Install `gcc-7` and `g++-7` (The default compiler is gcc-9 which is not compatible with nvcc)
- Create symbolic links to gcc-7 and g++-7 in the directory of nvcc: `ln -s /usr/bin/gcc-7 /usr/local/cuda/bin/gcc`  `ln -s /usr/bin/g++-7 /usr/local/cuda/bin/g++`

## Setup of environment variables in CMAKE:
- Set `CUDACXX` environment variable to `/usr/local/cuda-10.0/bin/nvcc`

## Test
- Test the install by compiling the nbody CUDA example

## Visual profiler start
Path to JAVA has to be defined: `nvvp -vm /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java`
