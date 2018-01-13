# GPU_project
### Pengxin Cheng Yuxiang Huan Fengkai Wan


## Github usage:
### Clone the git of our project
git clone https://github.com/wanfengkai/GPU_project.git
cd GPU_project

## Usage on final working version
The final working version that you can use is folder: ```GPU```.
The way to run the code is:

***./scripts/build_libs.sh ./scripts/run.sh***

scripts in  ``` /scripts``` are used to excute the program.

```common.hpp``` and ```part.h``` are header files.

```fragment.glsl ``` and ```vertex.glsl``` are shader file for OpenGL

```part_kernel.cu``` is the cuda calculate kernel file where we include tile computation and so on.

```main.cpp``` is the main frame and it calls all the function in order.






 
