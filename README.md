# GPU_project
### Pengxin Cheng Yuxiang Huan Fengkai Wan


Github usage:

first : git clone https://github.com/wanfengkai/GPU_project.git

cd GPU_project


after init

git add file_you_want

git commit -m 'message why you add'

git push

Cuda reference :https://developer.nvidia.com/gpugems/GPUGems3/gpugems3_ch31.html

Cuda intergration:

export PATH=/usr/local/cuda-9.0/bin:$PATH

nvcc -arch=sm_30 -c part.cu -o part.o

g++ main.cpp -o main.out part.o -L/usr/local/cuda-9.0/lib64 -lcudart

latex public:
https://www.overleaf.com/12785258qvjpxpnpjpyx

Here is our final report.

We can edit it online.



 
