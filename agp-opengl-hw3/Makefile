
CPP     = g++ -std=c++11 -Wall --pedantic
NVCC    = /usr/local/cuda-9.1/bin/nvcc -arch=sm_30
INCDIR  = -I./ -I./util -I./modules/build/include
LIBDIR  = -L./ -L./modules/build/lib 
BINDIR  = ./bin
OBJDIR  = $(BINDIR)/obj
LIBS    = -lglad -lGL -lEGL -lXrandr -lXext -lX11 -lrt -ldl -lglut -pthread 
CFLAGS  = $(INCDIR) $(LIBDIR)
MKDIR   = mkdir -p

all: setup main.out

main.out: util.o cuda.o
	@$(CPP) $(CFLAGS) main.cpp -o $(BINDIR)/main.out $(OBJDIR)/*.o $(LIBS) -L/usr/local/cuda-9.1/lib64 -lcudart

# The following line is an example on how to compile your own clases. Use
# the $(OBJDIR) directory to output your object files:
#
# example_class.o:
#	@$(CPP) $(CFLAGS) -c example_class.cpp \
#					  -o $(OBJDIR)/example_class.o
#
# Do not forget to specify example_class.o alongside the main.out label. This
# way, we tell make that main.out depends on example_class.o, triggering the
# build beforehand.
cuda.o:
	@$(NVCC) $(CFLAGS) -c initial.cu -o $(OBJDIR)/cuda.o 

util.o:
	@$(CPP) $(CFLAGS) -c util/util.cpp -o $(OBJDIR)/util.o

setup:
	@$(MKDIR) $(BINDIR) $(OBJDIR)

clean:
	@$(RM) -rf $(BINDIR) $(OBJDIR)
	
rebuild: clean all

