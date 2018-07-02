# makefile for SkyUAM v1.0 

.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

GSL_INC = /home2/srodriguez/data_tools/dev/include
GSL_LIB = /home2/srodriguez/data_tools/dev/lib

HDF5 = /home2/srodriguez/data_tools/hdf5-1.8
LIBSHDF = -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -L$(HDF5)/lib \
	$(HDF5)/lib/libhdf5_hl.a $(HDF5)/lib/libhdf5.a -lz -ldl -lm \
	-Wl,-rpath -Wl,$(HDF5)/lib
INCLHDF = $(HDF5)/include
LGSL = -L$(GSL_LIB) -lgsl -lgslcblas
INCLUDECOM = -I./src -I$(GSL_INC) 
LIBCPU = $(LGSL) -lm $(LIBSHDF)


CC = gcc
CFLAGS = -g -fopenmp -Wall -I./src/lib -L$(GSL_LIB) -I$(GSL_INC) \
	$(DEFINEOPTIONS) -I$(INCLHDF)
SRC = src/reader.c src/utils.c src/main.c src/define.c src/cosmo_tool.c  \
	src/checkparam.c src/alloc_mem.c src/lc_constructor.c src/selection.c \
	src/mask_rand.c src/cosmo_tool.h src/common.h src/define.h
OBJ = src/reader.o src/utils.o src/main.o src/checkparam.o src/alloc_mem.o \
	src/define.o src/cosmo_tools.o src/selection.o src/lc_constructor.o \
	src/mask_rand.o

all: $(OBJ)
	$(CC) $(CFLAGS) -o SUGAR $(OBJ) $(INCLUDECOM) $(LIBCPU)

clean:
	$(RM) $(OBJ)

src/utils.o: src/utils.c src/common.h src/define.h src/cosmo_tools.h
src/reader.o: src/reader.c src/common.h src/define.h src/cosmo_tools.h
src/src/checkparam.o: src/src/checkparam.c src/common.h src/define.h 
src/src/alloc_mem.o: src/src/alloc_mem.c src/common.h src/define.h  
src/main.o: src/main.c src/common.h src/define.h src/cosmo_tools.h  
src/define.o: src/define.c src/common.h src/define.h src/cosmo_tools.h 
src/cosmo_tools.o: src/cosmo_tools.c src/common.h src/define.h src/cosmo_tools.h 
src/lc_constructor.o: src/lc_constructor.c src/common.h src/define.h src/cosmo_tools.h  
src/selection.o: src/selection.c src/common.h src/define.h src/cosmo_tools.h  
src/mask_rand.o: src/mask_rand.c src/common.h src/define.h src/cosmo_tools.h  




