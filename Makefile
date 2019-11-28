# compiler settings

COMPILER 		:= gfortran
FLAGS	 	:= -c -Mcuda=8.0 -mp=allcore -Mpreprocess -Mmpi=mpich

DEBUG_FLAGS		:= -g
INCLUDES		:= -I /opt/cgnslib_3.2.1_modified/include -I /usr/include/python2.6
INCLUDES_UTIL		:= -I /opt/cgnslib_3.2.1_modified/include
LINKS			:= -L /opt/hdf5/lib -lhdf5 -L /opt/cgnslib_3.2.1_modified/lib -lcgns
# C compiler settings (use /usr/local/cuda-5.5/bin/nvcc on K40)
C_COMPILER		:= nvcc
C_COMPILER_Rock		:= nvcc --compiler-bindir /usr/bin/gcc
C_COMPILER_Pascal	:= nvcc
C_COMPILER_K40		:= /usr/local/cuda-6.5/bin/nvcc
C_FLAGS			:= -arch=sm_35
C_FLAGS_K10		:= -arch=sm_30

HOME			:= $(shell pwd)
# directories
BINARY_DIR		:= bin
BDIR 			:= build
SOURCE_DIR		:= src

all : $(BDIR)
	find src -type f -iregex ".*\.\(cuf\|f90\)"  > .srcfiles
	@echo "mkdep warnings redirected to /dev/null"
	mkdep .srcfiles > /dev/null
	mv .mkdep_dependencies .deps
	python build.py
# create bin and build directories if they don't already exist
$(BDIR) :
	@mkdir -p $(BINARY_DIR)
	@mkdir -p $(BDIR)

clean:
	rm -rf $(BDIR) $(BINARY_DIR) .srcfiles .mkdep_* .deps
