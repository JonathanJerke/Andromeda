gcc   = icc

MKLPATH = $(MKLROOT)/lib/intel64
PYTHONPATH = /opt/ohpc/pub/libs/intel/python3/3.6.4/include/python3.6m

all: olympics
lib =  -Wl,--start-group ${MKLPATH}/libmkl_intel_ilp64.a ${MKLPATH}/libmkl_intel_thread.a ${MKLPATH}/libmkl_core.a ${MKLPATH}/libmkl_blacs_openmpi_ilp64.a  -Wl,--end-group  -liomp5 -lpthread -lm -ldl
flags =-std=c99 -DMKL_ILP64 -m64 -I${MKLROOT}/include

##GSL:: /opt/ohpc/pub/libs/intel/gsl/2.4/lib/libgsl.so /opt/ohpc/pub/libs/intel/gsl/2.4/lib/libgslcblas.so
##GSL:: -I/opt/ohpc/pub/libs/intel/gsl/2.4/include 

olympics: ioPrint.c ioPrint.h eigen.c eigen.h jobs.c input.c input.h constants.h  coreMath.c  coreUtil.h   mAls.c   Model.h coreForce.c  coreMath.h  interfaceMath.c  mAls.h    saUtil.c coreForce.h  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(gcc) $(flags)  -fopenmp jobs.c input.c eigen.c coreMath.c mAls.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c $(lib) -o olympics

adams: ioPrint.c ioPrint.h eigen.c eigen.h main.c input.c input.h constants.h  coreMath.c  coreUtil.h  mAls.c   Model.h coreForce.c  coreMath.h  interfaceMath.c  mAls.h    saUtil.c coreForce.h  input.c  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(gcc) $(flags) main.c input.c eigen.c coreMath.c  mAls.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c $(lib) -o adams
