gcc   = gcc

LAPACK = /home/jjerke/lapack-3.8.0/
GSL = /home/jjerke/gsl
OMPFLAG = -fopenmp

flags = -I${GSL}/include -I/${LAPACK}/LAPACKE/include
lib = -lm -lgfortran ${LAPACK}/liblapacke.a ${LAPACK}/liblapack.a  ${GSL}/lib/libgslcblas.a ${GSL}/lib/libgsl.a -lgfortran

andromeda: ioPrint.c ioPrint.h eigen.c eigen.h main.c input.c input.h constants.h  coreMath.c  coreUtil.h   mAls.c   Model.h coreForce.c  coreMath.h  interfaceMath.c  mAls.h    saUtil.c coreForce.h  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(gcc) $(flags) $(OMPFLAG) main.c input.c eigen.c coreMath.c mAls.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c $(lib) -o andromeda
