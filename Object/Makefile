cc   = gcc

LAPACKlib          = /usr/lib/x86_64-linux-gnu       
LAPACKinc       = /usr/include/x86_64-linux-gnu

flags           =-I${LAPACKinc}   -L${LAPACKlib}
lib             = -lm -lgfortran -latlas -llapack_atlas  -llapacke -lcblas

andromeda: system.h ioPrint.c ioPrint.h eigen.c eigen.h jobs.c input.c input.h constants.h  coreMath.c  coreUtil.h  mAls.c   Model.h coreForce.c  coreMath.h  interfaceMath.c  mAls.h  saUtil.c coreForce.h  input.c  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(cc) $(flags) -fopenmp jobs.c input.c eigen.c coreMath.c mAls.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c $(lib) -o andromeda

