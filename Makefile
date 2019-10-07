cc   = gcc

LAPACK          = 
LAPACKlib       = ${LAPACK}
LAPACKinc       = ${LAPACK}/LAPACKE/include

GSL             = /usr/local/include/gsl
GSLlib          = ${GSL}/lib
GSLinc          = ${GSL}/include

flags           = -I${GSLinc} -I${LAPACKinc} -L${GSLlib}  -L${LAPACKlib}
lib             = -lm -lgfortran -llapacke -llapack -lgslcblas  -lgsl -lgfortran


andromeda: system.h ioPrint.c ioPrint.h eigen.c eigen.h jobs.c input.c input.h constants.h  coreMath.c  coreUtil.h  mAls.c   Model.h coreForce.c  coreMath.h  interfaceMath.c  mAls.h  saUtil.c coreForce.h  input.c  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(cc) $(flags) -Wall input.c eigen.c coreMath.c mAls.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c -fopenmp jobs.c $(lib) -o andromeda

