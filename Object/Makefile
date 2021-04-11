### from ubuntu
#sudo apt-get update 
#sudo apt-get install -y apt-utils build-essential libatlas-base-dev liblapack-dev liblapacke-dev gfortran libhdf5-serial-dev python3-pip hdf5-tools
#git clone https://github.com/JonathanJerke/Andromeda 
#git clone https://github.com/TACC/launcher.git
###
#

cc   = h5cc

LAPACKlib          = /usr/lib/x86_64-linux-gnu       
LAPACKinc       = /usr/include/x86_64-linux-gnu

flags           =-I${LAPACKinc}   -L${LAPACKlib}
lib             = -lm -lgfortran -latlas -llapack_atlas  -llapacke -lcblas

olympics: Decompose.c Decompose.h Compression.c Compression.h system.h ioPrint.c ioPrint.h eigen.c eigen.h jobs.c input.c input.h constants.h  coreMath.c  coreUtil.h   Model.h coreForce.c  coreMath.h  interfaceMath.c saUtil.c coreForce.h  input.c  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(cc) $(flags) -fopenmp jobs.c input.c eigen.c coreMath.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c Decompose.c Compression.c $(lib) -o olympics

#cp olympics andromeda

libandromeda.so: Decompose.c Decompose.h Compression.c Compression.h system.h ioPrint.c ioPrint.h eigen.c eigen.h jobs.c input.c input.h constants.h  coreMath.c  coreUtil.h  Model.h coreForce.c  coreMath.h  interfaceMath.c  saUtil.c coreForce.h  input.c  coreUtil.c  interfaceMath.h  Model.c  saUtil.h
	$(cc) $(flags) -fopenmp jobs.c input.c eigen.c coreMath.c coreForce.c interfaceMath.c saUtil.c coreUtil.c Model.c ioPrint.c  Decompose.c Compression.c $(lib) -shlib -shared -Wall -fPIC -o libandromeda.so

