##after building an ubuntu Virtual Machine 
#git clone https://github.com/JonathanJerke/Andromeda 
# then execute this first...

sudo apt-get update 
sudo apt-get install -y apt-utils build-essential libatlas-base-dev liblapack-dev liblapacke-dev gfortran libhdf5-serial-dev python3-pip hdf5-tools hwloc util-linux
pip3 install pandas
pip3 install numpy
pip3 install Cython
pip3 install h5py
pip3 install scipy
##compile Andromeda

source compile.sh
