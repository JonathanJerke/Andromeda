from h5py import File
from sys import argv

File(argv[1],"w").close()
