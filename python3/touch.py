#!/usr/bin/env python3


from h5py import File
from sys import argv

File(argv[1],'a').close()
