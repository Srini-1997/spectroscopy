#------------------Replica of zerocombine------------------------------
import numpy as np 
import os
from glob import glob
from astropy.io import fits
import user_input

bias_list = input("Enter the bias files")
bias_files = sorted(glob(os.path.join(user_input.path,bias_list)))
print(f"You have {len(bias_files)} bias files")

#reject- none, minmax, ccdclip, crrject, sigclip, avsigclip, pclip
def mimax(nlow, nhigh):
    for bfiles in bias_files:
        data = fits.open(bfiles)
        header = data[0].header
        image = data[0].header
        nl = header['NAXIS1']*header['']