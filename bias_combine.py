#------------------Replica of zerocombine------------------------------
import numpy as np 
import os
from glob import glob
from astropy.io import fits
import statistics

with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
bias_list = lines[1].strip()
bias_files = sorted(glob(os.path.join(path,bias_list)))
print(f"You have {len(bias_files)} bias files")

#reject- none, minmax, ccdclip, crrject, sigclip, avsigclip, pclip
#combine - average, median 
data = fits.open(bias_files[0])
header = data[0].header
xpixels = header['NAXIS1']
ypixels = header['NAXIS2']

final_image = np.zeros((ypixels, xpixels))
len_bias_files = len(bias_files)

def median(data_list):
    middle_point = len_bias_files // 2
    for i in range(ypixels):
        for j in range(xpixels):
            pixels = [data[i, j] for data in data_list]
            final_image[i,j] = statistics.median(pixels)
    return final_image

def average(data_list):
    for i in range(0,ypixels):
        for j in range(0,xpixels):
            pixels = [data[i,j] for data in data_list]
            final_image[i,j] = np.average(pixels)
    return final_image

data_list = [fits.open(bfile)[0].data for bfile in bias_files]
header['NCOMBINE'] = len_bias_files
output_filename = os.path.join(path,'mbias.fits')

        