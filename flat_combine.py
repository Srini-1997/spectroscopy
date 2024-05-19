#------------------Replica of flatcombine------------------------------
import numpy as np 
import os
from glob import glob
from astropy.io import fits


with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
flat_list = lines[2].strip()
flat_files = sorted(glob(os.path.join(path,flat_list)))
print(f"You have {len(flat_files)} flat files")

data = fits.open(flat_files[0])
header = data[0].header
xpixels = header['NAXIS1']
ypixels = header['NAXIS2']

final_image = np.zeros((ypixels, xpixels))
len_flat_files = len(flat_files)
mbias = os.path.join(path,'mbias.fits')
mbias_data = fits.open(mbias)[0].data
    
def median(data_list):
    middle_point = len_flat_files // 2
    for i in range(0,ypixels):
        for j in range(0,xpixels):
            pixels = [data[i, j] for data in data_list]
            pixels = pixels - mbias_data[i, j]
            pixels.sort()
            if len_flat_files % 2 == 1 :
                final_image[i, j] = pixels[middle_point]
            else:
                final_image[i, j] = (pixels[middle_point-1] + pixels[middle_point])/2
    return final_image

def average(data_list):
    for i in range(0,ypixels):
        for j in range(0,xpixels):
            pixels = [data[i, j] for data in data_list]
            pixels = pixels - mbias_data[i,j]
            final_image[i,j] = np.average(pixels)
    return final_image

data_list = [fits.open(ffile)[0].data for ffile in flat_files]
header['NCOMBINE'] = len_flat_files
output_filename = os.path.join(path,'mflat.fits')