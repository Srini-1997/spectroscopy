#---------------------Calling all the packaged for cleaning---------------------------
import numpy as np 
from astropy.io import fits
from glob import glob
import os
import user_call
import bias_combine

with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
src = os.path.join(path,lines[3].strip())
std = os.path.join(path,lines[4].strip())
src_lamp = os.path.join(path,lines[5].strip())
std_lamp = os.path.join(path,lines[6].strip())


src_data = fits.open(src)[0].data
src_header = fits.open(src)[0].header
xpixels = src_header['NAXIS1']
ypixels = src_header['NAXIS2']

std_data = fits.open(std)[0].data
std_header = fits.open(std)[0].header

src_lamp_data = fits.open(src_lamp)[0].data
src_lamp_header = fits.open(src_lamp)[0].header

std_lamp_data = fits.open(std_lamp)[0].data
std_lamp_header = fits.open(std_lamp)[0].header

#----------------------------Bias Combine--------------------------------

bias_combine_dict = {
    'median': bias_combine.median,
    'average': bias_combine.average
}
while True:
    while True:
        interactive_mode = input("Do you want to proceed in parameter changing mode? (y/n): ")
        if interactive_mode=='y':
            bias_combine_mode_input = input("Enter the desired combining mode(median/average): ")
            break
        elif interactive_mode=='n':
            bias_combine_mode_input = "median"
            break
        else:
            print("Wrong Choice")
            continue
    bias_function_to_call = bias_combine_dict.get(bias_combine_mode_input)
    if bias_function_to_call:
        bias_function_to_call(bias_combine.data_list)
        break 
    else:
        print("Invalid input!")
        continue
fits.writeto(bias_combine.output_filename, bias_combine.final_image, header= bias_combine.header, overwrite=True, output_verify='ignore')


#----------------------------Flat combine---------------------------------
import flat_combine
flat_combine_dict = {
    'median': flat_combine.median,
    'average': flat_combine.average
}
while True:
    if interactive_mode=='y':
        flat_combine_mode_input = input("Enter the desired combining mode(median/average) : ")
    else :
        flat_combine_mode_input = "median"
    flat_function_to_call = flat_combine_dict.get(flat_combine_mode_input)
    if flat_function_to_call:
        flat_function_to_call(flat_combine.data_list)
        break 
    else:
        print("Invalid input!")
        continue
fits.writeto(flat_combine.output_filename, flat_combine.final_image, header= flat_combine.header, overwrite=True, output_verify='ignore')

#-----------------------Normalising the flat-----------------------------------
import flat_normalise
flat_normalise.response(flat_normalise.data, flat_normalise.ypixels)
fits.writeto(flat_normalise.output_filename, flat_normalise.final_image, header = flat_normalise.header, overwrite = True, output_verify='ignore')


#---------------------Dividing the normalised flat------------------------------
def bf_cor(file):
    global ypixels, xpixels
    final_image = np.zeros((ypixels, xpixels))
    mbias = os.path.join(path,'mbias.fits')
    mbias_data = fits.open(mbias)[0].data
    nmflat = os.path.join(path,'nmflat.fits')
    nmflat_data = fits.open(nmflat)[0].data
    header = fits.open(mbias)[0].header
    xpixels = header['NAXIS1']
    ypixels = header['NAXIS2']
    for i in range(0,ypixels):
        for j in range(0,xpixels):
            final_image[i,j] = (file[i,j] - mbias_data[i,j])/nmflat_data[i,j]
    return final_image

src_bf = bf_cor(src_data)
std_bf = bf_cor(std_data)
src_lamp_bf = bf_cor(src_lamp_data)
std_lamp_bf = bf_cor(std_lamp_data)

output_filename_src = os.path.join(path, src.replace('.fits','bf.fits'))
fits.writeto(output_filename_src, src_bf, header = src_header, overwrite = True, output_verify='ignore')

output_filename_std = os.path.join(path, std.replace('.fits','bf.fits'))
fits.writeto(output_filename_std, std_bf, header = std_header, overwrite = True, output_verify='ignore')

output_filename_src_lamp = os.path.join(path, src_lamp.replace('.fits','f.fits'))
fits.writeto(output_filename_src_lamp, src_lamp_bf, header = src_lamp_header, overwrite = True, output_verify='ignore')

output_filename_std_lamp = os.path.join(path, std_lamp.replace('.fits','f.fits'))
fits.writeto(output_filename_std_lamp, std_lamp_bf, header = std_lamp_header, overwrite = True, output_verify='ignore')



