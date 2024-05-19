from astropy.io import fits
import os 
import numpy as np

with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
mflat = os.path.join(path, 'mflat.fits')

header = fits.open(mflat)[0].header
xpixels = header['NAXIS1']
ypixels = header['NAXIS2']

#parameters= naverage, function, order, lreject, hreject, iter
def response(data):
    normalised_data = np.mean(data, axis=1)
    pixels = np.arange(0,ypixels,1)
            

    

data = fits.open(mflat)[0].data
while True:
    interactive_node = input("Do you want to change the parameters (y/n): ")
    if interactive_node == 'y':
    