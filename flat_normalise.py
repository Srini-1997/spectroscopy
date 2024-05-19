from astropy.io import fits
import os 
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
mflat = os.path.join(path, 'mflat.fits')

header = fits.open(mflat)[0].header
xpixels = header['NAXIS1']
ypixels = header['NAXIS2']

#parameters= naverage, function, order, lreject, hreject, iter
def response(data, parameters):
    normalised_data = np.mean(data, axis=1)
    pixels = np.arange(0,ypixels,1)
    naverage = parameters['naverage']
    function = parameters['function']
    order = parameters['order']
    lreject = parameters['lreject']
    hreject = parameters['hreject']
    iter = parameters['iter']
    normalised_data = np.mean(data, axis=1)
    
    while True:
        interactive_mode = input("Do yo want to fit inetractively ? (y/n): ")
        if interactive_mode == 'y':
            
            

    

data = fits.open(mflat)[0].data
parameters = {
        "naverage": 1,
        "function": "spline1",
        "order": 2,
        "lreject": 3.0,
        "hreject": 3.0,
        "iter": 5
    }
function_options = ["spline1", "spline3", "chebyshev", "legendre"]

print("Please enter the values for the parameters you want (default values are shown):")
while True:
    for param, default in parameters.items():
        print(f"{param} (default: {default})")
    param = input("Enter the parameter you want to input (or 'done' to finish): ")
    if param.lower() == 'done':
        break
    if param in parameters:
        value = input(f"{param} (default: {parameters[param]}): ")
        if ' ' in value:
            print("Invalid input. This parameter should not contain a space.")
            continue
        if param == 'function':
            # Check if the input is one of the four options
            if value in function_options:
                parameters[param] = value
            else:
                print("Invalid input. 'function' parameter should be one of the following: " + ", ".join(function_options))
        else:
            # Check if the input can be converted to a number
            try:
                value = float(value)
                parameters[param] = value
            except ValueError:
                print("Invalid input. This parameter should be a number.")
    else:
        print("Invalid parameter. Please try again.")

        