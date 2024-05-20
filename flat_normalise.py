from astropy.io import fits
import os 
import numpy as np
import statistics
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from numpy.polynomial import Legendre, Chebyshev
import matplotlib.pyplot as plt
plt.style.use('dark_background')

with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()
mflat = os.path.join(path, 'mflat.fits')

header = fits.open(mflat)[0].header
xpixels = header['NAXIS1']
ypixels = header['NAXIS2']

def spline(x, y, order):
    s = UnivariateSpline(x, y, k=order)
    return s(x)
def legendre(x, y, order):
    l = Legendre.fit(x, y, order)
    return l(x)

def chebyshev(x, y, order):
    c = Chebyshev.fit(x, y, order)
    return c(x)

function_map = {
    'spline': spline,
    'legendre': legendre,
    'chebyshev': chebyshev
}

#parameters= naverage, function, order, lreject, hreject, iter
def response(data, naverage, function,order, lreject, hreject, iter):
    normalised_data = np.mean(data, axis=1)
    pixels = np.arange(0,ypixels,1)
    sample_points = []
    for i in range(0,len(normalised_data),abs(int(naverage))):
        while True:
            if naverage > 0:
                sample_points.append(np.mean(normalised_data[i:i+int(naverage)]))
                break
            elif naverage < 0 :
                sample_points.append(statistics.median(normalised_data[i:i+abs(int(naverage))]))
                break
            else :
                print("Enter either a positive or negative number")
                naverage = input("naverage: ")
                continue
    sample_pixels = np.arange(0,ypixels,abs(int(naverage)))
    popt, pcov = curve_fit(function_map[function], sample_pixels, sample_points)

    
    #while True:
    #    interactive_mode = input("Do yo want to fit inetractively ? (y/n): ")
    #    if interactive_mode == 'y':
            
            

    

data = fits.open(mflat)[0].data
parameters = {
        "naverage": 1,
        "function": "spline",
        "order": 2,
        "lreject": 3.0,
        "hreject": 3.0,
        "iter": 5
    }
function_options = ["spline", "chebyshev", "legendre"]

naverage = parameters['naverage']
function = parameters['function']
order = parameters['order']
lreject = parameters['lreject']
hreject = parameters['hreject']
iter = parameters['iter']

print("Do you want to change the parameters (default values are shown):")
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

        