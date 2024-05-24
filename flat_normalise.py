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

#parameters= naverage, function, order, lreject, hreject, iter
def response(data, ypixels):
    parameters = {
        "naverage": 0,
        "function": "spline",
        "order": 3,
        "lreject": 3.0,
        "hreject": 3.0,
        "iter": 1
    }
    function_options = ["spline", "chebyshev", "legendre"]

    #---------------------Entering the parameters--------------------
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
    naverage = parameters['naverage']
    function = parameters['function']
    order = parameters['order']
    lreject = parameters['lreject']
    hreject = parameters['hreject']
    iter = parameters['iter']
    #---------------------------Fitting the data------------------------
    sample_points = []
    while True:
        if abs(int(naverage)) == 0:
            print("naverage cannot be zero. Enter either a positive or negative number")
            naverage = int(input("naverage: "))
            continue
        else:
            break
    sample_pixels = np.arange(0,ypixels,abs(int(naverage)))
    for i in range(0,len(normalised_data),abs(int(naverage))):
        if naverage > 0:
            sample_points.append(np.mean(normalised_data[i:i+int(naverage)]))
        else:
            sample_points.append(statistics.median(normalised_data[i:i+abs(int(naverage))]))
    for j in  range(0,iter):
        if function == "spline":
            spline = UnivariateSpline(sample_pixels, sample_points, k=order)
            y_fit = spline(sample_pixels)
        elif function == "legendre":
            coefficients = np.polynomial.legendre.Legendre.fit(sample_pixels, sample_points, deg=order).convert().coef
            y_fit = np.polynomial.legendre.legval(sample_points, coefficients)
        else:
            coefficients = np.polynomial.chebyshev.Chebyshev.fit(sample_points, sample_pixels, deg=order).convert().coef
            y_fit = np.polynomial.chebyshev.chebval(sample_points, coefficients)
        residuals = sample_points - y_fit
        rms_error = np.sqrt(np.mean(residuals**2))
        std_residuals = np.std(residuals)
        rejected_indices = np.where(np.logical_or(residuals > (np.mean(residuals) + (hreject * std_residuals)), residuals < (np.mean(residuals) - (lreject * std_residuals))))[0]
        sample_pixels = np.delete(sample_pixels, rejected_indices)
        sample_points = np.delete(sample_points, rejected_indices)
        y_fit = np.delete(y_fit, rejected_indices)

    #--------------------------------Plotting the fit--------------------------------------
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_axes([0,0,1,1])
    ax.plot(pixels, normalised_data, c ='red')
    ax.plot(sample_pixels,y_fit, ls='--', c="white")
    ax.grid(True, color='grey')
    ax.set_xlabel("Dispersion axis pixels")
    ax.set_ylabel("Counts")
    ax.set_title("naverage: {}, function: {}, order: {}, lreject:{}, hreject: {}, iter: {}, rms = {}".format(naverage, function, order, lreject, hreject, iter, rms_error))
    plt.show()

    #-----------------------------Checking the fit--------------------------------
    while True:
        fit_check = input("Are you satisfied with the fit ? (y/n): ")
        if fit_check == 'y':
            break
        elif fit_check == 'n':
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
                        if value in function_options:
                            parameters[param] = value
                        else:
                            print("Invalid input. 'function' parameter should be one of the following: " + ", ".join(function_options))
                    else:
                        try:
                            value = float(value)
                            parameters[param] = value
                        except ValueError:
                            print("Invalid input. This parameter should be a number.")
                else:
                    print("Invalid parameter. Please try again.")
            naverage = parameters['naverage']
            function = parameters['function']
            order = parameters['order']
            lreject = parameters['lreject']
            hreject = parameters['hreject']
            iter = parameters['iter']
    
            sample_points = []
            while True:
                if abs(int(naverage)) == 0:
                    print("naverage cannot be zero. Enter either a positive or negative number")
                    naverage = int(input("naverage: "))
                    continue
                else:
                    break
            sample_pixels = np.arange(0,ypixels,abs(int(naverage)))
            for i in range(0,len(normalised_data),abs(int(naverage))):
                if naverage > 0:
                    sample_points.append(np.mean(normalised_data[i:i+int(naverage)]))
                else:
                    sample_points.append(statistics.median(normalised_data[i:i+abs(int(naverage))]))
            for j in  range(0,iter):
                if function == "spline":
                    spline = UnivariateSpline(sample_pixels, sample_points, k=order)
                    y_fit = spline(sample_pixels)
                elif function == "legendre":
                    coefficients = np.polynomial.legendre.Legendre.fit(sample_pixels, sample_points, deg=order).convert().coef
                    y_fit = np.polynomial.legendre.legval(sample_points, coefficients)
                else:
                    coefficients = np.polynomial.chebyshev.Chebyshev.fit(sample_points, sample_pixels, deg=order).convert().coef
                    y_fit = np.polynomial.chebyshev.chebval(sample_points, coefficients)
                residuals = sample_points - y_fit
                rms_error = np.sqrt(np.mean(residuals**2))
                std_residuals = np.std(residuals)
                rejected_indices = np.where(np.logical_or(residuals > (np.mean(residuals) + (hreject * std_residuals)), residuals < (np.mean(residuals) - (lreject * std_residuals))))[0]
                sample_pixels = np.delete(sample_pixels, rejected_indices)
                sample_points = np.delete(sample_points, rejected_indices)
                y_fit = np.delete(y_fit, rejected_indices)
    
            fig = plt.figure(figsize=(7,7))
            ax = fig.add_axes([0,0,1,1])
            ax.plot(pixels, normalised_data, c ='red')
            ax.plot(sample_pixels,y_fit, ls='--', c="white")
            ax.grid(True, color='grey')
            ax.set_xlabel("Dispersion axis pixels")
            ax.set_ylabel("Counts")
            ax.set_title("naverage: {}, function: {}, order: {}, lreject:{}, hreject: {}, iter: {}, rms = {}".format(naverage, function, order, lreject, hreject, iter, rms_error))
            plt.show()

            continue
        else:
            print("wrong Choice")
            continue
    if abs(int(naverage)) != 1:
        y_fit = np.interp(pixels, sample_pixels, y_fit) 
    final_image = np.divide(data, y_fit[:, np.newaxis]) 
    return final_image

data = fits.open(mflat)[0].data
normalised_data = np.mean(data, axis=1)
pixels = np.arange(0,ypixels,1)


output_filename = os.path.join(path,'nmflat.fits')
