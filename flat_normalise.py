#---------------------Replica of response task----------------
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
final_image = np.zeros((ypixels, xpixels))
#parameters= naverage, function, order, lreject, hreject, iter
def response(data, ypixels):
    parameters = {
        "naverage": 1,
        "function": "spline",
        "order": 3,
        "lreject": 0,
        "hreject": 0,
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
        if lreject > 0 and hreject > 0:
            rejected_indices = np.where(np.logical_or(residuals > (np.mean(residuals) + (hreject * std_residuals)), residuals < (np.mean(residuals) - (lreject * std_residuals))))[0]
        elif lreject > 0 and hreject ==0:
            rejected_indices = np.where(residuals < (np.mean(residuals) - (lreject * std_residuals)))[0]
        elif lreject == 0 and hreject > 0:
            rejected_indices = np.where(residuals > (np.mean(residuals) + (hreject * std_residuals)))[0]
        elif lreject == 0 and hreject ==0:
            rejected_indices = []    
        for index in sorted(rejected_indices, reverse=True):
            if index < len(sample_points) - 1 and index > 0:
                next_index = index + 1
                while next_index in rejected_indices and next_index < len(sample_points) - 1:
                    next_index += 1
                prev_index = index - 1
                while prev_index in rejected_indices and prev_index > 0:
                    prev_index -= 1
                mean_point = (sample_points[prev_index] + sample_points[next_index]) / 2
                sample_points[index] = mean_point
            elif index == 0:
                next_index = index + 1
                while next_index in rejected_indices and next_index < len(sample_points) - 1:
                    next_index += 1
                sample_points[index] = sample_points[next_index]
            elif index == len(sample_points) - 1:
                prev_index = index - 1
                while prev_index in rejected_indices and prev_index > 0:
                    prev_index -= 1
                sample_points[index] = sample_points[prev_index]
    
    #--------------------------------Plotting the fit--------------------------------------
    plt.figure(figsize=(7,7), facecolor='black')
    plt.plot(pixels, normalised_data, c ='red')
    plt.plot(sample_pixels,y_fit, ls='--', c="white")
    plt.grid(True, color='grey')
    plt.xlabel("Dispersion axis pixels", color='white')
    plt.ylabel("Counts", color='white')
    plt.title("naverage: {}, function: {}, order: {}, lreject:{}, hreject: {}, iter: {}, rms = {}".format(naverage, function, order, lreject, hreject, iter, rms_error, color='white'))
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
                if lreject > 0 and hreject > 0:
                    rejected_indices = np.where(np.logical_or(residuals > (np.mean(residuals) + (hreject * std_residuals)), residuals < (np.mean(residuals) - (lreject * std_residuals))))[0]
                elif lreject > 0 and hreject ==0:
                    rejected_indices = np.where(residuals < (np.mean(residuals) - (lreject * std_residuals)))[0]
                elif lreject == 0 and hreject > 0:
                    rejected_indices = np.where(residuals > (np.mean(residuals) + (hreject * std_residuals)))[0]
                elif lreject == 0 and hreject ==0:
                    rejected_indices = [] 
                for index in sorted(rejected_indices, reverse=True):
                    if index < len(sample_points) - 1 and index > 0:
                        next_index = index + 1
                        while next_index in rejected_indices and next_index < len(sample_points) - 1:
                            next_index += 1
                        prev_index = index - 1
                        while prev_index in rejected_indices and prev_index > 0:
                            prev_index -= 1
                        mean_point = (sample_points[prev_index] + sample_points[next_index]) / 2
                        sample_points[index] = mean_point
                    elif index == 0:
                        next_index = index + 1
                        while next_index in rejected_indices and next_index < len(sample_points) - 1:
                            next_index += 1
                        sample_points[index] = sample_points[next_index]
                    elif index == len(sample_points) - 1:
                        prev_index = index - 1
                        while prev_index in rejected_indices and prev_index > 0:
                            prev_index -= 1
                        sample_points[index] = sample_points[prev_index]
    
    
            plt.figure(figsize=(7,7), facecolor='black')
            plt.plot(pixels, normalised_data, c ='red')
            plt.plot(sample_pixels,y_fit, ls='--', c="white")
            plt.grid(True, color='grey')
            plt.xlabel("Dispersion axis pixels", color='white')
            plt.ylabel("Counts", color='white')
            plt.title("naverage: {}, function: {}, order: {}, lreject:{}, hreject: {}, iter: {}, rms = {}".format(naverage, function, order, lreject, hreject, iter, rms_error, color='white'))
            plt.show()


            continue
        else:
            print("wrong Choice")
            continue
    if abs(int(naverage)) != 1:
        y_fit = np.interp(pixels, sample_pixels, y_fit) 
    
    for i in range(0,ypixels):
        for j in range(0,xpixels):
            final_image[i,j] = data[i,j]/y_fit[i] 
    return final_image

data = fits.open(mflat)[0].data
normalised_data = np.mean(data, axis=1)
pixels = np.arange(0,ypixels,1)


output_filename = os.path.join(path,'nmflat.fits')





'''Spline fit can be done only upto order 5. Legendre fit is not working and have not yet tested the chebyshev fit. Incorporate odd naverage'''