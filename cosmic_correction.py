import astroscrappy
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.visualization import ZScaleInterval, ImageNormalize, SqrtStretch
plt.style.use('dark_background')


with open('user_input.txt', 'r') as file:
    lines = file.readlines()

path = lines[0].strip()

def cosmic(data):
    parameters = {
        "sigclip": 4.5,
        "sigfrac": 0.3,
        "objlim": 5.0,
        "satlevel": np.inf,
    }
    sigclip = parameters["sigclip"]
    sigfrac = parameters["sigfrac"]
    objlim = parameters["objlim"]
    satlevel = parameters["satlevel"]

    with open('user_input.txt', 'r') as file:
        lines = file.readlines()

    readnoise = float(lines[7].strip())
    gain = float(lines[8].strip())
    cr_corrected = astroscrappy.detect_cosmics(data, sigclip=sigclip, sigfrac = sigfrac, objlim = objlim, gain = gain, readnoise=readnoise, satlevel= satlevel)[1]

    interval = ZScaleInterval()
    zscale_min, zscale_max = interval.get_limits(cr_corrected)
    norm = ImageNormalize(vmin=zscale_min, vmax=zscale_max, stretch=SqrtStretch())

    #--------------------------------Plotting the fit--------------------------------------
    plt.figure(figsize=(3,7), facecolor='black')
    plt.subplot(1,2,1)
    plt.imshow(data, norm=norm, cmap='gray', origin='lower')
    plt.tick_params(axis='both', labelleft= False, labelbottom = False, left=False, bottom=False)
    plt.subplot(1,2,2)
    plt.imshow(cr_corrected, norm=norm, cmap='gray', origin='lower')
    plt.tick_params(axis='both', labelleft= False, labelbottom = False, left=False, bottom=False)
    plt.show()

    while True:
        cosmic_check = input("Are you satisfied with the cosmic ray correction(y/n): ")
        if cosmic_check == 'y':
            break
        elif cosmic_check == 'n':
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
                    try:
                        value = float(value)
                        parameters[param] = value
                    except ValueError:
                        print("Invalid input. This parameter should be a number.")
                else:
                    print("Invalid parameter. Please try again.")
            sigclip = parameters["sigclip"]
            sigfrac = parameters["sigfrac"]
            objlim = parameters["objlim"]
            satlevel = parameters["satlevel"]

            cr_corrected = astroscrappy.detect_cosmics(data, sigclip=sigclip, sigfrac = sigfrac, objlim = objlim, gain = gain, readnoise=readnoise, satlevel= satlevel)[1]

            interval = ZScaleInterval()
            zscale_min, zscale_max = interval.get_limits(cr_corrected)
            norm = ImageNormalize(vmin=zscale_min, vmax=zscale_max, stretch=SqrtStretch())

            #--------------------------------Plotting the fit--------------------------------------
            plt.figure(figsize=(3,7), facecolor='black')
            plt.subplot(1,2,1)
            plt.imshow(data, norm=norm, cmap='gray', origin='lower')
            plt.tick_params(axis='both', labelleft= False, labelbottom = False, left=False, bottom=False)
            plt.subplot(1,2,2)
            plt.imshow(cr_corrected, norm=norm, cmap='gray', origin='lower')
            plt.tick_params(axis='both', labelleft= False, labelbottom = False, left=False, bottom=False)
            plt.show()

            continue
        else:
            print("wrong choice")
            continue
    
    return cr_corrected
