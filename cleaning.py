#---------------------Calling all the packaged for cleaning---------------------------
import numpy as np 
from astropy.io import fits
from glob import glob
import os
import user_call
import bias_combine
import flat_combine



'''flat_list = input("Enter the list of flat frames")
standard = input("Enter the name of standard star file")

while True:
    arc = input("Do you have separate arc lamp for source and standard star(y/n)")
    if arc=='y':
        src_arc_lamp = input("Enter the arc lamp for source")
        standard_arc_lamp = input("Enter the arc lamp for standard")
        break
    elif arc=='n':
        src_arc_lamp = standard_arc_lamp = input("Enter the lamp")
        break
    else:
        print("Wrong Choice")
        continue
'''

#----------------------------Bias Combine--------------------------------
bias_combine_dict = {
    'median': bias_combine.median,
    'average': bias_combine.average
}
while True:
    while True:
        interactive_mode = input("Do you want to proceed in default mode? (y/n): ")
        if interactive_mode=='y':
            bias_combine_mode_input = input("Enter the desired combining mode(median/average) : ")
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


