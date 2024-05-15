#---------------------Calling all the packaged for cleaning---------------------------
import numpy as np 
from astropy.io import fits
from glob import glob
import os



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
import bias
