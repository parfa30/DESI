#!/usr/bin/env python
#####################
#  Written by P. Fagrelius     
#  June, 2015
#
#  PURPOSE: Takes output from IDL data reduction to evaluate the FRD performance
#  at different steps in the fiber production process
#  OUTPUTS: F_in vs. FWHM plots and trend of data throughout a process.
#  INPUTS: quick analysis data from IDL
######################

import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def frd_data(directory):
    """
    Takes output of IDL data reduction and plots
    the FWHM as a function of f_in.

    Returns data for a given directory.
    """
    
    #Import data
    data_file = glob.glob(directory+'/quick*.dat')
    data = np.genfromtxt(data_file[0])
    
    print(data_file)
    names = ['f_in', 'FWHM', 'f85', 'f90', 'f95', 'file']

    #Save data in a dictionary and display in prompt
    Data = {}
    for i,name in enumerate(names):
        Data[name] = data[:,i]

    Data = pd.DataFrame(Data)
    print(Data)

    #Remove data point if it seems off.
    #This should be done by looking at the data printed above
    remove = str(input('Take out data point?(y/n) '))
    if remove == 'y':
        rem_ind1 = int(input('line#1 (start at 0): '))
        Data = Data.drop(Data.index[rem_ind1])
        print(Data)

        remove_again = str(input('Another?(y/n)'))
        if remove_again == 'y':
            rem_ind2 = int(input('line#2: '))
            Data = Data.drop(Data.index[rem_ind2])

        elif remove_again == 'n':
            pass
        else:
            print("You didn't select y/n")
        
    elif remove == 'n':
        pass
    else:
        print("You didn't select y/n")

    #Plot FRD
    plt.plot(Data['f_in'],Data['FWHM'],'.',label='Data')
    plt.axvline(x=3.9,color='r',label='f_in = 3.9')
    #plt.xlim(0,2)
    plt.ylabel('FWHM')
    plt.xlabel('f_in')
    plt.title('FRD')

    return Data


def frd_fit(directory):
    """
    Polyfit data to 3rd order to find FRD at 3.9.
    Will possibly have to remove points to get fit to work.

    Saves plot of data and fit.
    Returns fit
    """
    Data = frd_data(directory)

    fit3 = np.poly1d(np.polyfit(Data['f_in'],Data['FWHM'],3))

    #Plot fit
    xp = np.linspace(np.amin(Data['f_in'],axis=0),np.amax(Data['f_in'],axis=0),100)
    plt.plot(xp,fit3(xp),'-',label = 'polyfit')
    plt.legend()

    #Save image of the plot in the directory
    #plt.savefig(directory+'/FWHM_'+os.path.basename(directory))
    
    plt.show()
    return fit3


    
def frd_trend(fin):
    """
    Plots the FWHM at a given f_in for various steps along the process
    """
    
    print("Input directory names in the order you took the data.\
    There are 6 open - if you don't use one of them type 'n'")
    
    dir1 = str(input('directory: '))
    dir2 = str(input('next directory: '))
    dir3 = str(input('next directory: '))
    dir4 = str(input('next directory: '))
    dir5 = str(input('next directory: '))
    dir6 = str(input('next directory: '))

    directories = [dir1, dir2, dir3, dir4, dir5, dir6]

    trend = []
    labels = []
    for dirn in directories:
        if dirn == 'n':
            pass
        else:
            name = os.path.split(dirn)[1]
            labels.append(name)
            fit = frd_fit(dirn)
            print(fit(fin))
            trend.append(fit(fin))
 
    x = np.arange(0,len(trend))
    plt.plot(x, trend,'.')
    plt.xticks(x, labels)
    plot_margin = 0.25
    x0, x1, y0, y1 = plt.axis()
    plt.axis((x0 - plot_margin,
          x1 + plot_margin,
          y0 - plot_margin,
          y1 + plot_margin))
    plt.title('FWHM trend at f_in = '+str(fin))
    plt.show()


#Uncomment if you want to look at one measurement
#plot_directory = str(input('directory: '))
#frd_fit(plot_directory)

#Uncomment if you want to look at the trend at a given f_in.
#It will prompt you for directories.
frd_trend(3.9)

    
    
