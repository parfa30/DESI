#!/usr/bin/env python

import pyfits
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys


pix_size = 4.54

dir = '/Users/parkerf/Documents/FRD_testing/jul_14/dataset_2/'
print('FWHM of the images from this file will be calculated......'+dir)

distance = str(input('do you know the distance from fiber tip to image? (y/n)'))


class Gauss_fit():
    def __init__(self,file_name,directory):
        self.dir = directory
        SE_file = np.genfromtxt(self.dir+'start_end.txt',names=True)
        self.start = SE_file[file_name][0]
        self.end = SE_file[file_name][1]
        self.file = self.dir+'hozslice_'+file_name+'.dat'
        self.data = np.genfromtxt(self.file)
        self.data = self.data[:,1]
        #print(len(self.data),self.end-self.start)
        self.middle_pix = int((self.end-self.start+1)/2.)
        #print(self.middle_pix)

    def y(self):
        yy = [self.data[0:self.middle_pix-1],self.data[self.middle_pix:]]
        return yy
    def x(self):
        xx = [np.arange(self.start,self.start+len(self.y()[0]),1),np.arange(self.start+len(self.y()[0])+1,self.start+len(self.y()[0])+1+len(self.y()[1]),1)]
        return xx

def gauss(xx,a,b,c):
        return a*np.exp(-(xx-b)**2/(2*c**2))    
    
def gauss_popt(x,y):
    FWHM = []
    x00 = []
    floor = np.mean(y[0][len(y[0])-20:len(y[0])-10])
    for i in range(len(y)):
        y[i] = y[i]-floor
        ind = np.where(y[i]==max(y[i]))
        #print(ind,max(y[i]))
        x0 = x[i][ind]
        #print(x0)
        popt = curve_fit(gauss,x[i],y[i],p0=[max(y[i]),x0,1])[0]
        FWHM.append(2*np.sqrt(2*np.log(2))*popt[2])
        plt.plot(x[i],y[i])
        plt.plot(x[i],gauss(x[i],*popt))
        #plt.show()
        x00.append(popt[1])
            
    avg_FWHM = np.average(FWHM)
    #print('Average_FWHM: '+str(avg_FWHM))
            
    return avg_FWHM, x00

def radii_pixels(x,y):
    popt = gauss_popt(x,y)
    x00 = popt[1]
    avg_FWHM = popt[0]
        
    left_out = int(x00[0]-avg_FWHM)
    left_in = int(x00[0]+avg_FWHM)
    right_in = int(x00[1]-avg_FWHM)
    right_out = int(x00[1]+avg_FWHM)

    print('Left: '+str(left_out)+', '+str(left_in)+', Right: '+str(right_in)+', '+str(right_out ))
    return left_out, left_in, right_in, right_out

def radius(right,left):
    tot_pix = right - left
    dia = pix_size*tot_pix
    dia = dia/1000.
    rad = dia/2.
    return(rad)
    
def xdist():
    y2 = z5_radius
    y1 = z0_radius
    xdist = 5*y1/(y2-y1)
    return xdist

def rtod(radian):
    return radian*(180/np.pi)

def FWHM(outer_rad,inner_rad):
    return rtod(np.arctan(outer_rad/xdist) - np.arctan(inner_rad/xdist))

def angle(radius):
    return rtod(np.arctan(radius/xdist))

def fin_radius(outer_rad,inner_rad):
    return (outer_rad-inner_rad)/2 + inner_rad
    
def F_IN(radius):
    return xdist/(2*radius)

###Calculate Distance
if distance == 'n':
    print('Z0 parameters............')
    z0 = Gauss_fit('z0',dir)
    left_out, left_in, right_in, right_out = radii_pixels(z0.x(),z0.y())
    z0_radius = radius(right_out,left_out)
    print('z0 radius: '+str(z0_radius))

    print('Z5 parameters............')
    z5 = Gauss_fit('z5',dir)
    left_out, left_in, right_in, right_out = radii_pixels(z5.x(),z5.y())
    z5_radius = radius(right_out,left_out)
    print('z5 radius: '+str(z5_radius))

    print('Distance from fiber to nearest image if moved 5mm between images: ')
    xdist = xdist()
    print(xdist)
    
elif distance == 'y':
    xdist = float(input('distance (mm): '))

#Now, with xdist we should be able to figure out the rest of the files


FWHM_list = []
F_IN_list = []

a_list = ['a0','a1','a2','a3','a4','a5','a6']

for angle in a_list:
    print(angle+' parameters............')
    A = Gauss_fit(angle,dir)
    lo, li, ri, ro = radii_pixels(A.x(),A.y())
    outer_rad, inner_rad = radius(ro,lo),radius(ri,li)
    FWHM_list.append(FWHM(outer_rad,inner_rad))
    mid_radius = fin_radius(outer_rad,inner_rad)
    F_IN_list.append(F_IN(mid_radius))

print('OUTPUT.................')


for i in range(len(FWHM_list)):
    print(a_list[i],'f_in,  FWHM')
    print(F_IN_list[i],FWHM_list[i])
        
#print(FWHM_list)
#print(F_IN_list)

plt.figure()
plt.plot(F_IN_list,FWHM_list,'.')
plt.xlabel('F_in')
plt.ylabel('FWHM')
plot_margin = 0.25
x0, x1, y0, y1 = plt.axis()
plt.axis((x0 - plot_margin,
          x1 + plot_margin,
          y0 - plot_margin,
          y1 + plot_margin))
plt.show()

#data1 = [FWHM_list,F_IN_list]
#Data = pd.DataFrame(data1,columns = ['FWHM','F_in'], index = {'a0','a1'})
#Data = Data.rename(index = {1:'a0',2:'a1',3:'a2',4:'a3',5:'a4',6:'a5'})
#print(Data)


