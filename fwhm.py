#!/usr/bin/env python

import pyfits
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import glob, os

pix_size = 4.54
def fits_file(file):
    hdu = pyfits.open(file)
    #print(hdu.info())
    #print hdu[0].data
    #print(hdu[0].header)
    #plt.imshow(hdu[0].data)
    #plt.show()
    data = hdu[0].data
    hdu.close()
    return data

dir = '/Users/parkerf/Documents/FRD_testing/jul_14/dataset_2/'
#z0 = fits_file('/Users/parkerf/Documents/FRD_testing/jul_9/cleaved_fiber/z5.fit')
#bg = fits_file('/Users/parkerf/Documents/FRD_testing/jul_9/cleaved_fiber/bg/bg.fit')

#data = z0-bg

    
#plt.plot(data[800])
#plt.show()
#sys.exit()

class Gauss_fit_row():
    def __init__(self,data_row):
        self.start = 0
        self.end = len(data_row)
        self.data = data_row
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
    Amp = []
    floor = np.mean(y[0][len(y[0])-20:len(y[0])-10])
    for i in range(len(y)):
        y[i] = y[i]-floor
        ind = np.where(y[i]==max(y[i]))
        #print(ind,max(y[i]))
        x0 = max(x[i][ind])
        #print(x0,max(y[i]))
        popt = curve_fit(gauss,x[i],y[i],p0=[max(y[i]),x0,1])[0]
        FWHM.append(2*np.sqrt(2*np.log(2))*popt[2])
        #plt.plot(x[i],y[i])
        #plt.plot(x[i],gauss(x[i],*popt))
        #plt.show()
        x00.append(popt[1])
        Amp.append(popt[0])
            
    avg_FWHM = np.average(FWHM)
    #print('Average_FWHM: '+str(avg_FWHM))
            
    return avg_FWHM, x00, Amp

def radii_pixels(x,y):
    popt = gauss_popt(x,y)
    x00 = popt[1]
    avg_FWHM = popt[0]
        
    left_out = int(x00[0]-avg_FWHM)
    left_in = int(x00[0]+avg_FWHM)
    right_in = int(x00[1]-avg_FWHM)
    right_out = int(x00[1]+avg_FWHM)

    #print('Left: '+str(left_out)+', '+str(left_in)+', Right: '+str(right_in)+', '+str(right_out ))
    return left_out, left_in, right_in, right_out

def radius(right,left):
    tot_pix = right - left
    dia = pix_size*tot_pix
    dia = dia/1000.
    rad = dia/2.
    return(rad)


def find_max_ring(data,mid_row):
    Radii = []
    Amp = []
    

    for i, row in enumerate(data):
        #slice = data[i]
        #params = Gauss_fit_row(slice)
        #left_out, left_in, right_in, right_out, amp = radii_pixels(params.x(),params.y())
        ##if amp > 8:
          #  rad = radius(right_out,left_out)
           # Radii.append(rad)
        #else:
        #    Radii.append(0)
        #HOW TO GET RID OF OUTLIERS?
        if ((mid_row-10 <= i <= mid_row+10)):
            slice = data[i]
            params = Gauss_fit_row(slice)
            left_out, left_in, right_in, right_out = radii_pixels(params.x(),params.y())
            Radii.append(radius(right_out,left_out))
        elif i < mid_row-10:
            Radii.append(0)
        elif i > mid_row+10:
            Radii.append(0)
    
    max_rad = max(Radii)
    max_rad_index = Radii.index(max_rad)
    #print(max_rad_index)
    plt.plot(data[max_rad_index])

    return data[max_rad_index], max_rad_index

def find_max_row(data):
    mid = int(data.shape[0]/2)

    diff = .001
    while (diff > 0):
        new_row, new_mid = find_max_ring(data,mid)
        diff = abs(new_mid-mid)
        mid = new_mid

    return new_row



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


os.chdir(dir)


slices = {}
for file in glob.glob("*.fit"):
    name = os.path.splitext(file)[0]
    print(name)
    bg = fits_file(dir+'bg.fit')
    full_image = fits_file(file)
    img_data = full_image-bg
    slice = find_max_row(img_data)
    slices[name] = slice

distance = str(input('do you know the distance from fiber tip to image? (y/n)'))

if distance == 'n':
    print('Z0 parameters............')
    z0 = Gauss_fit_row(slices['z0'])
    left_out, left_in, right_in, right_out = radii_pixels(z0.x(),z0.y())
    z0_radius = radius(right_out,left_out)
    print('z0 radius: '+str(z0_radius))

    print('Z5 parameters............')
    z5 = Gauss_fit_row(slices['z5'])
    left_out, left_in, right_in, right_out = radii_pixels(z5.x(),z5.y())
    z5_radius = radius(right_out,left_out)
    print('z5 radius: '+str(z5_radius))

    print('Distance from fiber to nearest image if moved 5mm between images: ')
    xdist = xdist()
    print(xdist)
elif distance == 'y':
    xdist = float(input('distance (mm): '))

FWHM_list = []
F_IN_list = []

a_list = ['a0','a1','a2','a3','a4','a5']

for angle in a_list:
    print(angle+' parameters............')
    A = Gauss_fit_row(slices[angle])
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
    
    
    
    
    
    # find distance between two peaks 
    # choose slice which has a max value between peaks
    # take that as the mid slice and then can do the fwhm measurements

sys.exit()
row = z0_new.shape[0]/2
plt.plot(range(z0_new.shape[1]),z0_new[row])
plt.show()


sys.exit()


dir = '/Users/parkerf/Documents/FRD_testing/jul_8/20150708_DESI_1/'
print('FWHM of the images from this file will be calculated......'+dir)

distance = str(input('do you know the distance from fiber tip to image? (y/n)'))




#Now, with xdist we should be able to figure out the rest of the files


FWHM_list = []
F_IN_list = []

a_list = ['a0','a1','a2','a3','a4','a5']

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


