import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.io import fits
from astropy.visualization import ZScaleInterval
from scipy.interpolate import interp1d

def tryMakePath(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        print('Directory {} already exists'.format(path))
    else:
        print('Directory {} created'.format(path))

def makeDirectoryStructure(DATA_parent, obj_name, readmodes):
    """Create's Fletcher's favourite GHOSTDR directory structure for easy trouble-shooting. Does not add your data"""

    dirpath = os.path.join('./', DATA_parent,)
    tryMakePath(dirpath)
    for i in obj_name:
        objpath = os.path.join(dirpath,  i)
        tryMakePath(objpath)
        for j in readmodes:
            readmodepath = os.path.join(objpath,  j)
            tryMakePath(readmodepath)
            rawpath = os.path.join(readmodepath,  'raw')
            tryMakePath(rawpath)
            tryMakePath(os.path.join(readmodepath,  'intermediate'))
            tryMakePath(os.path.join(rawpath,  'packed'))
            tryMakePath(os.path.join(rawpath,  'obj'))    
            
def show2dCutout(wmin,wmax,path='',file='',flat='',wvl=True,pxl=True, weight=False, camera='blue',debug=False):
    orderref = {'blue':35, 'red':34}
    mrefdict = {'blue': 80, 'red': 50}
    ordermindict = {'blue': 64, 'red': 33}
    wavemodfile =  {'blue': '/Users/fwaller/Documents/GHOST/GHOSTDR/ghostdr/ghost/lookups/Polyfit/blue/std/220620/wavemod.fits',
                    'red': '/Users/fwaller/Documents/GHOST/GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/std/220620/wavemod.fits'}
    echelle = fits.open(path+file)
    xprimes = echelle[1].data.shape[0]
    orderwave=[]
    m_ref = mrefdict[camera]
    wavemodhdu = fits.open(wavemodfile[camera])
    wavemod = wavemodhdu[0].data
    yprimes = np.arange(echelle[1].data.shape[1])
    for i in range(orderref[camera]):
        m = i + ordermindict[camera]
        poly_w_model = np.poly1d([np.poly1d(wavemod[i,:])((m_ref/m)-1) for i in range(len(wavemod[:,0]))])
        wave = poly_w_model(yprimes-(yprimes[-1]+1)//2)
        orderwave.append(wave)
#     print(orderwave)
    ooi=[]
    ooilim = {'min':[], 'mid':[], 'max':[], 'num':[]}
    for i,a in enumerate(orderwave):
        if (a[0] < wmax and a[-1] > wmax): #lower portion is in range
            rang = a[(a < wmax) & (a>wmin)]
            ooi.append(rang)
            ooilim['min'].append(rang[0])
            ooilim['mid'].append(rang[len(rang)//2])
            ooilim['max'].append(rang[-1])
            ooilim['num'].append(i)
        elif a[0] > wmin and a[-1] < wmax: # entire order is in range
            rang = a
            ooi.append(rang)
            ooilim['min'].append(rang[0])
            ooilim['mid'].append(rang[len(rang)//2])
            ooilim['max'].append(rang[-1])
            ooilim['num'].append(i)
        elif a[0] < wmin and a[-1] > wmin: #upper portion is in range
            rang = a[(a > wmin) & (a<wmax)]
            ooi.append(rang)
            ooilim['min'].append(rang[0])
            ooilim['mid'].append(rang[len(rang)//2])
            ooilim['max'].append(rang[-1])
            ooilim['num'].append(i)
    if debug:
        print(ooilim)
    flatpath = '/Users/fwaller/Documents/GHOST/DATA/BPSCS31082/NORM/'
    flatfile = flat
    flat = fits.open(flatpath+flatfile)
    xmod = flat[4].data
    print(camera, 'camera')
    print('The wavelengths here are an approxmation at best!!! do not worry if they only vaguely line up.')
    if wvl:
        wfig = plt.figure(figsize=(15,len(ooi)))
    if pxl:
        pfig = plt.figure(figsize=(15,len(ooi)))
    echelle = fits.open(path+file)
    m_ref = mrefdict[camera]
    yprimes = np.arange(echelle[1].data.shape[1])
    for g,a in enumerate(ooi):
        if camera == 'red':
            j = -(g+1)
        if camera == 'blue':
            j = g
        k = ooilim['num'][j] 
        m = k + ordermindict[camera]
        poly_x_model = np.poly1d([np.poly1d(xmod[l,:])((m_ref/m)-1) for l in range(len(xmod[:,0]))])
        poly_w_model = np.poly1d([np.poly1d(wavemod[n,:])((m_ref/m)-1) for n in range(len(wavemod[:,0]))])
        wave = poly_w_model(yprimes-(yprimes[-1]+1)//2)
        xx = poly_x_model(yprimes- (yprimes[-1]+1)//2) +xprimes//2
        if camera == 'red':
            spatmin = xx[wave==ooilim['min'][j]]
            spatmax = xx[wave==ooilim['max'][j]]
        if camera == 'blue':
            spatmin = xx[wave==ooilim['max'][j]]
            spatmax = xx[wave==ooilim['min'][j]]
        ymin = yprimes[wave==ooilim['min'][j]]
        ymax = yprimes[wave==ooilim['max'][j]]
        vmin=-1
        vmax=6
        buf = 50
        if weight:
            if wvl:
                wax=wfig.add_axes([(ooilim['min'][j]-wmin)/(wmax-wmin), float(np.abs(j)/len(ooi)), (ooilim['max'][j]-ooilim['min'][j])/(wmax-wmin), 0.9/float(len(ooi))]) #left bottom width height
                wax.imshow(echelle[3].data[0,int(round(spatmin[0]))-buf:int(round(spatmax[0]))+buf,ymin[0]:ymax[0]],origin='lower',vmin=vmin,vmax=vmax,extent=[ooilim['min'][j],ooilim['max'][j],spatmin[0]-buf,spatmax[0]+buf],aspect='auto')
                wax.plot(wave[(wave > ooilim['min'][j]) & (wave < ooilim['max'][j])],xx[(wave > ooilim['min'][j]) & (wave < ooilim['max'][j])])
            if pxl:
                pax=pfig.add_axes([float(ymin/(yprimes[-1]+1)), float(np.abs(j)/len(ooi)), float((ymax-ymin)/(yprimes[-1]+1)), 0.9/float(len(ooi))]) #left bottom width height
                pax.imshow(echelle[3].data[0,int(round(spatmin[0]))-buf:int(round(spatmax[0]))+buf,ymin[0]:ymax[0]],origin='lower',vmin=vmin,vmax=vmax,extent=[ymin[0],ymax[0],spatmin[0]-buf,spatmax[0]+buf])
                pax.plot(yprimes[(yprimes<ymax) & (yprimes>ymin)],xx[(yprimes<ymax) & (yprimes>ymin)])
        else:
            if wvl:
                wax=wfig.add_axes([(ooilim['min'][j]-wmin)/(wmax-wmin), float(np.abs(j)/len(ooi)), (ooilim['max'][j]-ooilim['min'][j])/(wmax-wmin), 0.9/float(len(ooi))]) #left bottom width height
                wax.imshow(echelle[1].data[int(round(spatmin[0]))-buf:int(round(spatmax[0]))+buf,ymin[0]:ymax[0]],origin='lower',vmin=vmin,vmax=vmax,extent=[ooilim['min'][j],ooilim['max'][j],spatmin[0]-buf,spatmax[0]+buf],aspect='auto')
                wax.plot(wave[(wave > ooilim['min'][j]) & (wave < ooilim['max'][j])],xx[(wave > ooilim['min'][j]) & (wave < ooilim['max'][j])])
            if pxl:
                pax=pfig.add_axes([float(ymin/(yprimes[-1]+1)), float(np.abs(j)/len(ooi)), float((ymax-ymin)/(yprimes[-1]+1)), 0.9/float(len(ooi))]) #left bottom width height
                pax.imshow(echelle[1].data[int(round(spatmin[0]))-buf:int(round(spatmax[0]))+buf,ymin[0]:ymax[0]],origin='lower',vmin=vmin,vmax=vmax,extent=[ymin[0],ymax[0],spatmin[0]-buf,spatmax[0]+buf])
                pax.plot(yprimes[(yprimes<ymax) & (yprimes>ymin)],xx[(yprimes<ymax) & (yprimes>ymin)])
            
    if (pxl==False) & (wvl==False):
        print('Error: enable pxl=True or wvl=True')
    plt.show()
    #return
    