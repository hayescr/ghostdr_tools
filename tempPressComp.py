
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from math import log, floor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import glob


#Input path to arc files you wish to compare 
arcList = []
arcList.append(glob.glob('/home/matt/GHOSTData/BPSCS31082/calibrations/processed_arc/*1_arc.fits'))
arcList.append(glob.glob('/home/matt/GHOSTData/S0002/calibrations/processed_arc/*1_arc.fits'))
arcList.append(glob.glob('/home/matt/GHOSTData/S0054/calibrations/processed_arc/*1_arc.fits'))
arcList.append(glob.glob('/home/matt/GHOSTData/HIP016085/calibrations/processed_arc/*_arc.fits'))

arcName = "" ##INPUT OBJECT NAME
#Output files 
pp1 = PdfPages(f'TempPresComparisonRed.pdf')
pp2 = PdfPages(f'TempPresComparisonBlue.pdf')
pp = PdfPages(f'TempPresComparisonFull.pdf')

#Initiating lists
wavemodR = []
tempR = []
pressureR = []
wavemodB = []
tempB = []
pressureB = []
dir = []


#Taking in data for temperature, pressure, and wavemod value for each coefficient
for arc in arcList:
    for colour in arc:
        if fits.open(colour)[0].header['CAMERA']=='BLUE':
            for i in range(6):
                for j in range(6):
                    dir.append((fits.open(colour)[4].data[5-i,5-j]))
            wavemodB.append(dir)            
            tempB.append(fits.open(colour)[0].header['GRATMNTT'])
            pressureB.append(fits.open(colour)[0].header['RPIP'])
            dir = []
        else:
            for i in range(6):
                for j in range(6):
                    dir.append((fits.open(colour)[4].data[5-i,5-j]))
            wavemodR.append(dir)
            tempR.append(fits.open(colour)[0].header['GRATMNTT'])
            pressureR.append(fits.open(colour)[0].header['RPIP'])
            dir = []


wavemodArrR = np.array(wavemodR)
wavemodArrB = np.array(wavemodB)     

#Plotting 
for colour in range(2):
    fig, axs = plt.subplots(6, 6)
    for i in range(6):
        for j in range(6):
            if colour % 2 == 0:
                axs[5-i, 5-j].scatter(tempR, wavemodArrR[:,6*i+j], c='r')
                axs[5-i, 5-j].set_yticklabels([])
                axs[5-i, 5-j].set_xticklabels([])

                fig1 = plt.figure()
                ax = fig1.add_subplot(111)
                ax.set_title(f'Q{i}{j} Coefficinet vs Temperature (Red)')
                ax.set_xlabel('Temp (C)')
                ax.set_ylabel('Wavemod Coefficient')
                plot = ax.scatter(tempR, wavemodArrR[:,6*i+j], c=pressureR, cmap='jet')
                cb = plt.colorbar(plot, label='Pressure (hPA)')
                pp1.savefig()
                plt.close(fig1)
                cb.remove()
            else:
                axs[5-i, 5-j].scatter(tempB, wavemodArrB[:,6*i+j], c='b')
                axs[5-i, 5-j].set_yticklabels([])
                axs[5-i, 5-j].set_xticklabels([])

                fig2 = plt.figure()
                ax = fig2.add_subplot(111)
                ax.set_title(f'Q{i}{j} Coefficinet vs Temperature (Blue)')
                ax.set_xlabel('Temp (C)')
                ax.set_ylabel('Wavemod Coefficient')
                plot = ax.scatter(tempB, wavemodArrB[:,6*i+j], c=pressureB, cmap='jet')
                cb = plt.colorbar(plot, label='Pressure (hPA)')
                pp2.savefig()
                plt.close(fig2)

    pp.savefig()
    plt.close(fig)

plt.close(fig1)
plt.close(fig2)

pp.close()
pp1.close()
pp2.close()