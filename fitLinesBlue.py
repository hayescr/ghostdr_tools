#This program will allow you to check the fit of the lines for a particular BLUE arc file
#Running it will produce an excel sheet contianing the fitting parameters for all lines
#It will also indicate lines that may be causing issues and have been flagged (NOT DEFINITIVE)
#A pdf file showing the fit of each line is also made

#The 4 reasons a line may be flagged:
    #1 - There is not enough data around the lines expected posiiton to fit a curve
    #2 - The data is too noisy and the curve_fit function cannot fit a gaussian to it
    #3 - The FWHM of the gaussian is outside the range of the average
    #4 - The fit of the gaussian is outside the range of the average 

from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from math import log, floor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


def findRange(line, waveList, data):
    inRange = []
    for wavelength in waveList:
        if wavelength >= line-0.18 and wavelength <= line+0.18:
            #print(wavelength, data[np.where(waveList == wavelength)][0], np.where(waveList == wavelength))
            inRange.append([wavelength, data[np.where(waveList == wavelength)][0]])

    if len(inRange) != 0:
        return inRange
    else:
        return None

def gauss(x, a, b, sigma):
    return a*np.exp(-(x-b)**2/(2*sigma**2))

def mase(actual : np.ndarray, predicted : np.ndarray):

    forecast_error = np.mean(np.abs(actual - predicted))
    naive_forecast = np.mean(np.abs(np.diff(actual)))
    mase = forecast_error / naive_forecast

    return mase
    
#! FILE INPUT REQUIRED !#
baseData = fits.open('Path to /GHOSTDR/ghostdr/ghost/lookups/Polyfit/red/high/220620/wavemod.fits')
flatData = fits.open('Path to flat.fits file')
arcData = fits.open('Path to arc.fits file')
thardata_full = np.genfromtxt('Path to /GHOSTDR/ghostdr/ghost/lookups/Polyfit/ghost_thar_linelist_20220718.txt', usecols=(1,), dtype=None, encoding=None)
arcName = "OBJ" ##INPUT OBJECT NAME

baseWavemod = baseData[0].data
arcWavemod = arcData[4].data
xmod = flatData[4].data

pp1 = PdfPages(f'curveFitting_{arcName}_blue.pdf')

wavelengths = np.arange(3475, 5480, 20)
m_ref = 80
yprimes = np.arange(arcData[1].data.shape[1])


# Initializing lists 
plotData = []
flagged = []
print("Fitting lines...")

for index, wavelength in enumerate(wavelengths):
    plot_id = (index % 5) + 1
    presentLines = thardata_full[(thardata_full < wavelength + 10) & (thardata_full > wavelength-10)]  

    for i in range(arcData[1].data.shape[0]):
        m = i + 64
        poly_arc_model = np.poly1d([np.poly1d(arcWavemod[i,:])((m_ref/m)-1) for i in range(len(arcWavemod[:,0]))])
        poly_xmod_model = np.poly1d([np.poly1d(xmod[i,:])((m_ref/m)-1) for i in range(len(xmod[:,0]))])

        wave = poly_arc_model(yprimes- 4096//2)
        xval = poly_xmod_model(yprimes- 4096//2) + 4112//2

        data = arcData[1].data[i,:,0]
    
        for line in presentLines:
            arcMaxFlux = 0
            arcMaxWave = 0
            arcRange = np.array(findRange(line, wave, data))

            if (arcRange == None).any() or len(arcRange[:,0]) < 5:
                break
            else:  
                arcWaveData, arcFluxData = arcRange[:,0], arcRange[:,1]

                #Determining initial guess for line position 
                for num in arcWaveData:
                    arcTempFlux = float(data[np.where(wave==num)])
                    if arcMaxFlux < arcTempFlux or arcMaxFlux == 0:
                        arcMaxFlux = arcTempFlux
                        arcMaxWave = num   
                
                #Fitting gaussian to data
                yData = []
                xData = []
                for i in arcWaveData:
                    yData.append((yprimes[np.where(wave==i)])[0])
                    xData.append((xval[np.where(wave==i)])[0])
                    
                try:
                    yprimeParams, p_covy = curve_fit(gauss, yData, arcFluxData, p0=((arcMaxFlux, (yprimes[np.where(wave==arcMaxWave)])[0], 0.3)), maxfev=10000)
                    xvalParams, p_covx = curve_fit(gauss, xData, arcFluxData, p0=((arcMaxFlux, (xval[np.where(wave==arcMaxWave)])[0], 0.03)), maxfev=10000)
                    waveParams, p_covw = curve_fit(gauss, arcWaveData, arcFluxData, p0=((arcMaxFlux, arcMaxWave, 0.03)), maxfev=10000)

                except RuntimeError:
                    flagged.append([line, "Line could not be fit", None, None, None, None, "Flag 2"])
                    print(f"{line} has been flagged! (Could not be fit)")

                else:      
                    #Determine width of lines in wavelength and pixels
                    fwhm = abs(waveParams[2]*np.sqrt(8*log(2)))
                    fwhmPixel = abs(yprimeParams[2]*np.sqrt(8*log(2)))
                    
                    #Determine error in the fit 
                    yFitError = gauss(yData, yprimeParams[0], yprimeParams[1], yprimeParams[2])
                    arcfitError = mase(arcFluxData, yFitError)
                    
                    #Plotting the best fit compared to data (for visualizing fit errors) 
                    ym = gauss(np.arange(yData[0], yData[len(yData)-1], 0.005), yprimeParams[0], yprimeParams[1], yprimeParams[2])
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(np.arange(yData[0], yData[len(yData)-1], 0.005), ym, c='r', label='Best Fit')                       
                    ax.scatter(yData, arcFluxData, c='k', label='Arc Data')
                    ax.plot([], [], ' ', label=f"A1 Fit Error: {round(arcfitError,2)}")                        
                    ax.plot([], [], ' ', label=f"fwhm: {round(fwhm,2)}")                        
                    ax.plot([], [], ' ', label=f"Pos: ({floor(xvalParams[1])}, {floor(yprimeParams[1])})")        
                    ax.set_title(f"{line} Line Fitting")
                    ax.legend()
                    pp1.savefig()
                    plt.close(fig)

                    plotData.append([line, yprimeParams[1], xvalParams[1], fwhm, arcfitError])
            
pp1.close()

plotArr = np.array(plotData)
avgFitError = np.mean(plotArr[:,4])
fitErrorSTD = np.std(plotArr[:,4])

#Finding fwhm for lines that area fit well
findFWHM = []
for i in range(plotArr.shape[0]):
    if plotArr[i,4] < 0.1:
        findFWHM.append(plotArr[i,3])
avgFWHM = np.mean(findFWHM)
fwhmSTD = np.std(findFWHM)

#Filtering poor lines / noisy data
print("Checking fitted lines...")
filteredData = []
for i in range(plotArr.shape[0]):
    #Checks if fit error is less that 0.1 for both arcs
    if plotArr[i,4] < 0.1:
        filteredData.append([plotArr[i,0], plotArr[i,1], plotArr[i,2], plotArr[i,3], plotArr[i,4], plotArr[i,0]/plotArr[i,3]])
    else:
        #Checks if the fwhm is around the average
        if avgFWHM-(1-fwhmSTD/avgFWHM)/20 < plotArr[i,3] < avgFWHM+(1-fwhmSTD/avgFWHM)/8:
            #Checks if the fit error is around the average
            
            if plotArr[i,4] <= avgFitError+(fitErrorSTD/avgFitError)/2:    
                filteredData.append([plotArr[i,0], plotArr[i,1], plotArr[i,2], plotArr[i,3], plotArr[i,4], plotArr[i,0]/plotArr[i,3]])
            else:
                flagged.append([plotArr[i,0], plotArr[i,1], plotArr[i,2], plotArr[i,3], plotArr[i,4], plotArr[i,0]/plotArr[i,3], "Flag 2"])
                print(f"{plotArr[i,0]} was flagged! (FLG4)")
        else:
            flagged.append([plotArr[i,0], plotArr[i,1], plotArr[i,2], plotArr[i,3], plotArr[i,4], plotArr[i,0]/plotArr[i,3], "Flag 1"])
            print(f"{plotArr[i,0]} was flagged! (FLG3)")
            
plotArr = np.array(filteredData)
flagArr = np.array(flagged)

#Writing line parameters and wavemod coefficients to excel
df1 = pd.DataFrame(plotArr, columns=['Wavelenth', 'y position', 'x position', 'fwhm', 'arcfitError', 'resolution'])
df2 = pd.DataFrame(flagArr, columns=['Wavelenth', 'y position', 'x position', 'fwhm', 'arcfitError', 'resolution', 'Cause'])
df3 = pd.DataFrame(arcWavemod)
with pd.ExcelWriter(f'comparePlotChip_{arcName}_blue.xlsx') as writer:
    df1.to_excel(writer, sheet_name='line parameters')
    df2.to_excel(writer, sheet_name='flagged lines')
    df3.to_excel(writer, sheet_name='wavemod Coefficients')