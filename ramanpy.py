#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
#########################################################################
File        : ramanpy.py
Author      : U. Wolfram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Raman analyses

import ramanpy as rp

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Modules:

BackgroundCorrection() ... performs quick and drity background correction based on a brutal median filtering of the original spectrum

SpectralAnalyses() ... analses of individual raman peaks. Yields fwhm, integreated area, peak intensity, peak position


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example when used with in main:


if __name__ == '__main__':

    ################################################################
    # initialising classes to extract data
    bgc = BackgroundCorrection()
    anal = SpectralAnalyses()
    
    ################################################################
    stdout.write( '\n S T A R T   %s \n\n' % "EXTRACT RAMAN DATA" )
    stdout.flush()
    startTime = clock()
    ctime1    = time()
    
    ################################################################
    # begin raman background correction
    bgc.setFileName(argv[1])
    bgc.readShifts()
    bgc.medianFilter()
    bgc.backgroundCorrection()
    bgc.writeCorrected()    

    if (len(argv) > 2 and argv[2] == 'plot'):
        bgc.plotFunction(bgc.getShifts(), bgc.getCounts(), "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum")
        bgc.plotFunction(bgc.getShifts(), bgc.getCorrected(), "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum")
    # end raman background correction
    ################################################################
    
    
    
    ################################################################
    # begin raman analyses
    anal.setPlotFilename(argv[1].split('.')[0] + '_v1PO4.pdf')
    anal.setShifts(bgc.getShifts())
    anal.setCounts(bgc.getCorrected())
    anal.setBand(1000, 1150)        
    fwhm, integ, peak, pos = anal.peakAnalysis()    
        
    # end raman analyses
    ################################################################
    
    
    endTime = clock()
    ctime2 = time()
    stdout.write( '\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime-startTime, ctime2-ctime1) )
    stdout.flush()
    ################################################################


#########################################################################
"""

from sys import *
from time import *
from math import *
import numpy as np
from string import *
import os # os must be imported like that. otherwise open(...) cannot be used. 
import re
import glob

import matplotlib.pyplot as plt
from BaselineRemoval import BaselineRemoval

from scipy import ndimage
from scipy import signal
from scipy import stats
from scipy import integrate
from scipy.optimize import leastsq # Levenberg-Marquadt Algorithm

import builtins as buin


class BackgroundCorrection():
    
    __version__='V_June_2018'
    __author__='U. Wolfram'

    def __init__(self):

        # variables identifying the specimen
        self.filename = 'foo'
        self.shifts = []
        self.counts = []
        self.smoothed = np.array([])
        self.corrected = np.array([])
        
    def getShifts(self):
        return self.shifts

    def getCounts(self):
        return self.counts

    def getSmoothed(self):
        return self.smoothed
        
    def getCorrected(self):
        return self.corrected

    def setShifts(self, mylist):
        self.shifts = mylist

    def setCounts(self, mylist):
        self.counts = mylist
        
    def setFileName(self, filename):
        self.filename = filename
        
        
    def breakPoint(self):

        stdout.write("\n ... B R E A K  P O I N T \n")
        stdout.flush()
        stop = "n"
        print("\n     Going on (y/n): ")
        stop = input()            
        if (stop.lower() != "y"):
            exit(0)
        elif (stop.lower() == "y"):
            return            
    
    
    def readShifts(self):
        
        stdout.write("\n \t... read data \n\n"); stdout.flush()
        
        isExist = os.path.isfile(self.filename)
        if isExist:
                
            skip = 0
            for line in open(self.filename,'r'):
                line = re.sub(' +','',line)
                line = line.replace(",", ".")
                line = line.replace("\n", "")
                line = line.replace("\r", "")
    
                if skip >0:
                    self.shifts.append(float(line.split('\t')[0]))
                    self.counts.append(float(line.split('\t')[1]))
                    
                skip+=1
                
        elif not isExist:
            stdout.write("\n \t    file does not exist - exiting \n\n"); stdout.flush()
            exit(0)

    def medianFilter(self):
        
        stdout.write("\n \t... filter data \n\n"); stdout.flush()
        
        kernelsize = round(len(self.counts)/3)
        
        if (np.mod(kernelsize,2) == 0): 
            kernelsize = kernelsize + 1
        
        # median filtering of signal
        self.smoothed = signal.medfilt(np.array(self.counts), kernelsize)
        #kernel = signal.gaussian(11, 1)
        #smoothed = signal.fftconvolve(counts, kernel, mode='same')
        
        
    def backgroundCorrection(self):      
        
        # filter data
        self.medianFilter()
        
        stdout.write("\n \t... correct background \n\n"); stdout.flush()
        
        self.corrected = np.array(self.counts) - self.smoothed

        
    def backgroundCorrection2(self):
        
        
        baseObj=BaselineRemoval(self.counts)
        
        self.corrected=baseObj.ZhangFit()
                    
    
    def plotFunction(self, x, y, xlab, ylab, title, vlines = None, filename = None): 
           
        stdout.write("\n \t... plot data \n\n"); stdout.flush()
                    
        if (filename==None and vlines==None):
            
            plt.plot(x,y, 'wo', markersize=4,zorder=1)
            plt.scatter(x, y, color='b', s=10, facecolors='none', edgecolors='b',zorder=3)
            plt.grid(True)
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)
            #plt.axis(axislimits)            
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.show()
                        
        elif (filename):
            
            plt.ioff()
            fig = plt.figure(figsize=(5.5, 5.5), dpi=600)
            plt.plot(x,y, 'wo', markersize=4,zorder=1)
            plt.scatter(x, y, color='b', s=10, facecolors='none', edgecolors='b',zorder=3)#plt.plot(x, y, linewidth=2.0, color='b')
            plt.grid(True)
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)            
            #plt.axis(axislimits)
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.savefig(filename, dpi = 600, format='pdf')
            plt.close(fig)
        
        if (vlines):
            plt.plot(x,y, 'wo', markersize=4,zorder=1)
            plt.scatter(x, y, color='b', s=10, facecolors='none', edgecolors='b',zorder=3)
            plt.grid(True)
            for i in range(len(vlines)):
                plt.axvline(x=vlines[i], ymin=0, ymax = 250, linewidth=2, color='k')            
            
            plt.ylabel(ylab)
            plt.xlabel(xlab)
            plt.title(title)            
            #plt.axis(axislimits)
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.show()
                
    # write tabseparated file
    def writeCorrected(self):
        
        stdout.write("\n \t... write background correct data \n\n"); stdout.flush()
        
        outfile = self.filename.split('.')[0] + '_corrected.dat'
        f = open(outfile, 'w')
        f.write('Wave\tIntensity\n')
        
        
        for k in range(len(self.counts)):

            f.write('%s\t%s\n' % (str(self.shifts[k]), str(self.corrected[k]))) 



class SpectralAnalyses():
    
    __version__='V_May_2015'
    __author__='U. Wolfram and J. Schwiedrzik'

    def __init__(self):

        # variables identifying the specimen
        self.filename = 'foo'
        self.plotFileName = 'v1PO4.pdf'
        self.shifts = []
        self.counts = []
        self.bandind = []
        
    def getShifts(self):
        return self.shifts

    def getCounts(self):
        return self.counts 
        
    def setPlotFilename(self, plotFileName):
        self.plotFileName = plotFileName
        
        
    def setShifts(self, data):
        if isinstance(data, np.ndarray):
            self.shifts = data.tolist()
        elif isinstance(data, list):
            self.shifts = data

    def setCounts(self, data):
        if isinstance(data, np.ndarray):
            self.counts = data.tolist()
        elif isinstance(data, list):
            self.counts = data

    def setBand(self, lower, upper):
        
        # getting indizees to isolate peaks
        lb = self.shifts.index(min(self.shifts, key=lambda x:abs(x - lower)))
        ub = self.shifts.index(min(self.shifts, key=lambda x:abs(x - upper)))
        
        if lb > ub:
            self.shifts.reverse()
            self.counts.reverse()
            lb = self.shifts.index(min(self.shifts, key=lambda x:abs(x - lower)))
            ub = self.shifts.index(min(self.shifts, key=lambda x:abs(x - upper)))
        
        self.bandind.append(lb)
        self.bandind.append(ub)
        
        
    def breakPoint(self):

        stdout.write("\n ... B R E A K  P O I N T \n")
        stdout.flush()
        stop = "n"
        print("\n     Going on (y/n): ")
        stop = input()            
        if (stop.lower() != "y"):
            exit(0)
        elif (stop.lower() == "y"):
            return
        
        
    def fwhm2(self, x, y, plot = False, filename = None):
       
        # initial values [hwhm, peak center, intensity]
        # see http://en.wikipedia.org/wiki/Full_width_at_half_maximum for the hwhm estimation
        
                
        ############################
        # go on here, two step procedure seems necessary but greatly exaggerates the maximum
        # the second iteration somehow does not yield a good position for the maximum using x[y.index(max(y))] yields a Gaussian with a two narrow std

        
        # defining x and y new
        xnew = x[y.index(max(y)):len(x)]
        ynew = y[y.index(max(y)):len(y)]
        
        
        #get Gaussian parameters
        pbest = leastsq(self.residuals, [2.355*np.std(x)/2.0, x[y.index(max(y))], max(y)], args=(ynew,xnew),full_output=1)        
        best_parameters = pbest[0]


        # get fit to data and integral
        fit = self.lorentzian(x, best_parameters)
        integral = integrate.quad(lambda x: self.lorentzian(x, best_parameters), min(x), max(x))
                
        if plot:
            
            plt.ioff()
            fig = plt.figure(figsize=(4.5, 4.5), dpi=600)
            plt.plot(x,y, 'wo', markersize=4,zorder=1)
            plt.scatter(x, y, color='b', s=10, facecolors='none', edgecolors='b',zorder=3)
            plt.plot(x,fit,'r-', lw=2, zorder=2)
            plt.xlabel(r'$\omega$ in cm$^{-1}$', fontsize=18)
            plt.ylabel('Intensity (a.u.)', fontsize=18)
            plt.grid(True)
            ax = plt.gca()
            ax.set_autoscale_on(False)
            plt.savefig(filename, bbox_inches="tight", dpi = 600, format='pdf')
            plt.close(fig)
            
        # recast best parameters to [fwhm, peak center, intensity]        
        para = [2*best_parameters[0], best_parameters[1], best_parameters[2]]
        return para, integral
    
    def lorentzian(self, x, p):        
        # p = [hwhm, peak center, intensity]
        numerator =  ( p[0]**2 )
        denominator = ( x - (p[1]) )**2 + p[0]**2
        y = p[2]*(numerator/denominator)
        return y    

    def gaussian (self, x, p):
        
        numerator = -1.0*(x - p[1])**2
        denominator = 2.0*p[2]**2        
        y = p[0]*np.exp(numerator/denominator)
        return y
        
    def residuals(self, p, y, x):        
        err = (y - self.lorentzian(x, p))**2
        return err
    
    def residualsGaussian(self, p, y, x):
        err = (y - self.gaussian(x, p))**2
        return err        


    def peakAnalysis(self):
                
        stdout.write("\n \t... Analyse peaks \n\n"); stdout.flush()

        shifts = self.shifts[self.bandind[0]:self.bandind[1]]
        counts = self.counts[self.bandind[0]:self.bandind[1]]
    
        fit, integral = self.fwhm2(shifts, counts, plot=True, filename=self.plotFileName)
        
        fwhm = fit[0]
        integ = integral[0]
        peak = fit[2]
        pos = fit[1]
    
        return fwhm, integ, peak, pos
        

#if __name__ == '__main__':
#
#    ################################################################
#    # initialising classes to extract data
#    bgc = BackgroundCorrection()
#    anal = SpectralAnalyses()
#    
#    ################################################################
#    stdout.write( '\n S T A R T   %s \n\n' % "EXTRACT RAMAN DATA" )
#    stdout.flush()
#    startTime = clock()
#    ctime1    = time()
#    
#    ################################################################
#    # begin raman background correction
#    bgc.setFileName(argv[1])
#    bgc.readShifts()
#    bgc.medianFilter()
#    bgc.backgroundCorrection()
#    bgc.writeCorrected()    
#
#    if (len(argv) > 2 and argv[2] == 'plot'):
#        bgc.plotFunction(bgc.getShifts(), bgc.getCounts(), "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum")
#        bgc.plotFunction(bgc.getShifts(), bgc.getCorrected(), "Raman Shifts 1/cm", "Corrected Counts", "Raman Spectrum")
#    # end raman background correction
#    ################################################################
#    
#    
#    
#    ################################################################
#    # begin raman analyses
#    anal.setPlotFilename(argv[1].split('.')[0] + '_v1PO4.pdf')
#    anal.setShifts(bgc.getShifts())
#    anal.setCounts(bgc.getCorrected())
#    anal.setBand(1000, 1150)        
#    fwhm, integ, peak, pos = anal.peakAnalysis()    
#        
#    # end raman analyses
#    ################################################################
#    
#    
#    endTime = clock()
#    ctime2 = time()
#    stdout.write( '\n E N D E D  SUCCESSFULLY in  CPU/TOT : %8.1f/%8.1f sec \n\n' % (endTime-startTime, ctime2-ctime1) )
#    stdout.flush()
#    ################################################################
