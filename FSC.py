# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 15:33:00 2015

FOURIER SHELL CORRELATION
FSC3D(image1, image2, SNRt, ring_thick)

Computes the Fourier Shell Correlation between image1 and image2, and computes
the threshold funcion T of 1 or 1/2 bit.

"""
from __future__ import division, print_function
import numpy as np
from numpy import meshgrid
import matplotlib.pyplot as plt

__all__ = ['FourierShellCorr', 'FSCPlot', 'HannApod']

def printv(txt):
    #printvtxt)
    return

def _radtap(X,Y,tappix,zerorad):
    """
    Creates a central cosine tapering.
    It receives the X and Y coordinates, tappix is the extent of
    tapering, zerorad is the radius with no data (zeros).
    """
    tau = 2*tappix # period of cosine function (only half a period is used)

    R = np.sqrt(X**2+Y**2)
    taperfunc = 0.5*(1+np.cos(2*np.pi*(R-zerorad-tau/2.)/tau))
    taperfunc = (R>zerorad+tau/2.)*1.0 + taperfunc*(R<=zerorad+tau/2)
    taperfunc = taperfunc*(R>=zerorad)
    return taperfunc

class HannApod:
    def __init__(self,outputdim,filterdim,unmodsize):
        printv('Calling the class HannApod')
        self.outputdim = outputdim
        self.unmodsize = unmodsize
        self.filterdim = filterdim

    def fract_hanning(self):#outputdim,unmodsize):
        """
        fract_hanning(outputdim,unmodsize)
        out = Square array containing a fractional separable Hanning window with
        DC in upper left corner.
        outputdim = size of the output array
        unmodsize = Size of the central array containing no modulation.
        Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array
        will contain an array of ones in the center and cosine modulation on the
        edges, the array of ones will have DC in upper left corner.
        """
        printv('Calling fract_hanning')
        N = np.arange(0,self.outputdim)
        Nc,Nr = np.meshgrid(N,N)
        if self.unmodsize == 0:
            out = (1.+np.cos(2*np.pi*Nc/self.outputdim))*(1.+np.cos(2*np.pi*Nr/self.outputdim))/4.
        else:
            #columns modulation
            outc = (1.+np.cos(2*np.pi*(Nc-np.floor((self.unmodsize-1)/2))/(self.outputdim+1-self.unmodsize)))/2.
            if np.floor((self.unmodsize-1)/2.)>0:
                outc[:,:int(np.floor((self.unmodsize-1)/2.))]=1
            outc[:,int(np.floor((self.unmodsize-1)/2)+self.outputdim+3-self.unmodsize):len(N)] = 1
            #row modulation
            outr = (1.+np.cos(2*np.pi*(Nr-np.floor((self.unmodsize-1)/2))/(self.outputdim+1-self.unmodsize)))/2.
            if np.floor((self.unmodsize-1)/2.)>0:
                outr[:int(np.floor((self.unmodsize-1)/2.)),:]=1
            outr[int(np.floor((self.unmodsize-1)/2)+self.outputdim+3-self.unmodsize):len(N),:] = 1

            out=outc*outr

        return out

    def fract_hanning_pad(self):#outputdim,filterdim,unmodsize):#(N,N,np.round(N*(1-filtertomo))):
        """
        fract_hanning_pad(outputdim,filterdim,unmodsize)
        out = Square array containing a fractional separable Hanning window with
        DC in upper left corner.
        outputdim = size of the output array
        filterdim = size of filter (it will zero pad if filterdim<outputdim)
        unmodsize = Size of the central array containing no modulation.
        Creates a square hanning window if unmodsize = 0 (or ommited), otherwise the output array
        will contain an array of ones in the center and cosine modulation on the
        edges, the array of ones will have DC in upper left corner.
        """
        if self.outputdim < self.unmodsize:
            raise SystemExit('Output dimension must be smaller or equal to size of unmodulated window')
        if self.outputdim < self.filterdim:
            raise SystemExit('Filter cannot be larger than output size')
        if self.unmodsize<0:
            self.unmodsize = 0
            printv('Specified unmodsize<0, setting unmodsize = 0')
        printv('Calling fract_hanning_pad')
        out = np.zeros((self.outputdim,self.outputdim))
        auxindini = int(np.round(self.outputdim/2.-self.filterdim/2.))
        auxindend = int(np.round(self.outputdim/2.+self.filterdim/2.))
        hanning_window = self.fract_hanning()
        out[auxindini:auxindend, auxindini:auxindend]=np.fft.fftshift(hanning_window)
        #out[auxindini:auxindend, auxindini:auxindend]=np.fft.fftshift(self.fract_hanning(filterdim,unmodsize))
        #return np.fft.fftshift(out)
        return out

class FourierShellCorr(HannApod):
    """
    FOURIER SHELL CORRELATION
    FSC3D(image1, image2, SNRt, ring_thick)

    Computes the Fourier Shell Correlation between image1 and image2, and computes
    the threshold funcion T of 1 or 1/2 bit.
    It can handle non-cube arrays, but it assumes that the voxel is isotropic.
    It applies a Hanning window of the size of the data to the data before the
    Fourier transform calculations to attenuate the border effects.

    INPUTS:
        image1 = image 1
        image2 = image 2
        SNRt = power SNR for threshold computation. Options:
                SNRt = 0.5 -> 1 bit threshold for average
                SNRt = 0.2071 -> 1/2 bit threshold for average
        ring_thick = thickness of the frequency rings.
                    Normally the pixels get assined to the closest integer pixel ring
                    in Fourier Domain. With ring_thick, each ring gets more pixels and
                    more statistics.
    Reference: M. van Heel, M. Schatzb, "Fourier shell correlation threshold
    criteria," Journal of Structural Biology 151, 250-262 (2005)

    @author: Julio Cesar da Silva (jdasilva@esrf.fr)
    default rad_apod 60 axial_apod 20
    """
    def __init__(self,img1,img2,snrt=0.2071,ring_thick=0,rad_apod=300,axial_apod=100):
        printv('Calling the class FourierShellCorr')
        self.snrt = snrt
        self.ring_thick = ring_thick
        self.img1 = np.array(img1)
        self.img2 = np.array(img2)
        self.rad_apod = rad_apod
        self.axial_apod = axial_apod
        printv('Input images have {} dimensions'.format(self.img1.ndim))
        if self.img1.shape != self.img2.shape:
            printv("Images must have the same size")
            raise SystemExit
        if ring_thick !=0:
            printv('Using ring_thick = {}'.format(ring_thick))
        printv('Using SNRt = %g' %snrt)

    def nyquist(self):
        """
        Evaluate the Nyquist Frequency
        """
        nmax = np.max(self.img1.shape)
        fnyquist = np.floor(nmax/2.0)
        f = np.arange(0,fnyquist+1)
        return f, fnyquist

    def ringthickness(self):
        """
        Define ring_thick
        """
        n = self.img1.shape
        nmax = np.max(n)
        x = np.arange(-np.fix(n[1]/2.0),np.ceil(n[1]/2.0))*np.floor(nmax/2.0)/np.floor(n[1]/2.0)
        y = np.arange(-np.fix(n[0]/2.0),np.ceil(n[0]/2.0))*np.floor(nmax/2.0)/np.floor(n[0]/2.0)
        if self.img1.ndim==3:
            z = np.arange(-np.fix(n[2]/2.0),np.ceil(n[2]/2.0))*np.floor(nmax/2.0)/np.floor(n[2]/2.0)
            X = meshgrid(x,y,z)
        elif self.img1.ndim==2:
            X = np.meshgrid(x,y)
        else:
            printv('Number of dimensions is different from 2 or 3.Exiting...')
            raise SystemExit('Number of dimensions is different from 2 or 3.Exiting...')
        sumsquares = np.zeros_like(X[-1])
        for ii in np.arange(0,self.img1.ndim):
            sumsquares += X[ii]**2
        index = np.round(np.sqrt(sumsquares))
        return index

    def apodization(self):
        """
        Compute the Hanning window of the size of the data for the apodization
        """
        n = self.img1.shape
        if self.img1.ndim==2:
            window = np.outer(np.hanning(n[0]),np.hanning(n[1]))
        elif self.img1.ndim==3:
            window1 = np.hanning(n[0])
            window2 = np.hanning(n[1])
            window3 = np.hanning(n[2])
            windowaxial = np.outer(window2,window3)
            windowsag = np.array([window1 for ii in range(n[1])]).swapaxes(0,1)
            #win2d = np.rollaxis(np.array([np.tile(windowaxial,(1,1)) for ii in range(n[0])]),1).swapaxes(1,2)
            win2d = np.array([np.tile(windowaxial,(1,1)) for ii in range(n[0])])
            window = np.array([np.squeeze(win2d[:,:,ii])*windowsag for ii in range(n[2])]).swapaxes(0,1).swapaxes(1,2)
        else:
            printv('Number of dimensions is different from 2 or 3. Exiting...')
            raise SystemExit('Number of dimensions is different from 2 or 3. Exiting...')
        return window

    def circle(self):
        printv('Calculating the axial apodization')
        if self.img1.ndim ==2:
            shape_x = self.img1.shape[1]
            shape_y = self.img1.shape[0]
        elif self.img1.ndim ==3:
            shape_x = self.img1.shape[2]
            shape_y = self.img1.shape[1]
        x_array = np.arange(0,shape_x)
        y_array = np.arange(0,shape_y)
        self.X,self.Y = np.meshgrid(x_array-np.round(shape_x/2.),y_array-np.round(shape_y/2.))
        circular_region=1-_radtap(self.X,self.Y,self.rad_apod,np.round(shape_x/2.)-self.rad_apod)
        return circular_region

    def transverse_apodization(self):
        """
        Compute the Hanning window of the size of the data for the apodization
        """
        printv('Calculating the transverse apodization')
        n = self.img1.shape
        HannApod.__init__(self,n[0],n[0],n[0]-2*self.axial_apod)
        filters=HannApod.fract_hanning_pad(self)
        window1d=filters[:,int(filters.shape[0]/2)]
        window2d=np.array([window1d for ii in range(n[1])]).swapaxes(0,1)
        return window2d

    def fouriercorr(self):
        """
        Compute FSC and threshold
        """
        # Apodization
        n = self.img1.shape
        circular_region = self.circle()
        if self.img1.ndim ==2:
            self.window = circular_region
        elif self.img1.ndim==3:
            window2D = self.transverse_apodization()
            circle3D = np.asarray([circular_region for ii in range(n[0])])
            self.window = np.array([np.squeeze(circle3D[:,:,ii])*window2D for ii in range(n[2])]).swapaxes(0,1).swapaxes(1,2)
            printv('Apodization in 3D')

        # FSC computation
        F1 = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(self.img1*self.window)))
        #F1 = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(self.img1)))
        #printvF1.shape)
        F2 = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(self.img2*self.window)))
        #F2 = np.fft.ifftshift(np.fft.fftn(np.fft.fftshift(self.img2)))

        C,C1,C2,npts = [[],[],[],[]]
        printv('Calling method fouriercorr from the class FourierShellCorr')
        index = self.ringthickness()
        f,fnyquist = self.nyquist()
        for ii in f:
            if self.ring_thick ==0:
                auxF1 = F1[np.where(index==ii)]
                auxF2 = F2[np.where(index==ii)]
            else:
                auxF1 = F1[(np.where( (index>=(ii-self.ring_thick//2)) & (index<=(ii+self.ring_thick//2)) ))]
                auxF2 = F2[(np.where( (index>=(ii-self.ring_thick//2)) & (index<=(ii+self.ring_thick//2)) ))]
            C.append(np.sum(auxF1*np.conj(auxF2)))
            C1.append(np.sum(auxF1*np.conj(auxF1)))
            C2.append(np.sum(auxF2*np.conj(auxF2)))
            npts.append(auxF1.shape[0])
            # The correlation
        FSC = np.abs(np.asarray(C))/(np.sqrt(np.asarray(C1)*np.asarray(C2)))

        npts = np.asarray(npts)
        # Threshold computation
        Tnum = (self.snrt + (2*np.sqrt(self.snrt)/np.sqrt(npts+np.spacing(1)))+1/np.sqrt(npts))
        Tden = (self.snrt + (2*np.sqrt(self.snrt)/np.sqrt(npts+np.spacing(1)))+1)
        T= Tnum/Tden

        return FSC, T

class FSCPlot(FourierShellCorr):
    """
    Plot the FSC and threshold curves
    """
    def __init__(self,img1,img2,snrt=0.2071,ring_thick=0,rad_apod=300,axial_apod=100):
        printv('calling the class FSCplot')
        FourierShellCorr.__init__(self, img1, img2, snrt, ring_thick,rad_apod,axial_apod)
        self.FSC, self.T = FourierShellCorr.fouriercorr(self)
        self.f, self.fnyquist = FourierShellCorr.nyquist(self)
    def plot(self):
        printv('calling method plot from the class FSCplot')
        plt.close()
        plt.figure("FSC %s %s %s"%(self.ring_thick, self.rad_apod, self.axial_apod))
        plt.clf()
        plt.plot(self.f/self.fnyquist,self.FSC.real,'-b', label='FSC')
        plt.legend()
        i =self.get_intersect()
        plt.plot([i,i],[0,1],label="%s"%i)
        if self.snrt == 0.2071:
            plt.plot(self.f/self.fnyquist, self.T, '--r',label='1/2 bit threshold')
            plt.legend()
        elif self.snrt == 0.5:
            plt.plot(self.f/self.fnyquist, self.T, '--r',label='1 bit threshold')
            plt.legend()
        else:
            plotT = plt.plot(self.f/self.fnyquist, self.T)
            plt.legend(plotT,'Threshold SNR = %g ' %self.snrt, loc='center')
        plt.xlim(0,1)
        plt.ylim(0,1.1)
        plt.xlabel('Spatial frequency/Nyquist')
        plt.ylabel('Magnitude')
        #plt.show()
        if self.img1.ndim==2:
            plt.savefig('FSC_2D.png', bbox_inches='tight')
        elif self.img1.ndim==2:
            plt.savefig('FSC_3D_%s_%s_%s.png'%(self.ring_thick, self.rad_apod, self.axial_apod), bbox_inches='tight')

    def plot_curfig(self,name=""):
        printv('calling method plot from the class FSCplot')
        if name == "":
            name='FSC %s %s %s'%(self.ring_thick, self.rad_apod, self.axial_apod)
        plt.plot(self.f/self.fnyquist,self.FSC.real, label=name)
        plt.legend()
    def get_intersect(self):
        ro = np.argmax(self.FSC.real < self.T)
        dr = self.FSC.real - self.T
        frac = dr[ro-1]/(dr[ro-1] - dr[ro])
        df = (self.f/self.fnyquist)[ro] - (self.f/self.fnyquist)[ro-1]
        return (self.f/self.fnyquist)[ro-1] + frac * df

    def plot_nyquist(self):
        if self.snrt == 0.2071:
            plt.plot(self.f/self.fnyquist, self.T, '--r',label='1/2 bit threshold')
            plt.legend()
        elif self.snrt == 0.5:
            plt.plot(self.f/self.fnyquist, self.T, '--r',label='1 bit threshold')
            plt.legend()
        else:
            plotT = plt.plot(self.f/self.fnyquist, self.T)
            plt.legend(plotT,'Threshold SNR = %g ' %self.snrt, loc='center')
        plt.xlim(0,1)
        plt.ylim(0,1.1)
        plt.xlabel('Spatial frequency/Nyquist')
        plt.ylabel('Magnitude')
    def save_fig(self,prefix):
        if self.img1.ndim==2:
            plt.savefig('FSC_2D.png', bbox_inches='tight')
        elif self.img1.ndim==3:
            plt.savefig('%sFSC_3D_%s_%s_%s.png'%(prefix,self.ring_thick, self.rad_apod, self.axial_apod), bbox_inches='tight')
