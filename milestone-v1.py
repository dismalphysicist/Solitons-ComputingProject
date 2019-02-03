# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 10:15:11 2018

@author: Lois
"""

import numpy


def fourier_transform(f):
    F = numpy.fft.fft(f)
    return F

def inverse_fourier(F):
    f = numpy.fft.ifft(F)
    return f
 
def gaussian_analytical():
    #solve S.E. analytically for Gaussian in diffracting potential 
    wf = numpy.array(x,t,psi) #array of wavefunctions at each time and space point
    return wf
  	
def schrodinger_time_evolution(psi_initial, V, nonlinear):
    psi = psi_initial
    psis = numpy.zeros([len(psi_initial),len(t)]) #initialising 2d array
	
#    for all time points t + dt: 	#try using enumerate to get both index over x,t and x,t themselves
#        array[x,t] = psi
#
#        if nonlinear == False:
#            psi = exp(-iVdt/2) inverse_fourier(exp(-ik^2dt) fourier_transform(exp(-iVdt/2) psi))
#        else: 
#            psi = 0 #fill in later 
#        return 2D array of psis

#test
print fourier_transform(numpy.array([1,2,1]))
print inverse_fourier(numpy.array([1,2,1]))
print schrodinger_time_evolution(numpy.array([1,2,1]),0,False)

# #to compare to analytic, use gaussian first. later, sech soliton and put density proportional term in V.
# 
# psi_initial = array[x,psi] #gaussian distribution along x
# 
# psi_gaussian = schrodinger_time_evolution(psi_initial, V = 0)
# 
# psi_harmonic_oscillator(psi_initial, V = 1/2 kx^2) #parabolic confining potential added
# 
# 
# pyplot.subplots 	#3 so can compare to analytic? or just 2 and open two separate windows? 
# 
# pyplot.imshow(mod(psi_whichever)^2, extent = (x extremes, t extremes)) for each subplot psi_gaussian, psi_harmonic_oscillator and maybe analytic gaussian
# colorbar(mod(psi_whichever)^2, (0,1), black and white) #show probability as black/white gradient from 0 to 1
# labels etc
# show()
#==============================================================================
