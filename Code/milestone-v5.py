"""
Created on Sun Oct 28 10:15:11 2018

@author: Lois
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.cm
import matplotlib.colors 


def fourier_transform(f):
    F = np.fft.fft(f)
    return F

def inverse_fourier(F):
    f = np.fft.ifft(F)
    return f
	
def getks(xs, dx):
	ks = np.fft.fftfreq(len(xs), dx) 
	return ks
	
def schrodinger_time_evolution(psi_initial, V, nonlinear):
	
	ks = getks(xs,dx)
	
	psisquareds = np.zeros((len(psi_initial),n_times)) #initialising 2d array
	
	psi = psi_initial #prepare for loop
	psisq = abs(psi_initial)**2 
	
	for it in range(len(ts)):
		psisquareds[:,it] = psisq #fill array 
		
		if nonlinear == False: 
			kspace = fourier_transform(np.exp(-1j*V*dt/2) * psi)
			psi = np.exp(-1j*V*dt/2) * inverse_fourier(np.exp(-1j*ks**2*dt) * kspace) 
		else:    
			psi = 0 #fill in later 
		
		psisq = abs(psi)**2
		print psisq[101] #at x=0, debugging 
		
	return psisquareds

#global variables for extent of modelling 

n_times = 200 #number of timesteps 
dt = 0.5
ts = np.linspace(0, n_times*dt, n_times) 
dx = 0.1
xs = np.linspace(-10, 10, n_times)


#generate initial psi
def gaussian():
	psi_initial = np.exp(-xs**2) 
	return psi_initial


#tests

testpsi = schrodinger_time_evolution(gaussian(), 0, False)
testpsi2 = schrodinger_time_evolution(gaussian(), xs**2, False) #harmonic potential 

# #to compare to analytic, use gaussian first. later, sech soliton and nonlinear. 

pyplot.figure(figsize=(10,6)) 	

pyplot.subplot(121)
pyplot.imshow(np.transpose(testpsi), extent=(-10,10,0,20), origin='lower') 
pyplot.xlabel("x")
pyplot.ylabel("t")

pyplot.subplot(122)
pyplot.imshow(np.transpose(testpsi2), extent=(-10,10,0,20), origin='lower') 
pyplot.xlabel("x")
pyplot.ylabel("t")


pyplot.show()
# colorbar(mod(psi_whichever)^2, (0,1), black and white) #show probability as black/white gradient from 0 to 1
