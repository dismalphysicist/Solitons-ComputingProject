"""
Created on Sun Oct 28 10:15:11 2018
Version 5 
@author: Lois
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot


def fourier_transform(f):
    F = np.fft.fft(f)
    return F

def inverse_fourier(F):
    f = np.fft.ifft(F)
    return f
	
def getks(xs, dx):
	ks = np.fft.fftfreq(len(xs), dx) 
	return ks
	
def schrodinger_time_evolution(psi_initial, V, g):
	
	ks = getks(xs,dx) #accessing global variables 
	
	psisquareds = np.zeros((len(psi_initial),n_times)) #initialising 2d array
	
	psi = psi_initial #prepare for loop
	psisq = abs(psi_initial)**2 
	
	for it in range(n_times):
		psisquareds[:,it] = psisq #fill array 
		
		interactions = g * abs(psi)**2
			
		step1 = fourier_transform(np.exp(-1j*(V + interactions)*dt/2) * psi)
		step2 = inverse_fourier(np.exp(-1j*ks**2*dt) * step1)
		psi = np.exp(-1j*(V + interactions)*dt/2) * step2 
		
		psisq = abs(psi)**2

		
	return psisquareds

#global variables for extent of modelling 

n_times = 2000 #number of timesteps (also number of x steps)
dt = 0.01
xs = np.linspace(-10, 10, n_times)
L = xs[-1] - xs[0] #box length 
dx = L/n_times


#generate initial psi
def gaussian():
	psi_initial = np.exp(-xs**2)
	return psi_initial

def soliton():
	family_param = 1
	velocity = L/4 #or L/6
	
	#psi_initial = 2**0.5 * np.exp(1j*L*xs/4) / np.cosh(xs)
	psi_initial = (2*family_param)**0.5 * np.exp(1j*velocity*xs) / np.cosh(xs * (family_param)**0.5)
	return psi_initial

def accuracy_test(g_test,t):
    psisquareds_test = schrodinger_time_evolution(soliton(), 0, g_test)
    t = t/dt #taking timestep into account 
    norm = np.sum(psisquareds_test[750:1250],0) # sum over axis 0, i.e. space between -2 and 2
    accuracy = (norm[0] - norm[t])*100/norm[0] #percentage change in norm 
    return accuracy

#tests

testpsi = schrodinger_time_evolution(soliton(), 0, 0) #no nonlinear term 
g_test = -1/(2*L)
g_test2 = -0.05
testpsi2 = schrodinger_time_evolution(soliton(), 0, g_test) #nonlinear interactions 

#testing accuracy
acc = accuracy_test(g_test, 19.99)
print "Norm is conserved to {} % accuracy".format(acc) 

#to compare to analytic, use gaussian first. later, sech soliton and nonlinear. 

pyplot.figure(figsize=(10,6)) 
pyplot.suptitle("The effect of inter-atom interactions " + r'$g |\psi|^2$' + " on a soliton model", fontsize=16)	

pyplot.subplot(121)
pyplot.imshow(np.transpose(testpsi), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title(r'$g=0$')

pyplot.subplot(122)
pyplot.imshow(np.transpose(testpsi2), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title(r'$g= -0.025$')

pyplot.savefig('milestonepic.png')
pyplot.show()
# colorbar(mod(psi_whichever)^2, (0,1), black and white) #show probability as black/white gradient from 0 to 1

