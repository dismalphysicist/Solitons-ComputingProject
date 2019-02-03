"""
Created on Sun Oct 28 10:15:11 2018
Version 5 
@author: Lois
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot

############ FUNCTIONS #############

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
	
	ks = 2*np.pi*getks(xs,dx) #accessing global variables. Factor of 2pi
	
	psisquareds = np.zeros((len(psi_initial),n_points)) #initialising 2d array
	
	psi = psi_initial #prepare for loop
	psisq = abs(psi_initial)**2 
	
	for it in range(n_points):
		psisquareds[:,it] = psisq #fill array 
		
		interactions = g * abs(psi)**2
			
		step1 = fourier_transform(np.exp(-1j*(V + interactions)*dt/2) * psi)
		step2 = inverse_fourier(np.exp(-1j*ks**2*dt) * step1)
		psi = np.exp(-1j*(V + interactions)*dt/2) * step2 
		
		psisq = abs(psi)**2

		
	return psisquareds


#global variables for extent of modelling 

n_points = 2000 #number of timesteps (also number of x steps)
dt = 0.01
ts = np.linspace(0,20,n_points) #only for phase 
xs = np.linspace(-10, 10, n_points)
L = xs[-1] - xs[0] #box length 
dx = L/n_points


#generate initial psi
def gaussian():
	psi_initial = np.exp(-xs**2)
	return psi_initial

def soliton(start, velocity, initial_phase, family_param): #frequency = 1
	psi_initial = (2*family_param)**0.5 * np.exp(1j*velocity*(xs-start)) * np.exp(1j*initial_phase) / np.cosh((xs-start) * (family_param)**0.5)
	#psi_initial = (2*family_param)**0.5 * np.exp(1j*velocity*(xs-start)) / np.cosh((xs-start) * (family_param)**0.5) #no phase
	return psi_initial
	
def two_solitons(x1,x2,v1,v2, phase1, phase2, family_param):
    psi_combined = soliton(x1,v1,phase1,family_param) + soliton(x2,v2,phase2,family_param) 
    return psi_combined

#test accuracy of modelling with a special case
def accuracy_test(g_test,t):
    psisquareds_test = schrodinger_time_evolution(soliton(0,0), 0, g_test)
    t = int(t/dt) #taking timestep into account and making sure it's an integer as it is used for indexing 
    norm = np.sum(psisquareds_test[750:1250,:],0) # sum over axis 0, e.g. space between -2 and 2
    accuracy = abs(norm[0] - norm[t])*100/norm[0] #percentage change in norm 
    return accuracy


######### TESTS #########

g_test = -0.025
g_test2 = -0.05
g_test3 = -1

rel_phase1 = 0
rel_phase2 = np.pi

testpsi = schrodinger_time_evolution(two_solitons(-6,6,L/2,-L/2,0,rel_phase1,1), 0, g_test3) 
testpsi2 = schrodinger_time_evolution(two_solitons(-6,6,L/2,-L/2,0,rel_phase2,1), 0, g_test3) 

testpsi3 = schrodinger_time_evolution(two_solitons(-6,6,L/2,-L/2,0,rel_phase1,50), 0, g_test3) 
testpsi4 = schrodinger_time_evolution(two_solitons(-6,6,L/2,-L/2,0,rel_phase2,50), 0, g_test3) 

difference1 = testpsi - testpsi2
difference2 = testpsi3 - testpsi4


########## PLOTS ############

pyplot.figure(figsize=(15,10)) 
#pyplot.suptitle("The effect of inter-atom interactions " + r'$g |\psi|^2$' + " on a soliton model", fontsize=16)	

pyplot.subplot(231)
pyplot.imshow(np.transpose(testpsi), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title("Relative phase {}".format(rel_phase1))
#pyplot.savefig('rel_phase_0.png')
 
pyplot.subplot(232)
pyplot.imshow(np.transpose(testpsi2), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title("Relative phase {}".format(rel_phase2))
#pyplot.savefig('rel_phase_pi.png')

#one minus the other to show where fringes are 
pyplot.subplot(233)
pyplot.plot(difference1[:,1200])
pyplot.xlabel("Space")
pyplot.title("Difference at t=12")
#pyplot.savefig('difference.png')

pyplot.subplot(234)
pyplot.imshow(np.transpose(testpsi3[600:1400,900:1600]), extent=(-3,3,9,16), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title("Relative phase {}".format(rel_phase1))
#pyplot.savefig('rel_phase_0.png')
 
pyplot.subplot(235)
pyplot.imshow(np.transpose(testpsi4[600:1400,900:1600]), extent=(-3,3,9,16), origin='lower', cmap='viridis') 
pyplot.xlabel("Space")
pyplot.ylabel("Time")
pyplot.title("Relative phase {}".format(rel_phase2))
#pyplot.savefig('rel_phase_pi.png')

#one minus the other to show where fringes are 
pyplot.subplot(236)
pyplot.plot(difference2[400:1600,1200])
pyplot.xlabel("Space")
pyplot.title("Difference at t=12")
#pyplot.savefig('difference.png')

pyplot.savefig('extensionpic2.png')
pyplot.show()