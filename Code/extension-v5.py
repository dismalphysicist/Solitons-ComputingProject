"""
Created on Sun Oct 28 10:15:11 2018
Version 5 
@author: Lois
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors 

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
	
	psisquareds = np.zeros((n_points,n_times)) #initialising 2d array
	
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
	

def particles(start, velocity, spring_constant, mass):
    
    frequency = np.sqrt(spring_constant/mass)
    x = np.zeros(n_times)
    
    for i,t in enumerate(ts):
        x[i] = start * np.cos(frequency*t) + velocity/frequency * np.sin(frequency*t) 
    
    return x


#global variables for extent of modelling 

n_points = 4000 #number of x steps 
n_times = 4000 #number of timesteps
dt = 0.001
ts = np.linspace(0,4,n_times) 
xs = np.linspace(-20, 20, n_points)
L = xs[-1] - xs[0] #box length 
dx = L/n_points


#generate initial psi
def gaussian():
	psi_initial = np.exp(-xs**2)
	return psi_initial

def soliton(start, velocity, initial_phase, family_param): #frequency = 1
	psi_initial = (family_param/2)**0.5 * np.exp(1j*velocity*(xs-start)) * np.exp(1j*initial_phase) / np.cosh((xs-start) * family_param)
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

g = -8

family_param = -g/4

rel_phase1 = 0
rel_phase2 = np.pi

#weak axial harmonic confining potential (assume radial confinement perfect, so 1D motion)
K = 1*(2*np.pi/4)**2 #m w^2 characterises strength of harmonic potential, mass of soliton = 1
V = 1/2 * K * xs**2

#zeroPot = schrodinger_time_evolution(two_solitons(-6,6,L/6,-L/6,0,rel_phase1,0.5), 0, g) 
#zeroPotPi = schrodinger_time_evolution(two_solitons(-6,6,L/6,-L/6,0,rel_phase2,0.5), 0, g) 

velocity = 10 #L/2 originally 

#particle of mass 0.5
p = particles(-4,20,K,0.5) #velocity scaling of 1/5 * 10 (10 from time scaling)


harmPot = schrodinger_time_evolution(two_solitons(-4,4,velocity,-velocity,0,rel_phase1,family_param), V, g) 

########## PLOTS ############
pyplot.figure() #figsize=(16,5)
#pyplot.suptitle("g={}".format(g))

pyplot.imshow(np.transpose(harmPot), extent=(-20,20,0,40), origin='lower', cmap='viridis', norm=colors.SymLogNorm(linthresh=0.3, vmin=harmPot.min(), vmax=harmPot.max()))
pyplot.xlabel("Space", fontsize=12)
pyplot.ylabel("Time", fontsize=12)
pyplot.title("Relative phase {}".format(rel_phase1), fontsize=16)
pyplot.plot(p,ts*10)

pyplot.savefig('particle.png')
pyplot.show() 