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
	
	psisquareds = np.zeros((len(psi_initial),n_points)) #initialising 2d array
	
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

#g = -5
gs = np.array([-7])
#family_param = -g/4
family_params = -gs/4

rel_phase1 = 0
rel_phase2 = np.pi

#weak axial harmonic confining potential (assume radial confinement perfect, so 1D motion)
w = 2*np.pi/4
mws = 1*w**2 #m w^2 characterises strength of harmonic potential 
V = 1/2 * mws * xs**2

#zeroPot = schrodinger_time_evolution(two_solitons(-6,6,L/6,-L/6,0,rel_phase1,0.5), 0, g) 
#zeroPotPi = schrodinger_time_evolution(two_solitons(-6,6,L/6,-L/6,0,rel_phase2,0.5), 0, g) 

velocity = 10 #L/2 originally 

for ig,g in enumerate(gs): 

    harmPot = schrodinger_time_evolution(two_solitons(-4,4,velocity,-velocity,0,rel_phase1,family_params[ig]), V, g) 
    harmPotPi = schrodinger_time_evolution(two_solitons(-4,4,velocity,-velocity,0,rel_phase2,family_params[ig]), V, g) 
    
    difference = harmPot - harmPotPi
    
    
    ########## PLOTS ############
    pyplot.figure(figsize=(16,5))
    #pyplot.suptitle("g={}".format(g))
    
    pyplot.subplot(131)
    pyplot.imshow(np.transpose(harmPot), extent=(-20,20,0,40), origin='lower', cmap='viridis', norm=colors.SymLogNorm(linthresh=0.3, vmin=harmPotPi.min(), vmax=harmPotPi.max()))
    pyplot.xlabel("Space")
    pyplot.ylabel("Time")
    pyplot.title("Relative phase {}".format(rel_phase1))
    
    pyplot.subplot(132)
    pyplot.imshow(np.transpose(harmPotPi), extent=(-20,20,0,40), origin='lower', cmap='viridis', norm=colors.SymLogNorm(linthresh=0.3, vmin=harmPotPi.min(), vmax=harmPotPi.max()))
    pyplot.xlabel("Space")
    pyplot.ylabel("Time")
    #pyplot.title("Relative phase {}".format(rel_phase2))
    pyplot.title("Relative phase " + r'$\pi$')
    
    pyplot.subplot(133)
    t=0
    if velocity == 10:
        t=16
        pyplot.plot(difference[1500:2500,1600], color="#000099")
    elif velocity == 20:
        t=15.16
        pyplot.plot(difference[1500:2500,1516], color="#000099")
    pyplot.xlabel("Space")
    pyplot.xticks(np.array([0,500,1000]),[-5.0,0.0,5.0])
    #pyplot.ylabel("Difference in " + r'$\psi^2$') # + "   " + r'$/10^{-5}$' if in units of 10^-5
    #pyplot.yticks(np.array([-0.00001,0,0.00001]),np.array([-1,0,1]))
    pyplot.title("Difference in "  + r'$\psi^2$' + " at t={}".format(t))
    
    pyplot.savefig('difference-g7.png')
    
    
    #pyplot.savefig("extension-v4-pic.png")
    pyplot.show() 

# pyplot.figure()
# pyplot.imshow(np.transpose(schrodinger_time_evolution(soliton(0,0,0,family_param), V, -4*family_param)), extent=(-20,20,0,40), origin='lower', cmap='viridis')
# pyplot.title("Harmonic potential on one soliton")
#pyplot.savefig("harmonic_confinement-1soliton.png")
# pyplot.show()