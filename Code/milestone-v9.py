"""
Created on Sun Oct 28 10:15:11 2018
Version 5 
@author: Lois
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as pyplot
import time 

start = time.time() #initialising time 


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
	
	ks = 2*np.pi*getks(xs,dx) #accessing global variables. MULTIPLYING BY 2PI 
	
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
	
	
	
def particle_SHM(start, velocity, frequency):
    
    x = np.zeros(n_times)
    
    for i,t in enumerate(ts):
        x[i] = start * np.cos(frequency*t) + velocity/frequency * np.sin(frequency*t)
        
    return x


#global variables for extent of modelling 

n_times = 2000
dt = 0.001
xs = np.linspace(-10,10,n_times)
ts = np.linspace(0,2,n_times)
L = xs[-1] - xs[0]
dx = L/n_times


#generate initial psi
def gaussian():
	psi_initial = np.exp(-xs**2)
	return psi_initial

def soliton(velocity):
	psi_initial = (family_param/2)**0.5 * np.exp(1j*velocity*xs) / np.cosh(xs * family_param)
	return psi_initial

def accuracy_test(g,t):
    #uses trapezium rule
    psisquareds_test = schrodinger_time_evolution(soliton(0), 0, g)
    t = int(t/dt) #taking timestep into account and making sure it's an integer as it is used for indexing 
    norm = np.trapz(psisquareds_test[750:1250,:],xs[750:1250],axis=0)
    #print "norm at 0: {}".format(norm[0])
    #print "norm at 1.999: {}".format(norm[t])
    accuracy = abs(norm[0] - norm[t])*100/norm[0] #percentage change in norm 
    return accuracy
    
def errors(gs):
    errors = np.zeros(len(gs))
    
    for ig, g in enumerate(gs):
        errors[ig] = accuracy_test(g, 1.999)
        
        if errors[ig] == 0:
            if g != 4:
                print "what the hell"
            else:
                print "noice work" 
        
    return errors

#tests

family_param = 1
velocity = 1 #was 5/3 
testpsi = schrodinger_time_evolution(soliton(velocity), 0, 0) #no nonlinear term 
#g_test = -0.025
g_test2 = -4*family_param 
testpsi2 = schrodinger_time_evolution(soliton(velocity), 0, g_test2) #nonlinear interactions 

p = particle_SHM(0,velocity*10/5,0.000001) #effectively zero frequency to get approximation of zero potential and test function 

#testing accuracy
acc = accuracy_test(g_test2, 1.999)
print "Norm is conserved to {} % accuracy".format(acc) 


do_errors = False

if do_errors:
    
    g_centre = -4
    range_g = 0.01
    n = 25
    gs = np.linspace(g_centre - range_g/2, g_centre + range_g/2, n)
    err_gs = errors(gs)

    # REDEFINING GLOBAL VARIABLES 
    n_times = 4000
    dt = 0.002
    xs = np.linspace(-10,10,n_times)
    ts = np.linspace(0,2,n_times)
    L = xs[-1] - xs[0]
    dx = L/n_times
    
    g_centre = -4
    range_g = 0.01
    n = 25
    gs2 = np.linspace(g_centre - range_g/2, g_centre + range_g/2, n)
    err_gs2 = errors(gs2)
    
    #plot with shared x axis and different y axes
    fig, ax1 = pyplot.subplots()
    ax1.plot(gs, err_gs, 'b-')
    ax1.set_xlabel(r'$g$')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('Probability leakage with 2000x2000 points', color='b')
    
    ax2 = ax1.twinx()
    ax2.plot(gs2, err_gs2, 'r-')
    ax2.set_ylabel('Probability leakage with 4000x4000 points', color='r')
    
    fig.tight_layout()
    pyplot.savefig("errors.png")
    pyplot.show()



####### PLOTS #######

pyplot.figure(figsize=(10,6)) 
#pyplot.suptitle("The effect of inter-atom interactions " + r'$g |\psi|^2$' + " on a soliton model", fontsize=16)	

pyplot.subplot(121)
pyplot.imshow(np.transpose(testpsi), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space", fontsize=12)
pyplot.ylabel("Time", fontsize=12)
pyplot.title(r'$g=0$', fontsize=16)

pyplot.subplot(122)
pyplot.imshow(np.transpose(testpsi2), extent=(-10,10,0,20), origin='lower', cmap='viridis') 
pyplot.xlabel("Space", fontsize=12)
pyplot.ylabel("Time", fontsize=12)
pyplot.title(r'$g= {}$'.format(g_test2), fontsize=16)
#pyplot.plot(p,ts*10) #adding particle 

pyplot.savefig('milestonepic.png')
pyplot.show()

print (time.time() - start), "s"