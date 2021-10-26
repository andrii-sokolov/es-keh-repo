from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def func(x, C0, ds):
	return( 2.0*C0*ds**2/(ds**2 - x**2) )

file = open('capac.csv','r')
dist = []
cap = []
#force = []

for line in file:
	if(line[0]!='%'):
		dist.append(float(line.split(',')[0]))
		cap.append(1e12*float(line.split(',')[1]))
		#force.append(0.5*float(line.split(',')[2]))

popt, pcov = curve_fit(func, dist, cap)

print(popt)

dist_fit = np.linspace(-0.052,0.052,200)


fffit = func(dist_fit, *popt)

plt.plot(np.array(dist_fit)*1000, fffit, label='Theoretical fit' )
plt.plot(np.array(dist)*1000, cap,'.', label = 'COMSOL')
plt.xlabel('Displacement, $\mu m$')
plt.ylabel('Capacitance, pF')
plt.legend()
plt.show()
