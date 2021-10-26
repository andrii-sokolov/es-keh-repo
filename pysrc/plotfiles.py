import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

exp0_05 = open("Data/Experimental_Data/20150909_P-f 0V 0.5g 6.67MOhm.csv",'r')
exp0_10 = open("Data/Experimental_Data/20150522_P-f_0V 1.00g 6.67MOhm.csv",'r')
exp0_20 = open("Data/Experimental_Data/20150522_P-f_0V 2.00g 6.67MOhm.csv",'r')
exp25_05 = open("Data/Experimental_Data/20150522_P-f 25V 0.5g 6.67MOhm.csv",'r')
exp25_10 = open("Data/Experimental_Data/20150522_P-f 25V 1.0g 6.67MOhm.csv",'r')
exp25_20 = open("Data/Experimental_Data/20150522_P-f 25V 2.0g 6.67MOhm.csv",'r')

the0_05 = open("export_0V_0.5g.csv",'r')
the0_10 = open("export_0V_1.0g.csv",'r')
the0_20 = open("export_0V_2.0g.csv",'r')
the25_05 = open("export_25V_0.5g.csv",'r')
the25_10 = open("export_25V_1.0g.csv",'r')
the25_20 = open("export_25V_2.0g.csv",'r')

exp0_05_f = []
exp0_10_f = []
exp0_20_f = []
exp0_05_w = []
exp0_10_w = []
exp0_20_w = []
exp25_05_f = []
exp25_10_f = []
exp25_20_f = []
exp25_05_w = []
exp25_10_w = []
exp25_20_w = []

the0_05_f = []
the0_10_f = []
the0_20_f = []
the0_05_w = []
the0_10_w = []
the0_20_w = []
the25_05_f = []
the25_10_f = []
the25_20_f = []
the25_05_w = []
the25_10_w = []
the25_20_w = []


for line in exp0_05:
	exp0_05_f.append(float(line.split(',')[0]))
	exp0_05_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the0_05:
	the0_05_f.append(float(line.split(',')[0]))
	the0_05_w.append(float(line.split(',')[1])*1e9)

for line in exp0_10:
	exp0_10_f.append(float(line.split(',')[0]))
	exp0_10_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the0_10:
	the0_10_f.append(float(line.split(',')[0]))
	the0_10_w.append(float(line.split(',')[1])*1e9)

for line in exp0_20:
	exp0_20_f.append(float(line.split(',')[0]))
	exp0_20_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the0_20:
	the0_20_f.append(float(line.split(',')[0]))
	the0_20_w.append(float(line.split(',')[1])*1e9)

for line in exp25_05:
	exp25_05_f.append(float(line.split(',')[0]))
	exp25_05_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the25_05:
	the25_05_f.append(float(line.split(',')[0]))
	the25_05_w.append(float(line.split(',')[1])*1e9)

for line in exp25_10:
	exp25_10_f.append(float(line.split(',')[0]))
	exp25_10_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the25_10:
	the25_10_f.append(float(line.split(',')[0]))
	the25_10_w.append(float(line.split(',')[1])*1e9)

for line in exp25_20:
	exp25_20_f.append(float(line.split(',')[0]))
	exp25_20_w.append(float(line.split(',')[2])/float(line.split(',')[0])*1e9)

for line in the25_20:
	the25_20_f.append(float(line.split(',')[0]))
	the25_20_w.append(float(line.split(',')[1])*1e9)


fig, axs = plt.subplots(3,2,figsize=(10,5))

axs[0][0].plot(exp0_05_f,exp0_05_w)
axs[0][0].plot(the0_05_f,the0_05_w, color = 'red')
axs[0][0].set_xlabel("Frequancy, Hz")
axs[0][0].set_ylabel("Energy per cycle ,nJ")
axs[0][0].set_ylim(0,3)
axs[0][0].annotate('$A_{ext} = 0.5g \\quad V_{bias} = 21 V$', xy = (0.6,0.8), xycoords = 'axes fraction')
axs[1][0].plot(exp0_10_f,exp0_10_w)
axs[1][0].plot(the0_10_f,the0_10_w, color = 'red')
axs[1][0].set_xlabel("Frequancy, Hz")
axs[1][0].set_ylabel("Energy per cycle ,nJ")
axs[1][0].set_ylim(0,4)
axs[1][0].annotate('$A_{ext} = 1.0g \\quad V_{bias} = 21 V$', xy = (0.6,0.8), xycoords = 'axes fraction')
axs[2][0].plot(exp0_20_f,exp0_20_w)
axs[2][0].plot(the0_20_f,the0_20_w, color = 'red')
axs[2][0].set_xlabel("Frequancy, Hz")
axs[2][0].set_ylabel("Energy per cycle ,nJ")
axs[2][0].set_ylim(0,5)
axs[2][0].annotate('$A_{ext} = 2.0g \\quad V_{bias} = 21 V$', xy = (0.6,0.8), xycoords = 'axes fraction')

axs[0][1].plot(exp25_05_f,exp25_05_w)
axs[0][1].plot(the25_05_f,the25_05_w, color = 'red')
axs[0][1].set_xlabel("Frequancy, Hz")
axs[0][1].set_ylabel("Energy per cycle ,nJ")
axs[0][1].set_ylim(0,15)
axs[0][1].annotate('$A_{ext} = 0.5g \\quad V_{bias} = 46 V$', xy = (0.6,0.8), xycoords = 'axes fraction')
axs[1][1].plot(exp25_10_f,exp25_10_w)
axs[1][1].plot(the25_10_f,the25_10_w, color = 'red')
axs[1][1].set_xlabel("Frequancy, Hz")
axs[1][1].set_ylabel("Energy per cycle ,nJ")
axs[1][1].set_ylim(0,25)
axs[1][1].annotate('$A_{ext} = 1.0g \\quad V_{bias} = 46 V$', xy = (0.6,0.8), xycoords = 'axes fraction')
axs[2][1].plot(exp25_20_f,exp25_20_w, label = "Experimental Data")
axs[2][1].plot(the25_20_f,the25_20_w, label = "Simmulations", color = 'red')
axs[2][1].set_xlabel("Frequancy, Hz")
axs[2][1].set_ylabel("Energy per cycle ,nJ")
axs[2][1].set_ylim(0,30)
axs[2][1].annotate('$A_{ext} = 2.0g \\quad V_{bias} = 46 V$', xy = (0.6,0.8), xycoords = 'axes fraction')

plt.legend()
plt.show()
