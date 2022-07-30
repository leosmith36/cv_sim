import matplotlib.pyplot as plt
import math
import csv
from tkinter import filedialog as fd
import numpy as np

F = 96485.0
R = 8.314

# Makes convection array
	
def conv_array(d,c):
	conv = []
	convect = lambda j,c : j/(1.0-0.51*math.sqrt((c[3]*(c[2]**3)*
		(c[0]**3))/c[1])*(j/math.sqrt(c[4]*c[5]**3)))
	for i in range(0,d-1,1):
		k = i + 1
		k2 = convect(k,c)
		k2-=1
		k2i = int(k2)
		if k2 >= d-1 or k2 < 0:
			k2 = d - 2
			k2i = d - 2
		dif = abs(k2-k2i)
		dif2 = 1-abs(dif)
		conv_t = (k2i,dif,dif2)
		conv.append(conv_t)
	return conv

def simulate(l,mech,conv,start):

	# Sets parameters
	
	E1 = l[0] # "E1"
	E2 = l[1] # "E2"
	v = l[2] # "v"
	n_t = l[3] # "t"
	n_d = l[4] # "d"
	n = l[5] # "n"
	E0 = l[6] # "E0"
	k0 = l[7] # "k0"
	Cox = l[8] # "c"
	k1 = l[9] # "k1"
	alpha = l[10] # "a"
	A = l[11] # "A"
	T = l[12] # "T"
	rot = l[13] # "rot"
	vis = l[14] # "vis"
	D = l[15] # "D"
	k2 = l[16] # "k2"
	
	d = int(n_d)
	t = int(n_t)
	
	t_tot = 2*(abs(E2-E1))/v
	t_inc = t_tot/n_t
	rot2 = (rot/60.0)*(2.0*math.pi)
	x_tot = 6*math.sqrt(D*t_tot)
	x_inc = x_tot/n_d
	lam = (D*t_inc)/(x_inc**2)
	calc_par = [t_tot,t_inc,x_tot,x_inc,lam]
	
	
	
	# Makes the time-potential array
	
	time_np = np.fromfunction(lambda j,i: i*t_inc, (1,t+1), dtype=float)
	pot_np = np.zeros((1,t+1),dtype=float)
	pot_np[0,0]=E1
	
	m = -1
	switched = False
	for i in range(t):
		if i == (t/2):
			m = m*-1
		if start == "Red" and switched == False:
			m = m*-1
			switched = True
		E = pot_np[0,i] + m*v*t_inc
		pot_np[0,i+1]=E

	dist_np = np.fromfunction(lambda i,j: j*x_inc, (1,d),dtype=float)
	
	# Makes the rate constant array
	
	kf = np.zeros((1,t+1))
	kb = np.zeros((1,t+1))
	count = 0
	for i in pot_np.flat:
		kfn = k0*math.exp((-1*alpha*n*F*(i-E0))/(R*T))
		kf[0,count]=kfn
		kbn = k0*math.exp(((1-alpha)*n*F*(i-E0))/(R*T))
		kb[0,count]=kbn
		count+=1
		
	# Makes the diffusion grids
	
	Cox_array = np.zeros((t+1,d))
	Cred_array = np.zeros((t+1,d))
	Cred_nr_array = np.zeros((t+1,d))
	
	Cox_array[0]=Cox
	
	Cred_array[0]=0.0

	Cred_nr_array[0]=0.0
	Jox = np.zeros((1,t))
	Jred = np.zeros((1,t))

	convt = conv_array(d,[rot2,vis,t_tot,D,lam,t])

	for j in range(t):
		
		for i in range(d):
			k = (d-1)-i
			p = k
			rt = float
			lt = float
			if i != 0:
				m = convt[k]
				k = m[0]
				rt = m[1]
				lt = m[2]
				
			if i == 0:
				Cox_array[j+1,p]=Cox
				Cred_array[j+1,p]=0.0
				Cred_nr_array[j+1,p]=0.0
				
			elif i != d-1:
			
				#Oxidized species
				Cio = Cox_array[j,k]
				Cip1 = Cox_array[j,k+1]
				Cin1 = Cox_array[j,k-1]
				Cip2 = float
				if k >= d-3 or k < 0:
					Cip2 = Cox
				else:
					Cip2 = Cox_array[j,k+2]
				Ci = lt*Cio + rt*Cip1
				Cip1 = lt*Cip1 + rt*Cip2
				Cin1 = lt*Cin1 + rt*Cio
				
				Cj = Ci + lam*(Cin1 - 2*Ci + Cip1)
				Cox_array[j+1,p]=Cj
				
				
				# Reduced species with chemical reaction
				Cio = Cred_array[j,k]
				Cip1 = Cred_array[j,k+1]
				Cin1 = Cred_array[j,k-1]
				Cip2 = float
				if k >= d-3 or k < 0:
					Cip2 = 0.0
				else:
					Cip2 = Cred_array[j,k+2]
				Ci = lt*Cio + rt*Cip1
				Cip1 = lt*Cip1 + rt*Cip2
				Cin1 = lt*Cin1 + rt*Cio
				Cj = float
				if mech == "E":		
					Cj = Ci + lam*(Cin1 - 2*Ci + Cip1) 
				elif mech == "EC":
					Cj = Ci + lam*(Cin1 - 2*Ci + Cip1)-k1*t_inc*Ci
				elif mech == "EC2":
					Cj = Ci + lam*(Cin1 - 2*Ci + Cip1)-2*k2*t_inc*(Ci**2) 
				Cred_array[j+1,p]=Cj
				
				
				# Reduced species without chemical reaction
				Cio = Cred_nr_array[j,k]
				Cip1 = Cred_nr_array[j,k+1]
				Cin1 = Cred_nr_array[j,k-1]
				Cip2 = float
				if k >= d-3 or k < 0:
					Cip2 = 0.0
				else:
					Cip2 = Cred_nr_array[j,k+2]
				Ci = lt*Cio + rt*Cip1
				Cip1 = lt*Cip1 + rt*Cip2
				Cin1 = lt*Cin1 + rt*Cio
				
				Cj = Ci + lam*(Cin1 - 2*Ci + Cip1) 
				Cred_nr_array[j+1,p]=Cj
				
			else:

				Coxi=Cox_array[j+1,p+1]
				Credi=Cred_array[j+1,p+1]
				
				if start == "Ox":
					Jox1 = (((kb[0,j+1]*Credi)-(kf[0,j+1]*Coxi)))/(1+((kf[0,j+1]*x_inc)/D)+((kb[0,j+1]*x_inc)/D))
					Jred1 = -1*Jox1
				elif start == "Red":
					Jox1 = -1*((kb[0,j+1]*Coxi)-(kf[0,j+1]*Credi))/(1+((kf[0,j+1]*x_inc)/D)+((kb[0,j+1]*x_inc)/D))
					Jred1 = -1*Jox1
				Cox0 = Coxi+(((Jox1)*x_inc)/D)
				Cred0 = Credi+(((Jred1)*x_inc)/D)
				Jox[0,j]=Jox1
				
				Jred[0,j]=Jred1

				Cox_array[j+1,0]=Cox0
				Cred_array[j+1,0]=Cred0
				Cred_nr_array[j+1,0]=Cred0
	
	Cprod_array = np.subtract(Cred_nr_array,Cred_array)
		
	# Makes current array

	
	
	I = np.zeros((1,t+1))
	I[0,0]=0.0
	count = 1
	for i in Jox.flat:
		In = -1*n*F*A*i*0.001
		if conv == "IUPAC":
			In = In*-1
		if start == "Red":
			In = In*-1
		I[0,count]=In
		count+=1
	
	
	
	# Plots the graph
	
	x = pot_np
	y = I
	graph = (pot_np,I,time_np,dist_np,Cox_array,Cred_array,calc_par,Cprod_array)
	
	return graph

def save_data(dat):
	headings = ["Time (s)","Potential (V)","Current (A)"]
	current = []
	for i in dat[1]:
		current.append("%.5e"%(i))
	filename = fd.asksaveasfilename()
	dat_r = [dat[2],dat[0],current]
	data = zip(*dat_r)
	with open(filename,"w",newline="") as csvfile:
		wr = csv.writer(csvfile,delimiter="	")
		wr.writerow(headings)
		wr.writerows(data)
	csvfile.close()

def rsd(a,b,weighted):
	if len(a)==len(b):
		summ = 0.0
		
		for i in range(len(a)):
			dif = (a[i]-b[i])**2

			if weighted == True:
				dif = dif/((a[i])**2)
			summ+=dif
		RSD = float(math.sqrt(summ/float(len(a))))

		return RSD	
	else:
		pass