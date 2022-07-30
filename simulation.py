from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
import graph
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import math
import mplcursors
from tkinter import filedialog as fd
import csv
from pathlib import Path
import time as time_
import numpy as np

global dt
dt = float
global tt
tt = float
global cv
cv = None

global all_data
all_data = {"Simulation":()}

root = Tk()
root.title("Cyclic Voltammetry Simulation")
root.state("zoomed")

main = ttk.Frame(root,padding="10 10 10 10")
main.grid(column=0,row=0,sticky=(N,W,E,S))
main2 = ttk.Frame(root,padding="10 10 10 10")
main2.grid(column=1,row=0,sticky=(N,W,E,S))



E1 = StringVar(value=0.2)
E2 = StringVar(value=-0.3)
v = StringVar(value=0.04)
datapts = StringVar(value=200)
dispts = StringVar(value=50)
D = StringVar(value=1.0E-5)
n = StringVar(value=1)
E0 = StringVar(value=0.013)
k0 = StringVar(value=1.0)
cox = StringVar(value=6.1e-5)
k1 = StringVar(value=0.075)
k2 = StringVar(value=1.6E3)
alpha = StringVar(value=0.5)
A = StringVar(value=2.54e-2)
T = StringVar(value=293.15)
rot = StringVar(value=0)
vis = StringVar(value=1.0)

comp_list = {
	"E1":(E1,3,"Starting Potential (E\N{SUBSCRIPT ONE})","V",False,
		"E\N{SUBSCRIPT ONE}"),
	"E2":(E2,4,"Switching Potential (E\N{SUBSCRIPT TWO})","V",False,
		"E\N{SUBSCRIPT TWO}"),
	"v":(v,9,"Scan Rate (\N{GREEK SMALL LETTER NU})","V/s",True,0.0,0.08,
		"\N{GREEK SMALL LETTER NU}"),
	"t":(datapts,1,"# Data Points (n\N{LATIN SUBSCRIPT SMALL LETTER T})",
		"",False,"n\N{LATIN SUBSCRIPT SMALL LETTER T}"),
	"d":(dispts,2,"# Distance Increments"+
		" (n\N{LATIN SUBSCRIPT SMALL LETTER X})","",False,
		"n\N{LATIN SUBSCRIPT SMALL LETTER X}"),
	"n":(n,5,"# Electrons (n)","",False,"n"),
	"E0":(E0,10,"Std. Reduction Potential"+
		" (E\N{SUBSCRIPT ZERO})","V",True,-1.987,2.013,"E\N{SUBSCRIPT ZERO}"),
	"k0":(k0,11,"Std. Rate Constant"+
		" (k\N{SUBSCRIPT ZERO})","cm/s",True,0,2,"k\N{SUBSCRIPT ZERO}"),
	"c":(cox,12,"Concentration (C)","mol/L",True,5.1e-5,7.1e-5,"C"),
	"k1":(k1,13,"1st Order Rate Constant"+
		" (k\N{SUBSCRIPT ONE})",
		"s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}",True,0,0.15,
		"k\N{SUBSCRIPT ONE}"),
	"a":(alpha,15,"Transfer Coef."+
		" (\N{GREEK SMALL LETTER ALPHA})","",True,0,1,
		"\N{GREEK SMALL LETTER ALPHA}"),
	"A":(A,16,"Electrode Area (A)","cm\N{SUPERSCRIPT TWO}",True,0,5.08e-2,"A"),
	"T":(T,17,"Temperature (T)","K",True,273.15,313.15,"T"),
	"rot":(rot,18,"Rotation Rate"+
		" (\N{GREEK SMALL LETTER OMEGA})","rpm",True,0,1500,
		"\N{GREEK SMALL LETTER OMEGA}"),
	"vis":(vis,19,"Kin. Viscosity (v)","cm\N{SUPERSCRIPT TWO}/s",True,0,
		2,"v"),
	"D":(D,8,"Diffusion Coef. (D)","cm\N{SUPERSCRIPT TWO}/s",
		True,0,2.0e-5,"D"),
	"k2":(k2,14,"2nd Order Rate Constant"+
		" (k\N{SUBSCRIPT TWO})","L/mol s",True,0,3.2e3,"k\N{SUBSCRIPT TWO}")
	}
	
entries = []
labels = []
scales = {}
units = []

for i,j in comp_list.items():
	var = j[0]
	rw = j[1]
	title = j[2]
	unit = j[3]
	scal = j[4]
	
	e = ttk.Entry(main,width=7,textvariable = var)
	e.grid(column=2,row=rw,sticky=(N,S,E,W),pady="5",padx="5")
	entries.append(e)
	
	l = ttk.Label(main,text = title)
	l.grid(column=1,row=rw,sticky=E,pady="5",padx="5")
	labels.append(l)
	
	u = ttk.Label(main,text=unit)
	u.grid(column=3,row=rw,sticky=W,pady="5",padx="5")
	units.append(u)
	
	if scal == True:
		f1 = j[5]
		t1 = j[6]
		abr = j[7]
		s = ttk.Scale(main,from_=f1,to=t1,variable = var)
		s.grid(column=4,row=rw,pady="5",padx="5")
		sca = [s,f1,t1]
		scales[abr]=sca
		
		
def scale_config(sc):

	

	edit_win = Toplevel(root)
	edit_win.geometry("200x200")
	edit_win.title("Edit")
	def close(*args):
		sca = sc[val.get()]
		fr = float(mn.get())
		tu = float(mx.get())
		print(fr)
		sca[0].configure(from_=fr)
		sca[0].configure(to=tu)
		sca[1] = fr
		sca[2] = tu
		
		edit_win.destroy()

	
	edit_fram = ttk.Frame(edit_win,padding="10 10 10 10")
	edit_fram.grid(column=0,row=0,sticky=(N,S,E,W))
	
	edit_lab = ttk.Label(edit_win,text = "Edit slider")
	edit_lab.grid(column=1,columnspan=2,row=1,padx="5",pady="5")
	
	mn = StringVar()
	mn_ent = ttk.Entry(edit_win,width=10,textvariable=mn)
	mn_ent.grid(column=2,row=2,padx="5",pady="5")
	mn_lab = ttk.Label(edit_win,text = "Minimum")
	mn_lab.grid(column=1,row=2,padx="5",pady="5")
	
	mx = StringVar()
	mx_ent = ttk.Entry(edit_win,width=10,textvariable=mx)
	mx_ent.grid(column=2,row=3,padx="5",pady="5")
	mx_lab = ttk.Label(edit_win,text = "Maximum")
	mx_lab.grid(column=1,row=3,padx="5",pady="5")
	
	edit_button = ttk.Button(edit_win,width=10,text="Enter",
		command=close)
	edit_button.grid(column=2,row=5,padx="5",pady="5")
	
	val = StringVar()
	vals = []
	for m,n in sc.items():
		vals.append(m)
	val_menu = ttk.OptionMenu(edit_win,val,vals[0],*vals)
	val_menu.grid(column=2,row=4,padx="5",pady="5")
	val_lab = ttk.Label(edit_win,text = "Parameter")
	val_lab.grid(column=1,row=4,padx="5",pady="5")
	
	mn_ent.focus_set()
	edit_win.bind("<Return>",close)
	
b = ttk.Button(main2,text = "Edit Scales",
	command = lambda:scale_config(scales))
b.grid(column=6,row=1,pady="5",padx="5")
b.configure(width=10)	

	

error = ttk.Label(main2, text = "",foreground="red")
error.grid(column=9,row=4,pady="5",padx="5")

gr = Figure(figsize=(5,4),dpi=100, tight_layout=True)
canvas = FigureCanvasTkAgg(gr,main2)  
canvas.get_tk_widget().grid(column=1,row=2,rowspan=2,columnspan = 8, 
	sticky=(N,S,E,W),padx="5",pady="5")
plot1 = gr.add_subplot(111)
plot1.set_xlabel("Potential (V)")
plot1.set_ylabel("Current (A)")
plot1.set_title("Cyclic Voltammogram")
plot1.plot()

gr2 = Figure(figsize = (3,3), dpi = 75, tight_layout=True)
canvas2 = FigureCanvasTkAgg(gr2,main2)
canvas2.get_tk_widget().grid(column=9,row=2,rowspan=1,sticky=NW,padx="5",
	pady="5")
plot2 = gr2.add_subplot(111)
plot2.set_xlabel("Time (s)")
plot2.set_ylabel("Potential (V)")
plot2.set_title("Potential Waveform")
plot2.plot()

gr3 = Figure(figsize = (3,3), dpi = 75, tight_layout=True)
canvas3 = FigureCanvasTkAgg(gr3,main2)
canvas3.get_tk_widget().grid(column=9,row=3,rowspan=1,sticky=SW,padx="5",
	pady="5")
plot3 = gr3.add_subplot(111)
plot3.set_xlabel("Distance (m)")
plot3.set_ylabel("Concentration (M)")
plot3.set_title("Conc. Profile")
plot3.plot()

time = StringVar(value="0.0")

canvas.draw_idle()
canvas2.draw_idle()
canvas3.draw_idle()

photo = PhotoImage(file = "duck.png")
root.iconphoto(False,photo)

mechanism = StringVar(value = "E")
mech = ttk.OptionMenu(main2,mechanism,"E",*["E","EC","EC2"])
mech.grid(column=2,row=4,padx="5",pady="5",sticky=W)
mech.configure(width=7)
mech_label = ttk.Label(main2,text="Mechanism:")
mech_label.grid(column=1,row=4,padx="5",pady="5",sticky=E)

conv = StringVar(value = "US")
convention = ttk.OptionMenu(main2,conv,"US",*["US","IUPAC"])
convention.grid(column=4,row=4,padx="5",pady="5",sticky=W)
convention.configure(width=7)
conv_label = ttk.Label(main2,text="Convention:")
conv_label.grid(column=3,row=4,padx="5",pady="5",sticky=E)

start = StringVar(value = "Ox")
start_spec = ttk.OptionMenu(main2,start,"Ox",*["Ox","Red"])
start_spec.grid(column=6,row=4,padx="5",pady="5",sticky=W)
start_spec.configure(width=7)
start_label = ttk.Label(main2,text="Starting:")
start_label.grid(column=5,row=4,padx="5",pady="5",sticky=E)

opt_lab = ttk.Label(main2,text = "")
opt_lab.grid(column=9,row=1,padx="5",pady="5",sticky=W)

def save(*args):
	try:
		xy = all_data["Simulation"][2]
		graph.save_data(xy[:3])
	except:
		pass
	
def name_prompt(x,y):
	
	name_win = Toplevel(root)
	name_win.geometry("200x100")
	name_win.title("Name")
	
	def close(*args):
		if name.get() != "":
			new_data = (x,y)
			
			name_win.destroy()
			
			p_name = name.get()
		
			all_data[p_name]=new_data
			
			
			run_cv()
		else:
			pass
	
	nam_fram = ttk.Frame(name_win,padding="10 10 10 10")
	nam_fram.grid(column=0,row=0,sticky=(N,S,E,W))
	
	nam_lab = ttk.Label(name_win,text = "Name plot")
	nam_lab.grid(column=1,row=1,padx="5",pady="5")
	
	name = StringVar()
	nam_ent = ttk.Entry(name_win,width=10,textvariable=name)
	nam_ent.grid(column=1,row=2,padx="5",pady="5")
	
	name_button = ttk.Button(name_win,width=10,text="Enter",
		command=close)
	name_button.grid(column=1,row=3,padx="5",pady="5")
	nam_ent.focus_set()
	name_win.bind("<Return>",close)
	
def add_data(*args):
	xy = all_data["Simulation"][2]
	name_prompt(xy[0],xy[1])
	


def clear_data(*args):

	sim = all_data["Simulation"]
	all_data.clear()
	all_data["Simulation"]=sim
	opt_lab.configure(text="")
	run_cv()

def upload_data(*args):
	try:
		filename = fd.askopenfilename()
		new_c = []
		new_p = []
		with open(filename,"r",newline="") as csvfile:
			read = csv.reader(csvfile,delimiter = "\t")
			for i in read:
				if  i[0][0].isdigit()==True:
					new_p.append(float(i[1]))
					new_c.append(float(i[2]))
				else:
					continue
		csvfile.close()
		
		
		name_prompt(new_p,new_c)
		
	except:
		pass
		
		
def fit_data(*args):
	

		
	fit_window = Toplevel(root)
	fit_window.title("Fitting Data")
	fit_main = ttk.Frame(fit_window,padding="10 10 10 10")
	fit_main.grid(column=0,row=0,sticky=(N,W,E,S))
	
	num = StringVar()
	num_ent = ttk.Entry(fit_main,width=7,textvariable=num)
	num_ent.grid(column=2,row=1,sticky=W,padx="5",pady="5")
	num_lab = ttk.Label(fit_main,text="# Points to Test")
	num_lab.grid(column=1,row=1,sticky=E,padx="5",pady="5")
	
	min_ = StringVar()
	min_ent = ttk.Entry(fit_main,width=7,textvariable=min_)
	min_ent.grid(column=2,row=2,sticky=W,padx="5",pady="5")
	min_lab = ttk.Label(fit_main,text="Minimum")
	min_lab.grid(column=1,row=2,sticky=E,padx="5",pady="5")
	
	max_ = StringVar()
	max_ent = ttk.Entry(fit_main,width=7,textvariable=max_)
	max_ent.grid(column=2,row=3,sticky=W,padx="5",pady="5")
	max_lab = ttk.Label(fit_main,text="Maximum")
	max_lab.grid(column=1,row=3,sticky=E,padx="5",pady="5")
	
	options = []
	for i,j in all_data.items():
		if i != "Simulation" and i != "Optimized":
			options.append(i)
	
	sel_plot = StringVar()
	plot_option = ttk.OptionMenu(fit_main,sel_plot,options[0],*options)
	plot_option.grid(column=2,row=4,sticky=W,padx="5",pady="5")
	plot_lab = ttk.Label(fit_main,text="Select Dataset")
	plot_lab.grid(column=1,row=4,sticky=E,padx="5",pady="5")
	
	pars = []
	for i,j in comp_list.items():
		if j[4]==True:
			pars.append(j[7])
	par_sel = StringVar()
	par_option = ttk.OptionMenu(fit_main,par_sel,pars[0],*pars)
	par_option.grid(column=2,row=5,sticky=W,padx="5",pady="5")
	par_lab = ttk.Label(fit_main,text="Select Parameter")
	par_lab.grid(column=1,row=5,sticky=E,padx="5",pady="5")
		
	prog_bar = ttk.Progressbar(fit_main,orient=HORIZONTAL,length=100,
		mode="determinate")
	prog_bar.grid(column=2,row=7,padx="5",pady="5",sticky=W)
	prog_lab = ttk.Label(fit_main,text = "Progress")
	prog_lab.grid(column=1,row=7,sticky=E,padx="5",pady="5")
	
	def fit(*args):
		try:
			exp = all_data[sel_plot.get()][1]
			
			mn = (float(min_.get()))
			mx = (float(max_.get()))
			np = (float(num.get()))
			inc = float(abs(mx-mn)/np)
			para = par_sel.get()
			
			best_rsd = 100000.0
			best_sim = tuple
			best_val = float
			par_nam = ""
			par_abr = ""
			prog = 0
			prog_bar["maximum"]=(np+1)
			for i in range(int(np)+1):
				te = {}
				test_val = (inc*float(i))+mn
				
				for k,j in comp_list.items():
					if j[4] == True:
						if j[7] == para: 
							te[k]=test_val
							par_abr = j[7]
							par_nam = k
						else:
							te[k]=float(j[0].get())
					else:
						te[k]=float(j[0].get())
				te1 = list(te.values())
				
				test_sim = graph.simulate(te1,mechanism.get(),conv.get(),
					start.get())
				new_sim = test_sim[1]
				
				rsd = graph.rsd(exp,new_sim,False)
				if rsd < best_rsd:
					best_rsd = rsd
					best_sim = test_sim[:2]
					best_val = test_val;
				prog += 1
				prog_bar["value"]=prog
				fit_window.update_idletasks()
			all_data["Optimized"] = best_sim
			
			unit = (comp_list[par_nam])[3]
			opt_lab.configure(text="Best: %s = %.3g %s"%(par_abr,best_val,unit))
			run_cv()
			fit_window.destroy()
		except:
			pass
	
	start_fit = ttk.Button(fit_main,text = "Fit", command = fit)
	start_fit.grid(column=2,row=6,sticky=W,padx="5",pady="5")
	
	num_ent.focus_set()
	fit_window.bind("<Return>",fit)
		

def simplex(*args):
	simp_win = Toplevel(root)
	simp_win.title("Simplex Optimization")
	simp_main = ttk. Frame(simp_win,padding = "10 10 10 10")
	simp_main.grid(column=0,row=0,sticky=(N,W,E,S))
	
	def optimize(*args):
		#try:
			mx1 = float(max_1.get())
			mn1 = float(min_1.get())
			mx2 = float(max_2.get())
			mn2 = float(min_2.get())

			

			x_range = mx1-mn1
			y_range = mx2-mn2
			
			par1 = par_sel1.get()
			par2 = par_sel2.get()
			

			
			dat = all_data[sel_plot.get()]
			dat_I = dat[1]
			
			x_var_key = None
			y_var_key = None
			x_var_unit = str
			y_var_unit = str
			x_var_abr = str
			y_var_abr = str
			
			
			for m,n in comp_list.items():
				if n[4] == True:
					if n[7] == par1:
						x_var_key = m
						x_var_unit = n[3]
						x_var_abr = n[7]
					elif n[7] == par2:
						y_var_key = m
						y_var_unit = n[3]
						y_var_abr = n[7]
			
			keys = (x_var_key,y_var_key)
			
			x_init = float(comp_list[x_var_key][0].get())
			y_init = float(comp_list[y_var_key][0].get())
			
			diffx = 100000.0
			diffy = 100000.0

			
			
			
			
			
			s_start = x_range/5.0
			h_start = ((y_range/10.0)*math.sqrt(3))/2
			
			
			x_init2 = x_init+s_start
			if x_init2 > mx1:
				x_init2 = x_init - s_start
			x_init3 = (x_init+x_init2)/2
			y_init2 = y_init+h_start
			if y_init2 > mx2:
				y_init2 = y_init - h_start
			

			
			vert1 = (x_init,y_init)
			vert2 = (x_init2,y_init)
			vert3 = (x_init3,y_init2)
			
			verts = {"A":vert1,"B":vert2,"C":vert3}
			
			def make_plot(verts):
				vals = []
				n_list = [(0,1),(0,2),(1,2)]
				vertsy = list(verts.values())
				for i in n_list:
					x = [vertsy[i[0]][0],vertsy[i[1]][0]]
					y = [vertsy[i[0]][1],vertsy[i[1]][1]]
					vals.append((x,y))
				return vals
			plt.ion()
			gr4.clear()
			plot5 = gr4.add_subplot(111)
			plot5.set_xlim(left=mn1,right=mx1)
			plot5.set_ylim(bottom=mn2,top=mx2)
			plot5.set_xlabel("Level of %s"%(x_var_abr))
			plot5.set_ylabel("Level of %s"%(y_var_abr))
			plot5.set_title("Simplex Optimization")
			init_plot = make_plot(verts)
			plot5.plot(init_plot[0][0],init_plot[0][1],"b")
			plot5.plot(init_plot[1][0],init_plot[1][1],"b")
			plot5.plot(init_plot[2][0],init_plot[2][1],"b")
			canvas4.draw()
			canvas4.draw_idle()
			
			
			def find_reject(verts,keys,dat_I):
				try:
					rsd_set = {}
					sim = None
					for m,n in verts.items():
						pars = []
						for i,j in comp_list.items():
							if i == keys[0]:
								pars.append(n[0])
							elif i == keys[1]:
								pars.append(n[1])
							else:
								pars.append(float(j[0].get()))
						
						sim = graph.simulate(pars,mechanism.get(),
							conv.get(),start.get())
						new_I = sim[1]
						rsd = graph.rsd(dat_I,new_I,False)
						rsd_set[m]=rsd
						
					new_graph = (sim[0],sim[1])	
					max_rsd = max(list(rsd_set.values()))
					rej_pt = str
					for m,n in rsd_set.items():
						if n == max_rsd:
							rej_pt = m
					return (rej_pt,max_rsd,rsd_set,new_graph)
				except:
					return "Error"
					
			def new_triangle(verts):
				mode = "R"
				rej = find_reject(verts,keys,dat_I)
				rej_name = rej[0]
				rsd = rej[1]
				sumx = 0.0
				rejx = 0.0
				for i,j in verts.items():
					if i != rej_name:
						sumx += j[0]
					elif i == rej_name:
						rejx = j[0]
				sumxn = sumx/2.0
				disx = sumxn-rejx
				newx = sumxn + disx
				newx_e = sumxn + 2*disx
				newx_c = sumxn - 0.5*disx
				
				sumy = 0.0
				rejy = 0.0
				for i,j in verts.items():
					if i != rej_name:
						sumy += j[1]
					elif i == rej_name:
						rejy = j[1]
				sumyn = sumy/2.0
				disy = sumyn-rejy
				newy = sumyn + disy
				newy_e = sumyn + 2*disy
				newy_c = sumyn - 0.5*disy
				newvert = (newx,newy)
				newvert_e = (newx_e,newy_e)
				newvert_c = (newx_c,newy_c)
				
				# Reflection
				verts2 = {}
				for m,n in verts.items():
					if m == rej_name:
						verts2[rej_name] = newvert
					else:
						verts2[m] = n
						
				rej_r = find_reject(verts2,keys,dat_I)
				
				if rej_r != "Error":
					rsd_r = rej_r[1]
					rsd_r = rej_r[2][rej_name]
					if rsd_r < rsd and newx < mx1 and newx > mn1 and newy < mx2 and newy > mn2:
						
						# Extension
						verts3 = {}
						for m,n in verts.items():
							if m == rej_name:
								verts3[rej_name] = newvert_e
							else:
								verts3[m] = n
						
						rej_e = find_reject(verts3,keys,dat_I)
						if rej_e != "Error":
							#rsd_e = rej_e[1]
							rsd_e = rej_e[2][rej_name]
							if rsd_e <= rsd_r and newx_e < mx1 and newx_e > mn1 and newy_e < mx2 and newy_e > mn2:
								mode = "E"
								new_rsd = rej_e[2][rej_name]
								per_diff = (abs(rsd-new_rsd)/rsd)*100
								return (verts3,mode,rsd_e,per_diff)
							else:
								mode = "R"
								new_rsd = rej_r[2][rej_name]
								per_diff = (abs(rsd-new_rsd)/rsd)*100
								return (verts2,mode,rsd_r,per_diff)
						else:
							mode = "R"
							new_rsd = rej_r[2][rej_name]
							per_diff = (abs(rsd-new_rsd)/rsd)*100
							return (verts2,mode,rsd_r,per_diff)
					else:
						# Contraction
						verts4 = {}
						for m,n in verts.items():
							if m == rej_name:
								verts4[rej_name] = newvert_c
							else:
								verts4[m] = n
						rej_c = find_reject(verts4,keys,dat_I)
						rsd_c = rej_c[1]
						rsd_c = rej_c[2][rej_name]
						mode = "C"
						new_rsd = rej_c[2][rej_name]
						per_diff = (abs(rsd-new_rsd)/rsd)*100
						return (verts4,mode,rsd_c,per_diff)
				else:
					# Contraction
					verts4 = {}
					for m,n in verts.items():
						if m == rej_name:
							verts4[rej_name] = newvert_c
						else:
							verts4[m] = n
					rej_c = find_reject(verts4,keys,dat_I)
					rsd_c = rej_c[1]
					mode = "C"
					new_rsd = rej_c[2][rej_name]
					per_diff = (abs(rsd-new_rsd)/rsd)*100
					return (verts4,mode,rsd_c,per_diff)
			
			
			
			best_coords = None
			contract = 0
			diff_count = 0
			
			#past_rsd = 10000.0
			simp_plot = None
			trial_num = 0
			while True:
				trial.configure(text="Trials: %d"%(trial_num))
				next = new_triangle(verts)
				next_verts = next[0]
				#new_rsd = next[2]
				#per_diff = (abs(past_rsd-new_rsd)/past_rsd)*100
				per_diff = next[3]
				print(next[1],per_diff)
				if next[1] == "C":
					contract += 1
				else:
					contract = 0
				if per_diff < 5:
					diff_count += 1
				else:
					diff_count = 0
				if contract >= 10 or diff_count >= 5 or trial_num > 40:
					best_verts = next_verts
					final = find_reject(best_verts,keys,dat_I)
					best_rsd_set = final[2]
					best_pt = min(best_rsd_set,key=best_rsd_set.get)
					best_coords = best_verts[best_pt]
					simp_plot = final[3]
					break
				verts = next_verts
				
				vertsy = list(verts.values())
				gr4.clear()
				cur_plot = make_plot(verts)
				plot6 = gr4.add_subplot(111)
				plot6.set_xlim(left=mn1,right=mx1)
				plot6.set_ylim(bottom=mn2,top=mx2)
				plot6.set_xlabel("Level of %s"%(x_var_abr))
				plot6.set_ylabel("Level of %s"%(y_var_abr))
				plot6.set_title("Simplex Optimization")
				plot6.plot(cur_plot[0][0],cur_plot[0][1],"b")
				plot6.plot(cur_plot[1][0],cur_plot[1][1],"b")
				plot6.plot(cur_plot[2][0],cur_plot[2][1],"b")
				
				canvas4.draw()
				canvas4.flush_events()
				
				trial_num += 1
				
				time_.sleep(0.2)
				
			all_data["Optimized"] = simp_plot
			opt_lab.configure(text = "Best: %s = %.3g %s, %s = %.3g %s"%
				(x_var_abr,best_coords[0],x_var_unit,y_var_abr,best_coords[1],
				y_var_unit))
			run_cv()
			
			
			
			simp_win.destroy()
		#except:
			#pass
	
	
	min_1 = StringVar()
	min_ent1 = ttk.Entry(simp_main,width=7,textvariable=min_1)
	min_ent1.grid(column=2,row=3,sticky=W,padx="5",pady="5")
	min_lab1 = ttk.Label(simp_main,text="Minimum")
	min_lab1.grid(column=1,row=3,sticky=E,padx="5",pady="5")
	
	max_1 = StringVar()
	max_ent1 = ttk.Entry(simp_main,width=7,textvariable=max_1)
	max_ent1.grid(column=2,row=4,sticky=W,padx="5",pady="5")
	max_lab1 = ttk.Label(simp_main,text="Maximum")
	max_lab1.grid(column=1,row=4,sticky=E,padx="5",pady="5")
	
	min_2 = StringVar()
	min_ent2 = ttk.Entry(simp_main,width=7,textvariable=min_2)
	min_ent2.grid(column=3,row=3,sticky=W,padx="5",pady="5")
	
	max_2 = StringVar()
	max_ent2 = ttk.Entry(simp_main,width=7,textvariable=max_2)
	max_ent2.grid(column=3,row=4,sticky=W,padx="5",pady="5")
	
	trial = ttk.Label(simp_main,text="")
	trial.grid(column=2,columnspan=2,row=6,padx="5",pady="5")
	
	pars = []
	for i,j in comp_list.items():
		if j[4]==True:
			pars.append(j[7])

	par_sel1 = StringVar()
	par_sel2 = StringVar()
	par_key1 = IntVar()
	kar_key2 = IntVar()
	par_option1 = ttk.OptionMenu(simp_main,par_sel1,pars[0],*pars)
	par_option1.grid(column=2,row=2,sticky=W,padx="5",pady="5")
	par_option2 = ttk.OptionMenu(simp_main,par_sel2,pars[0],*pars)
	par_option2.grid(column=3,row=2,sticky=W,padx="5",pady="5")
	par_lab1 = ttk.Label(simp_main,text="Select Parameters")
	par_lab1.grid(column=1,row=2,sticky=E,padx="5",pady="5")
	
	options = []
	for i,j in all_data.items():
		if i != "Simulation" and i != "Optimized":
			options.append(i)
	
	sel_plot = StringVar()
	plot_option = ttk.OptionMenu(simp_main,sel_plot,options[0],*options)
	plot_option.grid(column=2,columnspan=2,row=1,padx="5",pady="5")
	plot_lab = ttk.Label(simp_main,text="Select Dataset")
	plot_lab.grid(column=1,row=1,sticky=E,padx="5",pady="5")
	
	start_opt = ttk.Button(simp_main,text = "Optimize", command = optimize)
	start_opt.grid(column=2,columnspan=2,row=5,padx="5",pady="5")
	
	gr4 = Figure(figsize = (5,5), dpi = 100, tight_layout=True)
	canvas4 = FigureCanvasTkAgg(gr4,simp_main)
	canvas4.get_tk_widget().grid(column=4,row=1,rowspan=8,sticky=SW,padx="5",
		pady="5")
	plot4 = gr4.add_subplot(111)
	plot4.set_xlabel("Level of X")
	plot4.set_ylabel("Level of Y")
	plot4.set_title("Simplex Optimization")
	plot4.plot()
	canvas4.draw_idle()
	
	min_ent1.focus_set()
	simp_win.bind("<Return>",optimize)

def fit_choice(*args):
	if len(list(all_data)) > 1:
		choice_win = Toplevel(root)
		choice_win.title("Choose Fitting Method")
		choice_main = ttk.Frame(choice_win,padding = "10 10 10 10")
		choice_main.grid(column=0,row=0,sticky=(N,W,E,S))
		
		def close(*args):
			fit_data()
			choice_win.destroy()
		
		def close2(*args):
			simplex()
			choice_win.destroy()
		
		onevar = ttk.Button(choice_main,text="One Parameters",command=close)
		onevar.grid(column=1,row=1,padx="5",pady="5")
		simp = ttk.Button(choice_main,text="Two Parameters",command=close2)
		simp.grid(column=1,row=2,padx="5",pady="5")
	else:
		pass

save = ttk.Button(main2,text="Save", command = save)
save.grid(column=1,row=1,padx="5",pady="5")
save.configure(width=10)

overlay = ttk.Button(main2,text="Overlay", command = add_data)
overlay.grid(column=2,row=1,padx="5",pady="5")
overlay.configure(width=10)

clear = ttk.Button(main2,text="Clear", command = clear_data)
clear.grid(column=3,row=1,padx="5",pady="5")
clear.configure(width=10)

upload = ttk.Button(main2,text="Upload", command = upload_data)
upload.grid(column=4,row=1,padx="5",pady="5")
upload.configure(width=10)

fitting = ttk.Button(main2,text="Fitting", command = fit_choice)
fitting.grid(column=5,row=1,padx="5",pady="5")
fitting.configure(width=10)

def get_conc(tim,cox,cred,cprod):
	
	t=int(tim/dt)
	ox = cox[t]
	red = cred[t]
	prod = cprod[t]
	con = (ox,red,prod)
	return con

def run_conc(*args):
	#try:
		xy = all_data["Simulation"][2]
		
		conc = get_conc(float(time.get()),xy[4],xy[5],xy[7])
		
		dist = xy[3]
		gr3.clear()
		plot3 = gr3.add_subplot(111)
		plot3.set_xlabel("Distance (cm)")
		plot3.set_ylabel("Concentration (M)")
		plot3.set_title("Conc. Profile")
		
		ox, = plot3.plot(dist[0],conc[0],label="Ox")
		red, = plot3.plot(dist[0],conc[1],label="Red")
		if mechanism.get() != "E":
			prod, = plot3.plot(dist[0],conc[2],label="Prod")
			plot3.legend(handles=[ox,red,prod])
		else:
			plot3.legend(handles=[ox,red])
		canvas3.draw_idle()	
		
	#except:
		#pass

def cursor_note(sel):
	xy = all_data["Simulation"][2]
	pot_array = xy[0][0]
	I_array = xy[1][0]
	time_array = xy[2][0]
	pos = sel.target
	x_pos = pos[0]
	y_pos = pos[1]

	diff=1.0E30
	count = 0
	loc = 0
	for i in I_array.flat:
		dif = abs(i - y_pos) + abs(pot_array[count]-x_pos)
		if dif < diff:
			loc = count
			diff = dif
		count+=1
	ti = time_array[loc]
	time.set(ti)
	sel.annotation.set_text("Time: %.1f s"%(ti))

def run_cv(*args):
#	try:
		error.configure(text="")
		par = [E1,E2,v,datapts,dispts,n,E0,k0,cox,k1,alpha,A,T,rot,vis,D,k2]
		vl = []
		for i in par:
			j = float(i.get())
			vl.append(j)

		

		xy = graph.simulate(vl,mechanism.get(),conv.get(),start.get())
		x = xy[0]
		y = xy[1]
		t = xy[2]


		cv_plots = []
		all_data["Simulation"] = (x,y,xy)
		gr.clear()
		
		plot1 = gr.add_subplot(111)
		plot1.set_xlabel("Potential (V)")
		plot1.set_ylabel("Current (A)")
		plot1.set_title("Cyclic Voltammogram")

		if conv.get() == "IUPAC":
			plot1.invert_xaxis()
		for i,j in all_data.items():
			if i == "Optimized":
				p, = plot1.plot(j[0][0],j[1][0],"r--",label=i)
				cv_plots.append(p)
			elif i == "Simulation":
				p, = plot1.plot(j[0][0],j[1][0],label=i)
				cv_plots.append(p)
				cv = p
			else:
				p, = plot1.plot(j[0][0],j[1][0],label=i)
				cv_plots.append(p)
		plot1.legend(handles=cv_plots)
		canvas.draw_idle()


		gr2.clear()
		plot2 = gr2.add_subplot(111)
		plot2.set_xlabel("Time (s)")
		plot2.set_ylabel("Potential (V)")
		plot2.set_title("Potential Waveform")
		plot2.plot(t[0],x[0])
		canvas2.draw_idle()
		
		global tt
		tt = xy[6][0]
		
		global dt
		dt = xy[6][1]
		
		time.set(0)
		run_conc()
		
		
		cur = mplcursors.cursor(cv,hover=False)
		cur.connect("add",cursor_note)
		if mechanism.get() == "E":
			entries[9].configure(state = "disabled")
			entries[16].configure(state = "disabled")
		elif mechanism.get() == "EC":
			entries[9].configure(state = "enabled")
			entries[16].configure(state = "disabled")
		elif mechanism.get() == "EC2":
			entries[9].configure(state = "disabled")
			entries[16].configure(state = "enabled")
	#except:
	#	error.configure(text = "Error: could not update")
	#	pass
		
	
run_cv()

for var in list(comp_list.values()):
	var[0].trace_add("write",run_cv)
mechanism.trace_add("write",run_cv)
conv.trace_add("write",run_cv)
start.trace_add("write",run_cv)
time.trace_add("write",run_conc)

root.mainloop()