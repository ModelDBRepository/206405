# Charles Heller
# Hodgkin-Huxley conduction velocity
#"A Hodgkin-Huxely Model for conduction velocity in the medial giant fiber of the earthworm, Lumbricus Terrestris"

from numpy import *
from matplotlib.pyplot import *
import numpy as np
import matplotlib.pyplot as plt
import csv

ga = [] 
gc = []
successGc1 = []
successGa1 = []
successGc2 = []
successGa2 = []
successGc3 = []
successGa3 = []
failT = []
failA = []
latency1 = []
latency2 = []
latency3 = []
c = []

#choose number of iterations you want to test
numTrials = 50
for trial in range(0,numTrials):
	# Parameters
	T = 25
	dt = .001
	time = np.arange(0, T + dt, dt)

	# Conductances
	Gna = 120
	Gk = 36
	gl = .3
	d = np.random.normal(12.5,7.5)     #micrometers (axon diameter)
	ga.append(1.024*(d/2)**2)  # axoplasmic conductance
	gc.append(np.random.normal(40,20))  #gc is not dependent on diameter in this model
	gr = .01  # reversal current conductance
	# Reversal potentials
	Ena = 115
	Ek = -12
	El = 10.613

	# capacitance
	c.append(.16*(d/2))    #capactitance is dependent on diameter

	# stimulus (arbitrary stimulus)
	I = zeros(len(time))
	for i, t in enumerate(time):
	    if 5 < t < 7: I[i] = 10

	# number of compartments
	compart = 3

	class cell:
	    # Constructor
	    def __init__(self, i, numComp, Vrest):
		self.num = i
		self.compartments = []
		for j in range(0, numComp):
		    self.compartments.append(compartment(j, Vrest))


	class compartment:
	    # Constructor
	    def __init__(self, i, Vrest):
		self.num = i
		self.Vrest = Vrest
		self.Vm = zeros(len(time))
		self.Vm[0] = Vrest
		self.m = self.m_inf((Vrest))
		self.h = self.h_inf((Vrest))
		self.n = self.n_inf((Vrest))

	    # Methods

	    # rate constants
	    def Am(self, a):
		v1 = a
		if v1 != 25:
		    am = .1 * (-v1 + 25) / (exp((-v1 + 25) / 10) - 1)
		else:
		    am = 1
		return am

	    def Bm(self, a):
		v2 = a
		bm = 4 * exp((-v2) / 18)
		return bm

	    def m_inf(self, a):
		out = self.Am(a) / (self.Am(a) + self.Bm(a))
		return out

	    def Ah(self, a):
		v3 = a
		ah = 0.07 * exp((-v3) / 20)
		return ah

	    def Bh(self, a):
		v4 = a
		bh = 1 / (exp((-v4 + 30) / 10) + 1)  # 1/(exp((-v4+30)/10) + 1)
		return bh

	    def h_inf(self, a):
		out1 = self.Ah(a) / (self.Ah(a) + self.Bh(a))
		return out1

	    def An(self, a):
		v5 = a
		if v5 != 10:
		    an = .01 * (-v5 + 10) / (exp((-v5 + 10) / 10) - 1)
		else:
		    an = 0.1
		return an

	    def Bn(self, a):
		v6 = a
		bn = .125 * exp((-v6) / 80)
		return bn

	    def n_inf(self, a):
		out2 = self.An(a) / (self.An(a) + self.Bn(a))
		return out2

	    # step function
	    def step(self, i, prev, next):
		if self.num == 0:
		    self.gna = Gna * self.m ** 3 * self.h
		    self.gk = Gk * self.n ** 4
		    self.gl = gl

		    self.m += (self.Am(self.Vm[i - 1]) * (1 - self.m) - self.Bm(self.Vm[i - 1]) * self.m) * dt
		    self.h += (self.Ah(self.Vm[i - 1]) * (1 - self.h) - self.Bh(self.Vm[i - 1]) * self.h) * dt
		    self.n += (self.An(self.Vm[i - 1]) * (1 - self.n) - self.Bn(self.Vm[i - 1]) * self.n) * dt


		    v = self.Vm[i - 1]
		    vNext = next.Vm[i - 1]
		    self.Vm[i] = v + (I[i] - self.gna * (v - Ena) - self.gk * (v - Ek) - self.gl * (v - El) + gr * (
		    vNext - v)) * dt / c[trial]

		if self.num > 0 and self.num < compart - 1:
		    self.gna = Gna * self.m ** 3 * self.h
		    self.gk = Gk * self.n ** 4
		    self.gl = gl

		    self.m += (self.Am(self.Vm[i - 1]) * (1 - self.m) - self.Bm(self.Vm[i - 1]) * self.m) * dt
		    self.h += (self.Ah(self.Vm[i - 1]) * (1 - self.h) - self.Bh(self.Vm[i - 1]) * self.h) * dt
		    self.n += (self.An(self.Vm[i - 1]) * (1 - self.n) - self.Bn(self.Vm[i - 1]) * self.n) * dt

		    vPrev = prev.Vm[i - 1]
		    v = self.Vm[i - 1]
		    vNext = next.Vm[i - 1]
		    self.Vm[i] = v + (-self.gna * (v - Ena) - self.gk * (v - Ek) - self.gl * (v - El) + ga[trial] * (vPrev - v) + gr * (
		                     vNext - v)) * dt / c[trial]

		  

		if self.num == compart - 1:
		    self.gna = Gna * self.m ** 3 * self.h
		    self.gk = Gk * self.n ** 4
		    self.gl = gl

		    self.m += (self.Am(self.Vm[i - 1]) * (1 - self.m) - self.Bm(self.Vm[i - 1]) * self.m) * dt
		    self.h += (self.Ah(self.Vm[i - 1]) * (1 - self.h) - self.Bh(self.Vm[i - 1]) * self.h) * dt
		    self.n += (self.An(self.Vm[i - 1]) * (1 - self.n) - self.Bn(self.Vm[i - 1]) * self.n) * dt

		    v = self.Vm[i - 1]
		    vPrev = prev.Vm[i - 1]
		    self.Vm[i] = v + (-self.gna * (v - Ena) - self.gk * (v - Ek) - self.gl * (v - El) + ga[trial] * ( vPrev - v)) * dt / c[trial]

		# gapJunction step method

	    def gapJunc(self, i, next, prevCell):
		self.gna = Gna * self.m ** 3 * self.h
		self.gk = Gk * self.n ** 4
		self.gl = gl

		self.m += (self.Am(self.Vm[i - 1]) * (1 - self.m) - self.Bm(self.Vm[i - 1]) * self.m) * dt
		self.h += (self.Ah(self.Vm[i - 1]) * (1 - self.h) - self.Bh(self.Vm[i - 1]) * self.h) * dt
		self.n += (self.An(self.Vm[i - 1]) * (1 - self.n) - self.Bn(self.Vm[i - 1]) * self.n) * dt


		v = self.Vm[i - 1]
		vNext = next.Vm[i - 1]
		previousCell = prevCell.Vm[i - 1]
		self.Vm[i] = v + (-self.gna * (v - Ena) - self.gk * (v - Ek) - self.gl * (v - El) + gr * (vNext - v) + gc[trial] * (previousCell - v)) * dt / c[trial]


	cells = []

	# creating the cells

	for i in range(0, 2):
	    if i == 0:
		cells.append(cell(1, compart, 0))

	    elif i > 0:
		cells.append(cell(i + 1, compart, 0))


	# now updating the compartments for each cell


	for t in range(1, len(time)):

	    for cel in range(0, len(cells)):

		for comp in cells[cel].compartments:
	
		    if comp.num == 0:
		        if cel == 0:
		            comp.step(t, 0, cells[cel].compartments[comp.num + 1])
		        else:	
		            comp.gapJunc(t, cells[cel].compartments[comp.num + 1], cells[cel - 1].compartments[compart - 1])
		    elif comp.num == compart - 1:
		        comp.step(t, cells[cel].compartments[comp.num - 1], 0)
		    else:
		        comp.step(t, cells[cel].compartments[comp.num - 1], cells[cel].compartments[comp.num + 1])



	
	t1 = float(np.argmax(cells[1].compartments[0].Vm) - np.argmax(cells[0].compartments[0].Vm)) #latency between spikes
	
	
	t1 = t1/(1/dt)
	print(t1)	
	print(trial)	#printing latency and iteration number to terminals to monitor progress of code
	if .001 < t1 < .004:
		successGc1.append(gc[trial])
		successGa1.append(ga[trial])
		latency1.append(t1)
	elif .003 < t1 < .006:
		successGc2.append(gc[trial])
		successGa2.append(ga[trial])
		latency2.append(t1)
	elif .005 < t1 < .008:
		successGc3.append(gc[trial])
		successGa3.append(ga[trial])
		latency3.append(t1)		 	
	else:
		failT.append(gc[trial])

		failA.append(ga[trial])	

figure()
xlabel('Gap junction coupling')
ylabel('Axoplasmic conductance')
title('Sensitivity analysis')
scatter(successGc1, successGa1, color='b')
scatter(successGc2, successGa2, color= 'g')
scatter(successGc3, successGa3, color= 'purple')
scatter(failT, failA, color = 'r')
show()


output = np.column_stack((successGc1,successGa1,latency1))

out = open("output1_gRand.csv","wb")
data = csv.writer(out)
for row in output:
	data.writerow(row)
out.close()

output2 = np.column_stack((successGc2,successGa2,latency2))

out = open("output2_gcRand.csv","wb")
data = csv.writer(out)
for row in output2:
	data.writerow(row)
out.close()

output3 = np.column_stack((successGc3,successGa3,latency3))

out = open("output3_gcRand.csv","wb")
data = csv.writer(out)
for row in output3:
	data.writerow(row)
out.close()

		
