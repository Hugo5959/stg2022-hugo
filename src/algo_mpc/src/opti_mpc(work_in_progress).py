from sys import path
path.append(r"home/bot/casadi-linux-py38-v3.5.5-64bit")
from casadi import *
import casadi as c
import numpy as np
from numpy import random, array, linalg, matrix, zeros, ones, ndarray, eye
from linpy import *

def f(x,u):
	return [x[3]*cos(x[2])-x[4]*sin(x[2]),
        	x[3]*sin(x[2])+x[4]*cos(x[2]),
        	x[5],
        	(u[0]+u[1])/m11+m22/m11*x[4]*x[5]+X_u/m11*x[3],
        	-m11/m22*x[3]*x[5]+Y_v/m22*x[4],
        	b*(u[0]-u[1])/m33+(m22-m11)/m33*x[3]*x[4]+N_r/m33*x[5]]


def kappa(x0,Omega,T,N):
#Kappa is the control inputs
#x0 is the current state
#Omega is the target Set
#T discretization time
#N is the horizon of prediction

	#parameters of the model
	b = 1.83/2
	m11 = 180.05
	m22 = 179
	m33 = 248.2
	X_u = 10
	Y_v = -0.5
	N_r = -0.65
	
	
	#Omega
	Ox = Omega[0]
	Oy = Omega[1]
	A = np.array([[-1,0],[0,-1],[1,0],[0,1]])
	b = np.array([[1-Ox],[1-Oy],[Ox+1],[Oy+1]])
	
	#Control parameters
	F1_max = 250 
	F2_max = 250
	
	#casadi
	opti = casadi.Opti()
	#variables d'optimisation
	X = opti.variable(6,N+1)
	x = X[0][:]
	y = X[1][:]
	psi = X[2][:]
	u = X[3][:]
	v = X[4][:]
	r = X[5][:]
	Xa = opti.variable(2,N+1)
	U = opti.variable(2,N)
	F1 = U[0][:]
	F2 = U[1][:]
	
	V = 0 # cost function
	Q = zeros(2,2)
	Q[1][1] = 30
	Q[2][2] = 30 # weight parameters for state
	R = zeros(2,2)
	R[1][1] = 0.02
	R[2][2] = 0.02 # weight parameters for control
	
	for k in range(1,N+1):
    	V = V + np.array([[x[k]-Xa[1][k]],[y[k]-Xa[2][k]]]) @ Q @ np.array([x[k]-Xa[1][k],y[k]-Xa[2][k]]) @ np.transpose(U[:][k]) @ R @ np.array(U[:][k])
    		
    	
    opti.subject_to(0<=F1<=F1_max)
    opti.subject_to(0<=F2<=F2_max)
    opti.subject_to(X[:][1]==x0)
    
    for k in range(1,N+1):
    	opti.subject_to(A @ np.array([Xa[1][k];Xa[2][k]])<=b)
    	
    for k in range(1,N+1):
    	k1 = f(X[:][k],         U[:][k])
   		k2 = f(X[:][k]+T/2*k1, U[:][k])
   		k3 = f(X[:][k]+T/2*k2, U[:][k])
   		k4 = f(X[:][k]+T*k3,   U[:][k])
   		x_next = X[:][k] + T/6*(k1+2*k2+2*k3+k4)
   		opti.subject_to(X[:][k+1]==x_next)
   		
   		
   	opti.minimize(V)
   	opti.solver('ipopt') # set numerical 
	sol = opti.solve()   # actual solve
	
	v1=sol.value(F1)
	v2=sol.value(F2)
	
	return np.array([[v1[1]],[v2[1]]])
    		
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
    	
	

