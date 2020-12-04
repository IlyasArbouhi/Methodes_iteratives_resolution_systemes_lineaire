# -*- coding: utf-8 -*-
"""
Arbouhi Ilyas Geshanth Kamaleswaran
@author: arbou
"""
import numpy as np
import matplotlib.pyplot as pp 

###QUESTION 1###

def MIGenerale(M,N,b,x0,epsilon,Nitermax):
	x = x0
	Niter = 0
	erreur = epsilon +1
	while (erreur > epsilon and Niter < Nitermax):
		xp = x
		x= np.linalg.solve(M,np.dot(N,x)+b)
		x= x.reshape(-1,1)
		Niter += 1
		erreur=np.linalg.norm(x-xp)
		
	return(x,Niter,erreur)

A=np.array([[1,2,-2],[1,1,1],[2,2,1]])
M=np.array([[1,0,0],[0,1,0],[0,0,1]])
N=np.array([[0,-2,2],[-1,0,-1],[-2,-2,0]])
b=np.array([[1],[1],[1]])
x0=np.array([[0],[0],[0]])
epsilon=10**-4
Nitermax=200
omega = 1

#print(np.linalg.eig(A))
#print(MIGenerale(M,N,b,x0,epsilon,Nitermax))
"""
A=np.array([[10,-1,0],[-1,10,-2],[-2,0,10]])
b=np.array([[9],[10],[7]])
x0=np.array([[0],[0],[0]])
epsilon=10**-10
Nitermax=200
omega = 1
"""
###QUESTION 2###

def MIJacobi(A,b,x0,epsilon,Nitermax):
	M = np.diag(np.diag(A))
	N = M-A
	x,Niter,erreur = MIGenerale(M,N,b,x0,epsilon,Nitermax)
	return ([x,Niter,erreur])

x,Niter,erreur = MIJacobi(A,b,x0,epsilon,Nitermax)
#print(x,'\n',Niter,'\n',erreur)
###QUESTION 3 ###

def MIGaussSeidel(A,b,x0,epsilon,Nitermax):
	M = np.tril(A)
	N = M-A
	x,Niter,erreur = MIGenerale(M,N,b,x0,epsilon,Nitermax)
	return ([x,Niter,erreur])

x,Niter,erreur = MIGaussSeidel(A,b,x0,epsilon,Nitermax)
#print(x,'\n',Niter,'\n',erreur)

def MIRelaxation(A,b,omega,x0,epsilon,Nitermax):
    M = (1/omega)*(np.diag(np.diag(A)))+(np.tril(A)-(np.diag(np.diag(A))))
    N = ((1/omega)-1)*(np.diag(np.diag(A)))-(np.triu(A)-(np.diag(np.diag(A))))
    x,Niter,erreur = MIGenerale(M,N,b,x0,epsilon,Nitermax)
    return ([x,Niter,erreur])

x,Niter,erreur = MIRelaxation(A,b,omega,x0,epsilon,Nitermax)

def MIRichardson(A,b,omega,x0,epsilon,Nitermax):
    x,y = np.shape(A)
    M = (1/omega)*np.identity(x)
    N = M-A
    x,Niter,erreur = MIGenerale(M,N,b,x0,epsilon,Nitermax)
    return ([x,Niter,erreur])
x,Niter,erreur = MIRichardson(A,b,omega,x0,epsilon,Nitermax)
#print(x,'\n',Niter,'\n',erreur)
#print(x,'\n',Niter,'\n',erreur)
def Matrice(n):
    A = np.zeros((n,n))
    b = np.zeros((n,1))
    x0 = np.zeros((n,1))
    for i in range(n):
        for j in range(i):
            A[i,j] = 1/(12+(3*i-5*j)**2)
            b[i,0] = np.cos(i/8)
        for j in [i] :
            A[i,i] = 3
        for j in range(i+1,n):
            A[i,j] = 1/(12+(3*i-5*j)**2)
            b[i,0] = np.cos(i/8)
    return A,b,x0

def Matriceb(n):
    A = np.zeros((n,n))
    b = np.zeros((n,1))
    x0 = np.zeros((n,1))
    for i in range(n):
        for j in range(i):
            A[i,j] = 1/(1+3*(np.abs(i-j)))
            b[i,0] = np.cos(i/8)
        for j in [i] :
            A[i,i] = 3
        for j in range(i+1,n):
            A[i,j] = 1/(1+3*(np.abs(i-j)))
            b[i,0] = np.cos(i/8)
    return A,b,x0
"""
###Partie 2 
###Question 1
    
A1,b1,x0 = Matrice(100)

epsilon = []
Niterj = []
Nitergs = []
for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    J = MIJacobi(A1,b1,x0,eps,Nitermax)
    Niterj.append(J[1])
    GS = MIGaussSeidel(A1,b1,x0,eps,Nitermax)
    Nitergs.append(GS[1])
    
pp.plot(epsilon,Niterj,label = 'Jacobi')
pp.plot(epsilon,Nitergs,label = 'Gauss-Seidel')
pp.legend()
pp.title("Evolution du nb d'itérations en fonction de epsilon")
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.xscale('log')
pp.show()

###Question 2

A2,b2,x0 = Matriceb(100)

epsilon = []
Niter = []
Nitergs = []

for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    J = MIJacobi(A2,b2,x0,eps,Nitermax)
    Niter.append(J[1])
    GS = MIGaussSeidel(A2,b2,x0,eps,Nitermax)
    Nitergs.append(GS[1])
    
pp.plot(epsilon,Niter,label = 'Jacobi')
pp.plot(epsilon,Nitergs,label = 'Gauss-Seidel')
pp.title("Evolution du nb d'itérations en fonction de epsilon")
pp.legend()
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.xscale('log')
pp.show()
"""
### Question 3 
A1,b1,x0 = Matrice(100)

epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []


for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRelaxation(A1,b1,1,x0,eps,Nitermax)
    R2 = MIRelaxation(A1,b1,1.2,x0,eps,Nitermax)
    R3 = MIRelaxation(A1,b1,1.4,x0,eps,Nitermax)
    R4 = MIRelaxation(A1,b1,1.6,x0,eps,Nitermax)
    R5 = MIRelaxation(A1,b1,1.8,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])
    
pp.plot(epsilon,NiterR1,label = 'Omega = 1')
pp.plot(epsilon,NiterR2,label = 'Omega = 1.2')
pp.plot(epsilon,NiterR3,label = 'Omega = 1.4')
pp.plot(epsilon,NiterR4,label = 'Omega = 1.6')
pp.plot(epsilon,NiterR5,label = 'Omega = 1.8')

print(NiterR1[-1])
pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()

"""
A2,b2,x0 = Matriceb(100)
epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []


for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRelaxation(A2,b2,1,x0,eps,Nitermax)
    R2 = MIRelaxation(A2,b2,1.2,x0,eps,Nitermax)
    R3 = MIRelaxation(A2,b2,1.4,x0,eps,Nitermax)
    R4 = MIRelaxation(A2,b2,1.6,x0,eps,Nitermax)
    R5 = MIRelaxation(A2,b2,1.8,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])
    
pp.plot(epsilon,NiterR1,label = 'Omega = 1')
pp.plot(epsilon,NiterR2,label = 'Omega = 1.2')
pp.plot(epsilon,NiterR3,label = 'Omega = 1.4')
pp.plot(epsilon,NiterR4,label = 'Omega = 1.6')
pp.plot(epsilon,NiterR5,label = 'Omega = 1.8')

pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()

A1,b1,x0 = Matrice(100)

epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []


for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRelaxation(A1,b1,1,x0,eps,Nitermax)
    R2 = MIRelaxation(A1,b1,0.2,x0,eps,Nitermax)
    R3 = MIRelaxation(A1,b1,0.4,x0,eps,Nitermax)
    R4 = MIRelaxation(A1,b1,0.6,x0,eps,Nitermax)
    R5 = MIRelaxation(A1,b1,0.8,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])
    
pp.plot(epsilon,NiterR1,label = 'Omega = 1')
pp.plot(epsilon,NiterR2,label = 'Omega = 0.2')
pp.plot(epsilon,NiterR3,label = 'Omega = 0.4')
pp.plot(epsilon,NiterR4,label = 'Omega = 0.6')
pp.plot(epsilon,NiterR5,label = 'Omega = 0.8')


pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()


A2,b2,x0 = Matriceb(100)
epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []


for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRelaxation(A2,b2,0.7,x0,eps,Nitermax)
    R2 = MIRelaxation(A2,b2,0.75,x0,eps,Nitermax)
    R3 = MIRelaxation(A2,b2,0.8,x0,eps,Nitermax)
    R4 = MIRelaxation(A2,b2,0.85,x0,eps,Nitermax)
    R5 = MIRelaxation(A2,b2,0.9,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])
    
pp.plot(epsilon,NiterR1,label = 'Omega = 0.7')
pp.plot(epsilon,NiterR2,label = 'Omega = 0.75')
pp.plot(epsilon,NiterR3,label = 'Omega = 0.8')
pp.plot(epsilon,NiterR4,label = 'Omega = 0.85')
pp.plot(epsilon,NiterR5,label = 'Omega = 0.9')

pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()


A1,b1,x0 = Matrice(100)

epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []    
for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRichardson(A1,b1,0.2,x0,eps,Nitermax)
    R2 = MIRichardson(A1,b1,0.25,x0,eps,Nitermax)
    R3 = MIRichardson(A1,b1,0.3,x0,eps,Nitermax)
    R4 = MIRichardson(A1,b1,0.35,x0,eps,Nitermax)
    R5 = MIRichardson(A1,b1,0.4,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])
    
pp.plot(epsilon,NiterR1,label = 'Omega = 0.2')
pp.plot(epsilon,NiterR2,label = 'Omega = 0.25')
pp.plot(epsilon,NiterR3,label = 'Omega = 0.3')
pp.plot(epsilon,NiterR4,label = 'Omega = 0.35')
pp.plot(epsilon,NiterR5,label = 'Omega = 0.4')


pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()


A2,b2,x0 = Matriceb(100)
epsilon = []
NiterR1 = []
NiterR2 = []
NiterR3 = []
NiterR4 = []
NiterR5 = []


for i in range(0,12,1):
    eps = 10**(-i)
    epsilon.append(eps)
    
    R1 = MIRichardson(A2,b2,0.1,x0,eps,Nitermax)
    R2 = MIRichardson(A2,b2,0.15,x0,eps,Nitermax)
    R3 = MIRichardson(A2,b2,0.2,x0,eps,Nitermax)
    R4 = MIRichardson(A2,b2,0.25,x0,eps,Nitermax)
    R5 = MIRichardson(A2,b2,0.3,x0,eps,Nitermax)
    NiterR1.append(R1[1])
    NiterR2.append(R2[1])
    NiterR3.append(R3[1])
    NiterR4.append(R4[1])
    NiterR5.append(R5[1])

pp.plot(epsilon,NiterR1,label = 'Omega = 0.1')
pp.plot(epsilon,NiterR2,label = 'Omega = 0.15')
pp.plot(epsilon,NiterR3,label = 'Omega = 0.2')
pp.plot(epsilon,NiterR4,label = 'Omega = 0.25')
pp.plot(epsilon,NiterR5,label = 'Omega = 0.3')
pp.xscale('log')
pp.ylabel('Ninter')
pp.xlabel('Epsilon')
pp.legend()
pp.show()
"""