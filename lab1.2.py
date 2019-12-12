﻿#Метод Гаусса с частичным выбором ведущего элемента
import numpy as np

A=np.array([[3,17,10],[2,4,-2],[6,18,-12]],float)
x=np.array([1,3,4])
b=np.dot(A,x)
qo=np.zeros(len(A)-1,float)

def gss(x):
	x = np.array(x,float)
	return x[1:]/x[0]


def gss_a(C, t):
	C = np.array(C, float)
	
	t = np.array([[t[i]] for i in range(len(t))], float)
	C[1:,:] = C[1:,:] - t*C[0,:]
	return C


def LU(A):
	A=np.array(A,float)
	
	for k in range(0,len(A)-1,1):
		qo[k]=list(A[k:,k]).index(max(A[k:,k]))
		ch = np.array(A[k,k:len(A)],float)
		A[k,k:len(A)] = A[int(qo[k]),k:len(A)]
		A[int(qo[k]),k:len(A)] = ch
		
		if A[k,k]!=0:
			t=gss(A[k:,k])
			
		A[k+1:,k]=t
		A[k:,k+1:]=gss_a(A[k:,k+1:],t)
	return A

def solve(A,b):

    C=LU(A)
    print(qo)
    ch1=np.zeros(len(A),float)
    M=[]
    E=[]
    
    for i in range (len(A)+1):
        M.append(np.identity(len(A),float))
        E.append(np.identity(len(A),float))
    
    for j in range(0,len(A)-1,1):
        ch1[:]=E[j][j,:]
        E[j][j,:]=E[j][int(qo[j]),:]
        E[j][int(qo[j]),:]=ch1[:]
    
    for k in range(1,len(A)):
        M[k-1][k:,k-1]=-np.dot(E[k],C)[k:,k-1]
   
    U=np.array(A,float)
    y=np.array(b,float)
    
    for t in range(0,len(A)-1):
        U=np.dot(M[t],np.dot(E[t],U))
        y=np.dot(M[t],np.dot(E[t],y))
    
    print(U)
    print(y)
    
    r=np.zeros(len(A),float)
    
    for i in range (len(A),0,-1):
        temp=0
        for j in range (i,len(A),1):
            temp=temp+U[i-1,j]*r[j]
        r[i-1]=(y[i-1]-temp)/U[i-1,i-1]
    print(r)
    
solve(A,b)