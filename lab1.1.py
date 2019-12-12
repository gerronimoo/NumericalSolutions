# Решение системы линейных уравнений с трехдиагональной матрицей
import numpy as np

A=np.array([[160,2,0,0],[6,185,5,0],[0,3,193,11],[0,0,8,134]])
b=np.array([10,22,42,72])


def solve(A,b):
	d=np.array([A[i,i] for i in range(len(A))],float)
	d1=np.array([A[i,i+1] for i in range(len(A)-1)],float)
	d2=np.array([A[i+1,i] for i in range(len(A)-1)],float)
	
	P=np.zeros(len(A),float)
	P[1]=-d1[0]/d[0]
	for i in range(2,len(A)):
		P[i]=-d1[i-1]/(d[i-1]+d2[i-1]*P[i-1])
	
	Q=np.zeros(len(A),float)
	
	Q[1]=b[0]/d[0]
	for i in range(2,len(A)):
		Q[i]=(-d2[i-1]*Q[i-1]+b[i-1])/(d[i-1]+d2[i-1]*P[i-1])
	
	x=np.zeros(len(A),float)
	
	x[-1]=(-d2[-1]*Q[-1]+b[-1])/(d[-1]+d2[-1]*P[-1])
	for i in range(len(A)-2,-1,-1):
		x[i]=P[i+1]*x[i+1]+Q[i+1]
	
	return x
    
x = solve(A,b)
print(np.dot(A,x))

