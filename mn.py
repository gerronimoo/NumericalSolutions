#Метод Ньютона
import numpy as np

def F(x):
	y = np.zeros_like(x)
	y[0] = (3 + 2*x[0])*x[0] - 2*x[1] - 3
	y[1:-1] = (3 + 2*x[1:-1])*x[1:-1] - x[:-2] - 2*x[2:] -2
	y[-1] = (3 + 2*x[-1])*x[-1] - x[-2] - 4
	return y

def J(x):
	n = len(x)
	ja = np.zeros((n, n))
	ja[0, 0] = 3 + 4*x[0]
	ja[0, 1] = -2
	
	for i in range(n-1):
		ja[i, i-1] = -1
		ja[i, i] = 3 + 4*x[i]
		ja[i, i+1] = -2

	ja[-1, -2] = -1
	ja[-1, -1] = 3 + 4*x[-1]
	return ja

def Newton(F, J, pre):
	N = len(pre)
	delta = np.ones(N)
	x = np.array(pre, float)
	acc = 0.001
	k = 0
	
	while max(abs(delta)) > acc and k < 100:
		delta = np.linalg.solve(J(x), -F(x))
		x = x + delta
		k += 1
	
	return x, k

n = 10
pre = 3*np.ones(n)
sol, its = Newton(F, J, pre)

if its > 0:
	print("x = {}".format(sol))
else:
	print("Решение не найдено!")