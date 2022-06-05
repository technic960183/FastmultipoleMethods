import numpy as np 
import math
import scipy.special as sp

q = [1,1,2]
rho = [1,3,2]
alpha = [0.5,2,1.3]
beta = [3.1,2,1.5]
Theta = 0.5
Phi = 1.2
r = 5
def Y(m,n,theta,phi):
    y = np.sqrt(math.factorial(n-abs(m))/math.factorial(n+abs(m)))\
        *sp.lpmn(abs(m),n,np.cos(theta))[0][-1][n]*np.cos(m*phi)
    return y

def M(m,n):
    k = len(q)
    summ=0
    for i in range (k):
        summ += q[i]*rho[i]**n*Y(-m,n,alpha[i],beta[i])
    return summ

def potential(p,r):
    tot = 0
    for n in range (p+1):
        for m in range (-n,n+1):
            tot += M(m,n)*Y(m,n,Theta,Phi)/(r**(n+1))
            print (tot)
    return tot 
#print (M(0,0))
#print (Y(0,0,Theta,Phi))
print (potential(20,5))