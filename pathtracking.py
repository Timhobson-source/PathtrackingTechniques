# Path tracking attempt.
from cmath import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
import newton as newt

# solutions we know:
def g(z):
    return z**2-1 #roots 1,-1.

def dg(z):
    return 2*z

#solutions we "don't know":
def f(z):
    return 5*z**2 +z -37

def df(z):
    return 10*z +1

#calculating the solutions to f=0 to error 10e-4 by newton's method:
#space = np.linspace(-5,15,num=50)
#plt.plot(space,[f(s) for s in space],color="red")
#plt.axhline(0)
#plt.show() #one real solution

fsol_1=newt.iterate_newton(f,df,x_0=5,error=0.0001)
fsol_2=newt.iterate_newton(f,df,x_0=13,error=0.0001) #solutions to f

# Now for the actual path tracking.

#Need to define our H as:

gamma=0.44+0.77j

def H(z,t):
    return gamma*t*g(z)+(1-t)*f(z)

def JH(z,t):
    return gamma*t*dg(z)+(1-t)*df(z)

def Dt_H(z,t):
    return gamma*g(z)-f(z)

# Have to solve IVP: p'(t) = -(JH(p(t),t)**(-1))*Dt_H(p(t),t)=F(p(t),t), p(1)= 1, where 1 is a root of gamma. We want p(0).

def euler(x,t,step,func):
    return x-step*func(x,t)

def F(x,t):
    return -(JH(x,t)**(-1))*Dt_H(x,t)

def LAonestep(p,t,step):
    z_0=euler(p,t,step,F)
    z_1=z_0-H(z_0,t)/JH(z_0,t+step)
    p_1=z_1-H(z_1,t+step)/JH(z_1,t+2*step)
    return p_1

p=1
#q=-1
t=1
N=500
for i in range(1,N+1):
    #r=LAonestep(q,t,-1/N)
    s=LAonestep(p,t,-1/N)
    p=s
    #q=r
    t-=1/N
print(p)#,q)
print(fsol_1,fsol_2)
# It works!
