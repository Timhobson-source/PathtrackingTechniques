# numerically reproducing figure 14 (PS2 QN1d)
# We have equation S=1-exp(-k_mean * S).
# f(x)=1-exp(-u*x), u=k_mean and x=s.
# f'(x)= u*exp(-u*x)
from math import *
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

def newton(g,dg,x=1):
    # g: function
    # dg: derivative of g
    # error: error of solution
    sol = x- g(x)/dg(x)
    return sol

def g(x):
    return (x-1)**2
def dg(x):
    return 2*(x-1)

def iterate_newton(g,dg,x_0=1,error=10**(-3)):
    e=1
    x=x_0
    while e>error:
        y=newton(g,dg,x)
        e=abs(y-x)
        x=y
    return x
        
        
   

