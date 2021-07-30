# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:55:56 2020

@author: 20181046
"""
#
#v_t=[0,0,0]
#v_inf=1
#t=[1,2,3]
#e_1=2
#N=4
#v_sum=0
#
###tangential velocity calculation
#
#def tan_vel
#for i in range ():
#    for j in range (1,N):
#        v_b=sig[j]*L[i][j]+y*K[i][j]
#        v_sum=v_sum+v_b
#    v_t[i]=v_inf*t[i]*e_1+v_sum
#
##Pressure Coeffiecient
#
#def compute_pressure_coeff
#cp=1.0-(v_t/v_inf)**2
#
##downforce coeff
#
#for i in rang(0,N-1):
#    sum_h=sum_h+h[i]
#    
#cd=(2*y)/(c*v_inf)*sum_h
#mesh[0].get_coord1()
from myFEMlib import *
import math
import numpy

def get_tangential_velocity(Panel, freestream):
    """
    Computes the tangential velocity on the surface of the Panel.
    
    Parameters
    ---------
    Panel: 1D array of Panel objects
        The source Panel.
    freestream: Freestream object
        The freestream conditions.
    """
    N = len(Panel)
    A = numpy.empty((N, N), dtype=float)
    numpy.fill_diagonal(A, 0.0)
    
    for i, p_i in enumerate(Panel):
        for j, p_j in enumerate(Panel):
            if i != j:
                A[i, j] = 0.5 / math.pi * integral(p_i.xc, p_i.yc, p_j,
                                                   -math.sin(p_i.beta),
                                                   math.cos(p_i.beta))
    
    b = freestream.u_inf * numpy.sin([freestream.alpha - panel.beta 
                                      for panel in Panel])
    
    sigma = numpy.array([panel.sigma for panel in Panel])
    
    vt = numpy.dot(A, sigma) + b
    
    for i, panel in enumerate(Panel):
        panel.vt = vt[i]
        

##########
# compute the tangential velocity at the center-point of each panel
get_tangential_velocity(Panel, freestream)

##########

def get_pressure_coefficient(Panel, freestream):
    """
    Computes the surface pressure coefficients on the Panel.
    
    Parameters
    ---------
    Panel: 1D array of Panel objects
        The source Panel.
    freestream: Freestream object
        The freestream conditions.
    """
    for panel in Panel:
        panel.cp = 1.0 - (panel.vt / freestream.u_inf)**2
##############3
        


# computes the surface pressure coefficients
get_pressure_coefficient(Panel, freestream)

