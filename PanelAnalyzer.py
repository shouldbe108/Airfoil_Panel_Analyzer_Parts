from myFEMlib import*
from myPreProcessor import create_mesh
from PanelAssembler import PanelAssembler
import math
from scipy import linalg
import matplotlib.pyplot as plt
import os
from Solver import solver
import numpy
#/////////////////////////////////////////////////////////////////////////////
#Compute the tangential velscity using L_ij,K_ij,sigma and gamma, tangential vector and x-dir vector e1. Returns tang vel
def tangential_velocity(mesh,x,L_ij,K_ij,V_inf):
    e1 = [1,0]
    lksum=numpy.zeros(mesh.__len__())
    vt=numpy.zeros(mesh.__len__()) #initialize empty matrix vtangential
    for i in range(mesh.__len__()):
        for k in range(mesh.__len__()):
            lksum[i]=lksum[i]+x[k]*L_ij[i][k]+x[-1]*K_ij[i][k]

    for i in range(mesh.__len__()):
        tang_i = mesh[i].get_coord_tangent() #Get tangent vector for i th panel
        vt[i]= V_inf*numpy.dot(tang_i,e1)+lksum[i]
    return vt
#_______________________
#unittest 4.1
def circulation_test_IV1(mesh, x, vt):

    LHS_sum=0
    hj=0
    for i in range(mesh.__len__()):
        h_j = mesh[i].get_length()
        LHS_sum=LHS_sum+h_j*vt[i]
    hj=0 #zero sum length var
    for i in range(mesh.__len__()):
        h_j = mesh[i].get_length()
        hj=hj+h_j#Sum of lengths of panels
    RHS_sum=-1*x[-1]*hj
    diff=LHS_sum-RHS_sum
    return diff

#________________________
#compute the C_pressure using tang vel and V_inf. Returns array containing Cp values over midpoint of each panel.

def c_pressure(mesh,vt,V_inf):

    cp=numpy.zeros(mesh.__len__()) #initialize empty matrix pressure coeff

    for i in range(mesh.__len__()):
        cp[i]=1-((vt[i]/V_inf))**2
    return cp
#_________________________
#Compute the chord length for the downforce calculation and then use panel length and V_inf to return c_downforce
def c_downforce(mesh,x,V_inf):
    chord = numpy.linalg.norm(mesh[0].get_coord1()-mesh[int(mesh.__len__()/2)].get_coord1())
    hj=0 #zero sum length var
    for i in range(mesh.__len__()):
        h_j = mesh[i].get_length()
        hj=hj+h_j#Sum of lengths of panels
    cf=((2*x[-1])/(chord*V_inf))*hj  #downforce
    return cf
#_________________________
def hj(mesh):
    chord = numpy.linalg.norm(mesh[0].get_coord1()-mesh[int(mesh.__len__()/2)].get_coord1())
    hj=0 #zero sum length var
    for i in range(mesh.__len__()):
        h_j = mesh[i].get_length()
        hj=hj+h_j#Sum of lengths of panels
    return hj,
#_________________________

def cp_plotter(mesh,cp):
    scaling=1

    n = mesh.__len__()
    midpoints=numpy.zeros((n,2))
    coord1=numpy.zeros((n,2))
    for i in range(0,n):
        midpoints[i,:]=mesh[i].get_coord_midpoint()
        coord1[i,:]=mesh[i].get_coord1()

        vec_basis = mesh[i].get_coord_midpoint()

        normal_vec =mesh[i].get_coord_normal()

        tangent_vec =mesh[i].get_coord_tangent()

#        plt.arrow(vec_basis[0] , vec_basis[1] , cp[i]*scaling*normal_vec[0], cp[i]*scaling*normal_vec[1], head_width=0.01, head_length=0.01, fc='k', ec='g')
#        plt.arrow(vec_basis[0] , vec_basis[1] , scaling*tangent_vec[0], scaling*tangent_vec[1], head_width=scaling*0.2, head_length=scaling*0.2, fc='k', ec='b')

    coord1 = numpy.vstack([coord1, coord1[0,:]])#add first element at the end
    plt.plot(midpoints[:,0], midpoints[:,1],'ro')
    plt.plot(coord1[:,0], coord1[:,1],'k--')
    plt.plot(coord1[:,0], coord1[:,1],'ko')
    plt.plot(midpoints[:,0],cp)
    plt.ylabel('c_p')
    plt.axes().set_aspect('equal')
    plt.axes().set_aspect('equal', 'datalim')
#
    if os.path.isfile('mesh.PNG'):
        os.remove('mesh.PNG')
    plt.savefig('mesh.PNG')
    plt.show()
#
#    l=numpy.linspace(0,len(n))
#    plt.scatter(l,cp[0:(len(cp))//2])
#    rev=l
#    plt.scatter(l,cp[0:len(cp)//2])
#    plt.scatter(rev,cp[len(cp)//2:len(cp)])

#/////////////////////////////////////////////////////////////////////////////
#unit test 4.2
#The test is predefined for NACA0002 airfoil AOA 3degress
def NACA_Airfoil_test1(npanels):

    XYZZ=str("0002")
    X = int(XYZZ[0])
    Y = int(XYZZ[1])
    ZZ = int(XYZZ[2:])

    c=1.0
    t = c*ZZ/100
    p = Y/10
    m = X/100

    npanels_per_side = (npanels/2)
    nnodes_per_side = npanels_per_side +1

    nns=int(nnodes_per_side)
    nps=int(npanels_per_side)

    coord1 = numpy.zeros((2*nns-2,2))
    coord2 = numpy.zeros((2*nns-2,2))
    mesh = []
    N = len(coord1)

    node_density = 50

    j = -1
    for i in range(nps,-1,-1):
        j += 1

        x = c * (math.pow(node_density,(i/nps))-1)/(node_density - 1)
        if p==0:
            dyc_dx = 0
            yc = 0
        elif x >= 0 and x <= p*c:
            dyc_dx = 2*m*(c*p-x)/(c*p**2)
            yc = m*x*(2*p-x/c)/(p**2)
        elif x > p*c and x <= c:
            dyc_dx = -2*m*(x-c*p)/(c*(p-1)**2)# check this one
            yc = m*(c-x)*(1+x/c -2*p)/((1-p)**2)
        else:
            sys.exit("Something is wrong with the defenition of x.")
#
        yt = 5*t*(0.2969*math.sqrt(x/c) - 0.1260*x/c - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1036*(x/c)**4)
        theta = math.atan2(dyc_dx,1)

        x_upper = x - yt*math.sin(theta)
        y_upper = yc + yt*math.cos(theta)
        if j!=0:
            coord1[int(N/2)+i,0]= x_upper
            coord1[int(N/2)+i,1]= y_upper


        x_lower = x + yt*math.sin(theta)
        y_lower = yc - yt*math.cos(theta)

        coord1[j,0]= x_lower
        coord1[j,1]= y_lower

    coord1 = rotation_translation_to_mesh_test1(coord1)

    create_mesh(j,coord1,mesh)
    return mesh
#_________________________

def rotation_translation_to_mesh_test1(coord1):

    rotation=math.radians(float(3))
    translation_x=float(0)
    translation_y=float(0)

    for i in range(len(coord1)):
        radius = numpy.linalg.norm(coord1[i])
        init_angle = math.atan2(coord1[i,1],coord1[i,0])
#        print("radius = ", radius, "angle = ", init_angle)
        new_angle = init_angle + rotation

        coord1[i,0] = math.cos(new_angle)*radius + translation_x
        coord1[i,1] = math.sin(new_angle)*radius + translation_y

    return coord1
#_________________________

def NACA_Airfoil_rotplot(npanels,rotation):

    XYZZ=str("0002")
    X = int(XYZZ[0])
    Y = int(XYZZ[1])
    ZZ = int(XYZZ[2:])
    aoa=rotation
    c=1.0
    t = c*ZZ/100
    p = Y/10
    m = X/100

    npanels_per_side = (npanels/2)
    nnodes_per_side = npanels_per_side +1

    nns=int(nnodes_per_side)
    nps=int(npanels_per_side)

    coord1 = numpy.zeros((2*nns-2,2))
    coord2 = numpy.zeros((2*nns-2,2))
    mesh = []
    N = len(coord1)

    node_density = 50

    j = -1
    for i in range(nps,-1,-1):
        j += 1

        x = c * (math.pow(node_density,(i/nps))-1)/(node_density - 1)
        if p==0:
            dyc_dx = 0
            yc = 0
        elif x >= 0 and x <= p*c:
            dyc_dx = 2*m*(c*p-x)/(c*p**2)
            yc = m*x*(2*p-x/c)/(p**2)
        elif x > p*c and x <= c:
            dyc_dx = -2*m*(x-c*p)/(c*(p-1)**2)# check this one
            yc = m*(c-x)*(1+x/c -2*p)/((1-p)**2)
        else:
            sys.exit("Something is wrong with the defenition of x.")
#
        yt = 5*t*(0.2969*math.sqrt(x/c) - 0.1260*x/c - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1036*(x/c)**4)
        theta = math.atan2(dyc_dx,1)

        x_upper = x - yt*math.sin(theta)
        y_upper = yc + yt*math.cos(theta)
        if j!=0:
            coord1[int(N/2)+i,0]= x_upper
            coord1[int(N/2)+i,1]= y_upper


        x_lower = x + yt*math.sin(theta)
        y_lower = yc - yt*math.cos(theta)

        coord1[j,0]= x_lower
        coord1[j,1]= y_lower

    coord1 = rotation_translation_to_mesh_plot(coord1,aoa)

    create_mesh(j,coord1,mesh)
    return mesh
#_________________________

def rotation_translation_to_mesh_plot(coord1,rotation):

    rotation=math.radians(float(rotation))
    translation_x=float(0)
    translation_y=float(0)

    for i in range(len(coord1)):
        radius = numpy.linalg.norm(coord1[i])
        init_angle = math.atan2(coord1[i,1],coord1[i,0])
#        print("radius = ", radius, "angle = ", init_angle)
        new_angle = init_angle + rotation

        coord1[i,0] = math.cos(new_angle)*radius + translation_x
        coord1[i,1] = math.sin(new_angle)*radius + translation_y

    return coord1
#_________________________

#tutorial for 4.2 aoa Cf vs AOA for NACA airfoil.
def aoa_plot():
    aoa=numpy.linspace(-10,10,20)
    cf=numpy.zeros(len(aoa))
    for i in range (len(aoa)):
        V_inf = 1
        mesh=NACA_Airfoil_rotplot(50,aoa[[i]])
        #2 Panel assembler
        A, b, L_ij, K_ij, I_ij, J_ij, L_j, K_j = PanelAssembler(mesh,V_inf)
        #3 The panel method system solver
        x=solver(A,b)
        cf[i]=c_downforce(mesh,x,V_inf)
    plt.scatter(aoa,cf)
    plt.ylabel('c_f')
    plt.xlabel('alpha')
    return cf
#_________________________

def aoa_test_IV2():
    V_inf = 1
    mesh=NACA_Airfoil_test1(50)
    #2 Panel assembler
    A, b, L_ij, K_ij, I_ij, J_ij, L_j, K_j = PanelAssembler(mesh,V_inf)
    #3 The panel method system solver
    x=solver(A,b)
    cf=c_downforce(mesh,x,V_inf)
    aoa = math.radians(float(3))
    diff=cf-2*numpy.pi*aoa
    return diff,x
#/////////////////////////////////////////////////////////////////////////////
#unit test 4.3
#The test is predefined for NACA0008 airfoil AOA 5degress
def NACA_Airfoil_test(npanels):

    XYZZ=str("0008")
    X = int(XYZZ[0])
    Y = int(XYZZ[1])
    ZZ = int(XYZZ[2:])

    c=1.0
    t = c*ZZ/100
    p = Y/10
    m = X/100

    npanels_per_side = (npanels/2)
    nnodes_per_side = npanels_per_side +1

    nns=int(nnodes_per_side)
    nps=int(npanels_per_side)

    coord1 = numpy.zeros((2*nns-2,2))
    coord2 = numpy.zeros((2*nns-2,2))
    mesh = []
    N = len(coord1)

    node_density = 200

    j = -1
    for i in range(nps,-1,-1):
        j += 1

        x = c * (math.pow(node_density,(i/nps))-1)/(node_density - 1)
        if p==0:
            dyc_dx = 0
            yc = 0
        elif x >= 0 and x <= p*c:
            dyc_dx = 2*m*(c*p-x)/(c*p**2)
            yc = m*x*(2*p-x/c)/(p**2)
        elif x > p*c and x <= c:
            dyc_dx = -2*m*(x-c*p)/(c*(p-1)**2)# check this one
            yc = m*(c-x)*(1+x/c -2*p)/((1-p)**2)
        else:
            sys.exit("Something is wrong with the defenition of x.")
#
        yt = 5*t*(0.2969*math.sqrt(x/c) - 0.1260*x/c - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1036*(x/c)**4)
        theta = math.atan2(dyc_dx,1)

        x_upper = x - yt*math.sin(theta)
        y_upper = yc + yt*math.cos(theta)
        if j!=0:
            coord1[int(N/2)+i,0]= x_upper
            coord1[int(N/2)+i,1]= y_upper


        x_lower = x + yt*math.sin(theta)
        y_lower = yc - yt*math.cos(theta)

        coord1[j,0]= x_lower
        coord1[j,1]= y_lower

    coord1 = rotation_translation_to_mesh_test(coord1)

    create_mesh(j,coord1,mesh)
    return mesh
#_________________________

def rotation_translation_to_mesh_test(coord1):

    rotation=math.radians(float(5))
    translation_x=float(0)
    translation_y=float(0)

    for i in range(len(coord1)):
        radius = numpy.linalg.norm(coord1[i])
        init_angle = math.atan2(coord1[i,1],coord1[i,0])
#        print("radius = ", radius, "angle = ", init_angle)
        new_angle = init_angle + rotation

        coord1[i,0] = math.cos(new_angle)*radius + translation_x
        coord1[i,1] = math.sin(new_angle)*radius + translation_y

    return coord1
#_________________________

def refiner_test_IV3():
    npanels = [10,200,202,204]
    cf=numpy.zeros(4);

    for k in range (0,4):
        V_inf=1
        mesh=NACA_Airfoil_test(npanels[k])
        #2 Panel assembler
        A, b, L_ij, K_ij, I_ij, J_ij, L_j, K_j = PanelAssembler(mesh,V_inf)
        #3 The panel method system solver
        x=solver(A,b)
        #4Panel Analyzer
        cf[k]=c_downforce(mesh,x,V_inf)

    return cf
