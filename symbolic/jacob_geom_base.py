from funciones import *
from jacob_geom_sis0 import *
from fk_dh import *
import sympy as sp 
import numpy as np


sp.init_printing()

#Jacobiano geometrico (respecto a la base)

#Pata 1
RB_10 = T1_B0[0:3,0:3]
XB_10 = sp.Matrix.vstack(sp.Matrix.hstack(RB_10,sp.zeros(3,3)),sp.Matrix.hstack(sp.zeros(3,3),RB_10)) 
JB_1 = XB_10*JG1

#Pata 2
RB_20 = T2_B0[0:3,0:3]
XB_20 = sp.Matrix.vstack(sp.Matrix.hstack(RB_20,sp.zeros(3,3)),sp.Matrix.hstack(sp.zeros(3,3),RB_20)) 
JB_2 = XB_20*JG2


#Pata 3
RB_30 = T3_B0[0:3,0:3]
XB_30 = sp.Matrix.vstack(sp.Matrix.hstack(RB_30,sp.zeros(3,3)),sp.Matrix.hstack(sp.zeros(3,3),RB_30)) 
JB_3 = XB_30*JG3


#Pata 4
RB_40 = T4_B0[0:3,0:3]
XB_40 = sp.Matrix.vstack(sp.Matrix.hstack(RB_40,sp.zeros(3,3)),sp.Matrix.hstack(sp.zeros(3,3),RB_40)) 
JB_4 = XB_40*JG4


# Imprimir los valores

# print("Jacobiano geometrico respecto a la base de la pata 1: "); 
# sp.pprint(JB_1)

# print("Jacobiano geometrico respecto a la base de la pata 2: "); 
# sp.pprint(JB_2)

# print("Jacobiano geometrico respecto a la base de la pata 3: "); 
# sp.pprint(JB_3)

# print("Jacobiano geometrico respecto a la base de la pata 4: "); 
# sp.pprint(JB_4)