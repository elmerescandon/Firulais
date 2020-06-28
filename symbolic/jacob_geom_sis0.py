from funciones import *
import sympy as sp 
import numpy as np


sp.init_printing()

q11, q12, q13 = sp.symbols("q11 q12 q13")
q21, q22, q23 = sp.symbols("q21 q22 q23")
q31, q32, q33 = sp.symbols("q31 q32 q33")
q41, q42, q43 = sp.symbols("q41 q42 q43")

l1, l2, l3, l4 = sp.symbols("l1 l2 l3 l4")
d1 = sp.symbols("d1")
#Parametros de Denavit-Hartenberg
#Pata 
DH1 = [[l1,   q11, 0, sp.pi/2,'r'],
    [ l2, q12,  -l3,       0,'r'],
    [  0,q13, -l4, 0,'r']]

#Pata 2
DH2 = [[l1, sp.pi+q21, 0, sp.pi/2,'r'],
      [ l2, q22,  l3,       0,'r'],
      [ 0,q23, l4, 0,'r']]

#Pata 3
DH3 = [[l1, sp.pi+q31, 0, sp.pi/2,'r'],
      [ l2, q32,  l3,       0,'r'],
      [ 0,q33, l4, 0,'r']]

#Pata 4
DH4 = [[l1,   q41, 0, sp.pi/2,'r'],
      [ l2, q42,  -l3,       0,'r'],
      [  0,q43, -l4, 0,'r']]

#Jacobianos Geometricos
[T1,JG1] = jacob_g(DH1)
[T2,JG2] = jacob_g(DH2)
[T3,JG3] = jacob_g(DH3)
[T4,JG4] = jacob_g(DH4)

# Imprimir los valores

# print("Jacobiano geometrico de la pata 1 respecto a su sistema 0: "); 
# sp.pprint(JG1)
# print("Jacobiano geometrico de la pata 2 respecto a su sistema 0: "); 
# sp.pprint(JG2)
# print("Jacobiano geometrico de la pata 3 respecto a su sistema 0: "); 
# sp.pprint(JG3)
# print("Jacobiano geometrico de la pata 4 respecto a su sistema 0: "); 
# sp.pprint(JG4)

