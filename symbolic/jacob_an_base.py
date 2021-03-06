# -*- coding: utf-8 -*-
import numpy as np 
import sympy as sp 
import sympy.matrices as Matrix
from funciones import *


q11, q12, q13 = sp.symbols("q11 q12 q13")
q21, q22, q23 = sp.symbols("q21 q22 q23")
q31, q32, q33 = sp.symbols("q31 q32 q33")
q41, q42, q43 = sp.symbols("q41 q42 q43")

l1, l2, l3, l4 = sp.symbols("l1 l2 l3 l4")
d1 = sp.symbols("d1")

# Matrices de transformación homogénea i con respecto a i-1
# Pata 1
T1_01 = sTdh(l1, q11, 0, sp.pi/2)
T1_12 = sTdh(l2, q12, -l3, 0)
T1_23 = sTdh(0, q13, -l4, 0)

# Pata 2
T2_01 = sTdh(l1, sp.pi+q21, 0, sp.pi/2)
T2_12 = sTdh(l2, q22, l3, 0)
T2_23 = sTdh(0, q23, l4, 0)

# Pata 3
T3_01 = sTdh(l1, sp.pi+q31, 0, sp.pi/2)
T3_12 = sTdh(l2, q32, l3, 0)
T3_23 = sTdh(0, q33, l4, 0)

# Pata 4
T4_01 = sTdh(l1, q41, 0, sp.pi/2)
T4_12 = sTdh(l2, q42, -l3, 0)
T4_23 = sTdh(0, q43, -l4, 0)

# Patas Respecto a la base

# Pata 1
T1_03 = sp.simplify(T1_01*T1_12*T1_23)
T1_B0 = sTroty(-sp.pi/2)*sTrotx(sp.pi)*sTrasl(0,-d1,0)
T1_B3 = sp.simplify(T1_B0*T1_03)

# Pata 2
T2_03 = sp.simplify(T2_01*T2_12*T2_23)
T2_B0 = sTroty(-sp.pi/2)*sTrasl(0,d1,0)
T2_B3 = sp.simplify(T2_B0*T2_03)

# Pata 3
T3_03 = sp.simplify(T2_01*T3_12*T3_23)
T3_B0 = sTroty(-sp.pi/2)*sTrotx(sp.pi)*sTrasl(0,d1,0)
T3_B3 = sp.simplify(T3_B0*T3_03)

# Pata 4
T4_03 = sp.simplify(T4_01*T4_12*T4_23)
T4_B0 = sTroty(-sp.pi/2)*sTrasl(0,-d1,0)
T4_B3 = sp.simplify(T4_B0*T4_03)

# Jacobiano analitico
sp.init_printing()
#Posicion de la matriz de la Matriz de Transformacion
D1 =T1_B3[0:3,3]
D2 =T2_B3[0:3,3]
D3 =T3_B3[0:3,3]
D4 =T4_B3[0:3,3]

#Orientacion de la Matriz de Transformacion en cuaterniones
P1 =  sp.simplify(simquaterion(T1_B3[0:3,0:3]))
P2 =  sp.simplify(simquaterion(T2_B3[0:3,0:3]))
P3 =  sp.simplify(simquaterion(T3_B3[0:3,0:3]))
P4 =  sp.simplify(simquaterion(T4_B3[0:3,0:3]))

J1 = sp.simplify(sp.Matrix.vstack(T1_B3[0:3,3],P1).jacobian(sp.Matrix([q11, q12, q13])))
J2 = sp.simplify(sp.Matrix.vstack(T2_B3[0:3,3],P2).jacobian(sp.Matrix([q21, q22, q23])))
J3 = sp.simplify(sp.Matrix.vstack(T3_B3[0:3,3],P3).jacobian(sp.Matrix([q31, q32, q33])))
J4 = sp.simplify(sp.Matrix.vstack(T4_B3[0:3,3],P4).jacobian(sp.Matrix([q41, q42, q43])))


print("Jacobiano analitico de la Pata 1:"); sp.pprint(J1)
print("Jacobiano analitico de la Pata 2:"); sp.pprint(J2)
print("Jacobiano analitico de la P4ata 3:"); sp.pprint(J3)
print("Jacobiano analitico de la Pata 4:"); sp.pprint(J4)

