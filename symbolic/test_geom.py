# -- coding: utf-8 --
import numpy as np 
import sympy as sp 
import sympy.matrices as Matrix
from funciones import *


# Cinemática diferencial de base flotante-------------------------
sp.init_printing()
pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz = sp.symbols('pb_x pb_y pb_z omeg_b eps_bx eps_by eps_bz')
q11, q12, q13 = sp.symbols("q11 q12 q13")
l1, l2, l3, l4 = sp.symbols("l1 l2 l3 l4")
d1 = sp.symbols("d1")

cos = sp.cos 
sin = sp.sin


# Parametros de Denavit-Hartenberg

T1_B0 = sTroty(-sp.pi/2)*sTrotx(sp.pi)*sTrasl(0,-d1,0)

DH1 = [[l1,   q11, 0, sp.pi/2,'r'],
    [ l2, q12,  -l3,       0,'r'],
    [  0,q13, -l4, 0,'r']]

# Jacobianos Geometricos de la pata con respecto al sistema 0
[T1,JG1] = jacob_g(DH1)
T1_03 = T1[2]
T1_B0 = sTroty(-sp.pi/2)*sTrotx(sp.pi)*sTrasl(0,-d1,0)
T1_B3 = sp.simplify(T1_B0*T1_03)


# Jacobiano Geométrico de la pata con respecto a la base
RB_10 = T1_B0[0:3,0:3]
XB_10 = sp.Matrix.vstack(sp.Matrix.hstack(RB_10,sp.zeros(3,3)),sp.Matrix.hstack(sp.zeros(3,3),RB_10)) 
JB_1 = XB_10*JG1


# Traslación del sistema inercial al base (Transformación)
TT_IB = sTrasl(pb_x,pb_y,pb_z) 
quater_B = Matrix([[omeg_b],[eps_bx],[eps_by],[eps_bz]])

# Rotación del sistema inercial al base (Transformación)
R_IB = symrquater(quater_B) 


# Transformada homogénea de la base con respecto al sistema inercial
T_IB = Matrix([[R_IB[0,0],R_IB[0,1],R_IB[0,2],pb_x],
			   [R_IB[1,0],R_IB[1,1],R_IB[1,2],pb_y],
			   [R_IB[2,0],R_IB[2,1],R_IB[2,2],pb_z],
			   [0,0,0,1]])

# Matriz E0 de cuaterniones a velocidad angular
TQ = 2*Matrix([[-eps_bx,   omeg_b,  -eps_bz,    eps_by],
               [-eps_by,   eps_bz,   omeg_b,   -eps_bx],
               [-eps_bz,  -eps_by,   eps_bx,    omeg_b]])


# Actualización de sistemas de referncia con respecto a la base
T1_I3 = T_IB*T1_B3

# print("De la pata a la base")
# sp.pprint(sp.simplify(T1_B3))
# print("De la base a la referencia")
# sp.pprint(sp.simplify(T_IB))


# sp.pprint(T1_I3.subs([(l1,125),(l2,25),(l3,103),(l4,104),(d1,75),(q11,0),(q12,0),(q13,0),(eps_bx,0),(eps_by,0),(eps_bz,0),(omeg_b,1),(pb_x,0),(pb_y,0),(pb_z,207)]))


sk_d = symskew(T_IB[0:3,3] - T1_I3[0:3,3]) # Skew Matrix de la diferencia de pb - pi
vel_jb = sk_d*TQ
eye = Matrix([[1,0,0],
              [0,1,0],
              [0,0,1]])
zeros = Matrix([[0,0,0],
                [0,0,0],
                [0,0,0]])
Jb_base = Matrix([[eye,vel_jb],
             [zeros,TQ]])

# # sp.pprint(JB_1)
# # print("\n")
# # sp.pprint(JB_1[0:3,0:3])
# # print("\n")
# # sp.pprint(JB_1[3:6,0:3])
Jg_1 = Matrix.vstack(R_IB*JB_1[0:3,0:3],R_IB*JB_1[3:6,0:3])

Jg_1 = Matrix.hstack(Jb_base,Jg_1)

sp.pprint(Jg_1[0:3,:])
print("\n")
print(Jg_1[0:3,:])

print()
