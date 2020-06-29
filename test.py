# -- coding: utf-8 --
import numpy as np
import sympy as sp
import sympy.matrices as Matrix
from funciones import *


# Cinemática diferencial de base flotante-------------------------
sp.init_printing()
pb_x,pb_y,pb_z,wb,eb_x,eb_y,eb_z = sp.symbols(r'p_{bx} p_{by} p_{bz} \omega_{b} \epsilon_{bx} \epsilon_{by} \epsilon_{bz}')

# Traslación del sistema inercial al base (Transformación)
TT_IB = sTrasl(pb_x,pb_y,pb_z)

quater_B = Matrix([[wb],[eb_x],[eb_y],[eb_z]])

# Rotación del sistema inercial al base (Transformación)
R_IB = symrquater(quater_B)

# Rotación del sistema inercial al base
TR_IB = symtransmaxtrix(R_IB)

# Transformada homogénea de la base con respecto al sistema inercial
T_IB = TT_IB*TR_IB
TQ = 2*Matrix([[-eb_x,wb,-eb_z,eb_y],
               [-eb_y,eb_z,wb,-eb_x],
               [-eb_z,-eb_y,eb_x,wb]])


# Actualización de sistemas de referncia con respecto a la base
T1_I3 = sp.simplify(T_IB*T1_B3)
T2_I3 = sp.simplify(T_IB*T2_B3)
T3_I3 = sp.simplify(T_IB*T3_B3)
T4_I3 = sp.simplify(T_IB*T4_B3)

#Posicion de la matriz de la Matriz de Transformacion
D1 =T1_I3[0:3,3]
D2 =T2_I3[0:3,3]
D3 =T3_I3[0:3,3]
D4 =T4_I3[0:3,3]

#Orientacion de la Matriz de Transformacion en cuaterniones
P1 =  sp.simplify(simquaterion(T1_I3[0:3,0:3]))
P2 =  sp.simplify(simquaterion(T2_I3[0:3,0:3]))
P3 =  sp.simplify(simquaterion(T3_I3[0:3,0:3]))
P4 =  sp.simplify(simquaterion(T4_I3[0:3,0:3]))

J1 = sp.simplify(sp.Matrix.vstack(D1,P1).jacobian(sp.Matrix([pb_x,pb_y,pb_z,wb,eb_x,eb_y,eb_z,q11, q12, q13])))
J2 = sp.simplify(sp.Matrix.vstack(D2,P2).jacobian(sp.Matrix([pb_x,pb_y,pb_z,wb,eb_x,eb_y,eb_z,q21, q22, q23])))
J3 = sp.simplify(sp.Matrix.vstack(D3,P3).jacobian(sp.Matrix([pb_x,pb_y,pb_z,wb,eb_x,eb_y,eb_z,q31, q32, q33])))
J4 = sp.simplify(sp.Matrix.vstack(D4,P4).jacobian(sp.Matrix([pb_x,pb_y,pb_z,wb,eb_x,eb_y,eb_z,q41, q42, q43])))


print("Jacobiano analitico de la Pata 1:"); sp.pprint(J1)
print("Jacobiano analitico de la Pata 2:"); sp.pprint(J2)
print("Jacobiano analitico de la Pata 3:"); sp.pprint(J3)
print("Jacobiano analitico de la Pata 4:"); sp.pprint(J4)
