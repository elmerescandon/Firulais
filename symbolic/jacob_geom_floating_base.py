# -*- coding: utf-8 -*-
import numpy as np 
import sympy as sp 
import sympy.matrices as Matrix
from funciones import *
from fk_dh import *
from jacob_geom_base import *

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
T1_I3 = T_IB*T1_B3
T2_I3 = T_IB*T2_B3
T3_I3 = T_IB*T3_B3
T4_I3 = T_IB*T4_B3

# Skew Matrix de la diferencia de pb - pi
sk_d = symskew(T_IB[0:3,3]-T1_I3[0:3,3]) 
vel_jb = sk_d*TQ
eye = Matrix([[1,0,0],
              [0,1,0],
              [0,0,1]])
zeros = Matrix([[0,0,0],
                [0,0,0],
                [0,0,0]])

# Obtención del Jacobiano geométrico de la base
Jb_base = Matrix([[eye,vel_jb],
             [zeros,TQ]])
Jg_1 = Matrix.vstack(R_IB*JB_1[0:3,0:3],R_IB*JB_1[3:6,0:3])
Jg_1 = sp.simplify(Matrix.hstack(Jb_base,Jg_1))

sp.pprint(Jg_1)