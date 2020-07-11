# -- coding: utf-8 --
import numpy as np 
import sympy as sp 
import sympy.matrices as Matrix
from funciones import *


# sp.init_printing()
# pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz = sp.symbols('pb_x pb_y pb_z omeg_b eps_bx eps_by eps_bz')
# q11, q12, q13 = sp.symbols("q11 q12 q13")
# l1, l2, l3, l4 = sp.symbols("l1 l2 l3 l4")
# d1 = sp.symbols("d1")

# cos = sp.cos 
# sin = sp.sin

# fk_I1 = Matrix([[pb_x + 2*(eps_bx*eps_by - eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*(eps_bx*eps_bz + eps_by*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (-2*eps_bx**2 - 2*omeg_b**2 + 1)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13))],
# 				  [pb_y - 2*(eps_bx*eps_by + eps_bz*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b - eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))],
# 				  [pb_z + 2*(-eps_bx*eps_bz + eps_by*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b + eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (-2*eps_bz**2 - 2*omeg_b**2 + 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13))]])

# l1 = 125
# l2 = 25
# l3 = 103
# l4 = 104
# d1 = 75
# q11 = 0
# q12 = 0 
# q13 = 0 
# eps_by = 0.261 
# eps_bz = 0.521
# eps_bx = 0.272
# omeg_b = 0.766
# pb_x = 0
# pb_y = 0
# pb_z = 207

# cos = np.cos
# sin = np.sin
# sqrt = np.sqrt

# J_ANAL = np.array([[1, 0, 0, -2*eps_by*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*eps_bz*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 4*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), -4*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), 2*eps_bx*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bx*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (-2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(-l2*cos(q11) - l3*sin(q11)*cos(q12) - l4*sin(q11)*cos(q12 + q13)), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(-l3*sin(q11)*sin(q12) - l4*sin(q11)*sin(q12 + q13)) + (-2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(-l3*sin(q12)*cos(q11) - l4*sin(q12 + q13)*cos(q11)) + (l3*cos(q12) + l4*cos(q12 + q13))*(-2*eps_bx**2 - 2*omeg_b**2 + 1), -l4*(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) - l4*(-2*eps_bx*eps_bz - 2*eps_by*omeg_b)*sin(q12 + q13)*cos(q11) + l4*(-2*eps_bx**2 - 2*omeg_b**2 + 1)*cos(q12 + q13)], [0, 1, 0, 2*eps_bx*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*eps_bz*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 4*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), -2*eps_by*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 4*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_by*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), (2*eps_bx*omeg_b - 2*eps_by*eps_bz)*(-l2*cos(q11) - l3*sin(q11)*cos(q12) - l4*sin(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), (-2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(l3*cos(q12) + l4*cos(q12 + q13)) + (2*eps_bx*omeg_b - 2*eps_by*eps_bz)*(-l3*sin(q12)*cos(q11) - l4*sin(q12 + q13)*cos(q11)) + (-l3*sin(q11)*sin(q12) - l4*sin(q11)*sin(q12 + q13))*(2*eps_by**2 + 2*omeg_b**2 - 1), l4*(-2*eps_bx*eps_by - 2*eps_bz*omeg_b)*cos(q12 + q13) - l4*(2*eps_bx*omeg_b - 2*eps_by*eps_bz)*sin(q12 + q13)*cos(q11) - l4*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11)*sin(q12 + q13)], [0, 0, 1, 2*eps_bx*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + 2*eps_by*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) - 4*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bz*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), 2*eps_bz*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + 2*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), -2*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 4*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (-2*eps_bz**2 - 2*omeg_b**2 + 1)*(-l2*cos(q11) - l3*sin(q11)*cos(q12) - l4*sin(q11)*cos(q12 + q13)), (-2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l3*cos(q12) + l4*cos(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(-l3*sin(q11)*sin(q12) - l4*sin(q11)*sin(q12 + q13)) + (-l3*sin(q12)*cos(q11) - l4*sin(q12 + q13)*cos(q11))*(-2*eps_bz**2 - 2*omeg_b**2 + 1), l4*(-2*eps_bx*eps_bz + 2*eps_by*omeg_b)*cos(q12 + q13) - l4*(2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11)*sin(q12 + q13) - l4*(-2*eps_bz**2 - 2*omeg_b**2 + 1)*sin(q12 + q13)*cos(q11)]])

# print(J_ANAL[:,3:6])
# J1 = fk_I1.jacobian(sp.Matrix([pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11, q12, q13]))
# sp.pprint(J1)
# print(J1)
