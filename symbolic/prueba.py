import numpy as np 

# Comprobación del Jacobiano Analítico

def fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13): 
	l1 = 125
	l2 = 25
	l3 = 103
	l4 = 104
	d1 = 75
	cos = np.cos 
	sin = np.sin
	sqrt = np.sqrt	

	fk_I1 = np.array([[pb_x + 2*(eps_bx*eps_by - eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*(eps_bx*eps_bz + eps_by*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (-2*eps_bx**2 - 2*omeg_b**2 + 1)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13))],
				  [pb_y - 2*(eps_bx*eps_by + eps_bz*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b - eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))],
				  [pb_z + 2*(-eps_bx*eps_bz + eps_by*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b + eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (-2*eps_bz**2 - 2*omeg_b**2 + 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13))]])
	return np.array([fk_I1[0,0],fk_I1[1,0],fk_I1[2,0]])



l1 = 125
l2 = 25
l3 = 103
l4 = 104
d1 = 75
q11 = 0
q12 = 0 
q13 = 0
eps_by = 0
eps_bz = 0.182979
eps_bx = 0.182979
omeg_b = 0.9659386
pb_x = 0
pb_y = 0
pb_z = 207



cos = np.cos 
sin = np.sin
sqrt = np.sqrt

# ============================================================================
# Cinemática directa de la pata 1 con respecto al sistema de la base
# ============================================================================ 

fk_B1 = np.array([[l1 - l3*sin(q12) - l4*sin(q12 + q13)],
		[d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)],
		[l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)]])

# ============================================================================
# Cinemática directa de la pata con respecto al sistema de referencia Inercial
# ============================================================================


# Matriz de transformación homogénea
# Matrix([[2*(-eps_bx*eps_by + eps_bz*omeg_b)*sin(q11)*cos(q12 + q13) + 2*(eps_bx*eps_bz + eps_by*omeg_b)*cos(q11)*cos(q12 + q13) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*sin(q12 + q13), 2*(eps_bx*eps_by - eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) - 2*(eps_bx*eps_bz + eps_by*omeg_b)*sin(q12 + q13)*cos(q11) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*cos(q12 + q13), 2*(eps_bx*eps_by - eps_bz*omeg_b)*cos(q11) + 2*(eps_bx*eps_bz + eps_by*omeg_b)*sin(q11), pb_x + 2*(eps_bx*eps_by - eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*(eps_bx*eps_bz + eps_by*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (-2*eps_bx**2 - 2*omeg_b**2 + 1)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13))], 
# 	[2*(eps_bx*eps_by + eps_bz*omeg_b)*sin(q12 + q13) + 2*(-eps_bx*omeg_b + eps_by*eps_bz)*cos(q11)*cos(q12 + q13) + (-2*eps_by**2 - 2*omeg_b**2 + 1)*sin(q11)*cos(q12 + q13), 2*(eps_bx*eps_by + eps_bz*omeg_b)*cos(q12 + q13) + 2*(eps_bx*omeg_b - eps_by*eps_bz)*sin(q12 + q13)*cos(q11) + (2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11)*sin(q12 + q13), 2*(-eps_bx*omeg_b + eps_by*eps_bz)*sin(q11) + (2*eps_by**2 + 2*omeg_b**2 - 1)*cos(q11), pb_y - 2*(eps_bx*eps_by + eps_bz*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b - eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))], 
# 	[2*(eps_bx*eps_bz - eps_by*omeg_b)*sin(q12 + q13) - 2*(eps_bx*omeg_b + eps_by*eps_bz)*sin(q11)*cos(q12 + q13) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*cos(q11)*cos(q12 + q13), 2*(eps_bx*eps_bz - eps_by*omeg_b)*cos(q12 + q13) + 2*(eps_bx*omeg_b + eps_by*eps_bz)*sin(q11)*sin(q12 + q13) + (-2*eps_bz**2 - 2*omeg_b**2 + 1)*sin(q12 + q13)*cos(q11), 2*(eps_bx*omeg_b + eps_by*eps_bz)*cos(q11) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q11), pb_z + 2*(-eps_bx*eps_bz + eps_by*omeg_b)*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*(eps_bx*omeg_b + eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (-2*eps_bz**2 - 2*omeg_b**2 + 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13))], 
# 	[0, 0, 0, 1]])

# Cinemática directa para la posición


# ============================================================================
# Jacobiano Geométrico para la posición de la pata 1
# ============================================================================
# Matrix([[1, 0, 0, -2*eps_by*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_by*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*eps_bz*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), -2*eps_bx*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), -(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) + (-l3*cos(q12) - l4*cos(q12 + q13))*(2*eps_bx**2 + 2*omeg_b**2 - 1), -l4*(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bx*eps_bz + 2*eps_by*omeg_b)*sin(q12 + q13)*cos(q11) - l4*(2*eps_bx**2 + 2*omeg_b**2 - 1)*cos(q12 + q13)], 
# 	[0, 1, 0, -2*eps_bx*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_by*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_by*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), (2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) - (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11), -l4*(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*cos(q12 + q13) + l4*(-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q12 + q13)*cos(q11) - l4*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11)*sin(q12 + q13)], 
# 	[0, 0, 1, -2*eps_bx*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) - 2*eps_by*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), 2*eps_bz*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_bz*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_bx*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*eps_by*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), (2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_bz**2 + 2*omeg_b**2 - 1)*cos(q11), -l4*(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*cos(q12 + q13) - l4*(2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q12 + q13)*cos(q11)], 
# 	[0, 0, 0, -2*eps_bx, 2*omeg_b, -2*eps_bz, 2*eps_by, 2*eps_bx**2 + 2*omeg_b**2 - 1, (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*cos(q11) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*sin(q11), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*cos(q11) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*sin(q11)], 
# 	[0, 0, 0, -2*eps_by, 2*eps_bz, 2*omeg_b, -2*eps_bx, 2*eps_bx*eps_by + 2*eps_bz*omeg_b, (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11) + (2*eps_by**2 + 2*omeg_b**2 - 1)*cos(q11), (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11) + (2*eps_by**2 + 2*omeg_b**2 - 1)*cos(q11)], 
# 	[0, 0, 0, -2*eps_bz, -2*eps_by, 2*eps_bx, 2*omeg_b, 2*eps_bx*eps_bz - 2*eps_by*omeg_b, (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*cos(q11) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q11), (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*cos(q11) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q11)]])

# geometric_jacob_pos = np.array([[1, 0, 0, -2*eps_by*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_by*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*eps_bz*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), -2*eps_bx*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), -(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) + (-l3*cos(q12) - l4*cos(q12 + q13))*(2*eps_bx**2 + 2*omeg_b**2 - 1), -l4*(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bx*eps_bz + 2*eps_by*omeg_b)*sin(q12 + q13)*cos(q11) - l4*(2*eps_bx**2 + 2*omeg_b**2 - 1)*cos(q12 + q13)], 
# 	[0, 1, 0, -2*eps_bx*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_by*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_by*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), (2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) - (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11), -l4*(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*cos(q12 + q13) + l4*(-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q12 + q13)*cos(q11) - l4*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11)*sin(q12 + q13)], 
# 	[0, 0, 1, -2*eps_bx*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) - 2*eps_by*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), 2*eps_bz*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_bz*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_bx*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*eps_by*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), (2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_bz**2 + 2*omeg_b**2 - 1)*cos(q11), -l4*(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*cos(q12 + q13) - l4*(2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q12 + q13)*cos(q11)]])

geometric_jacob_pos = np.array([[1, 0, 0, -2*eps_by*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_by*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*eps_bz*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), -2*eps_bx*((2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), -(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) + (-l3*cos(q12) - l4*cos(q12 + q13))*(2*eps_bx**2 + 2*omeg_b**2 - 1), -l4*(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bx*eps_bz + 2*eps_by*omeg_b)*sin(q12 + q13)*cos(q11) - l4*(2*eps_bx**2 + 2*omeg_b**2 - 1)*cos(q12 + q13)], [0, 1, 0, -2*eps_bx*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) - 2*eps_bz*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_by*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_bx*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) - 2*eps_bz*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))), 2*eps_by*(-(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), (2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) - (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11), -l4*(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*cos(q12 + q13) + l4*(-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q12 + q13)*cos(q11) - l4*(2*eps_by**2 + 2*omeg_b**2 - 1)*sin(q11)*sin(q12 + q13)], [0, 0, 1, -2*eps_bx*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) - 2*eps_by*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), 2*eps_bz*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*omeg_b*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), -2*eps_bz*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))) + 2*omeg_b*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))), -2*eps_bx*(-(2*eps_bx*eps_by - 2*eps_bz*omeg_b)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - (2*eps_bx*eps_bz + 2*eps_by*omeg_b)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) - (2*eps_bx**2 + 2*omeg_b**2 - 1)*(l1 - l3*sin(q12) - l4*sin(q12 + q13))) + 2*eps_by*((2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) + (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13))), (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), (2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(-l3*cos(q12) - l4*cos(q12 + q13)) - (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_bz**2 + 2*omeg_b**2 - 1)*cos(q11), -l4*(2*eps_bx*eps_bz - 2*eps_by*omeg_b)*cos(q12 + q13) - l4*(2*eps_bx*omeg_b + 2*eps_by*eps_bz)*sin(q11)*sin(q12 + q13) + l4*(2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q12 + q13)*cos(q11)]])
# w_cont_geom = -2*eps_by*( (2*eps_bx*eps_bz - 2*eps_by*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) 
# 	+ (2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) 
# 	+ (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13))) 
#     - 2*eps_bz*(-(2*eps_bx*eps_by + 2*eps_bz*omeg_b)*(l1 - l3*sin(q12) - l4*sin(q12 + q13)) 
#     - (-2*eps_bx*omeg_b + 2*eps_by*eps_bz)*(l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)) 
#     - (2*eps_by**2 + 2*omeg_b**2 - 1)*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)))

# w_cont_anal = -2*eps_by*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) 
#     - 2*eps_bz*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) 
#     - 4*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13))

# ============================================================================
# Jacobiano Analítico para la posición de la pata 1
# ============================================================================
analytic_jacob_pos = np.array([[1, 0, 0, -2*eps_by*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*eps_bz*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 4*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), -4*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), 2*eps_bx*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bx*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), 2*(eps_bx*eps_by - eps_bz*omeg_b)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + 2*(eps_bx*eps_bz + eps_by*omeg_b)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), 2*(-eps_bx*eps_by + eps_bz*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + 2*(eps_bx*eps_bz + eps_by*omeg_b)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) + (l3*cos(q12) + l4*cos(q12 + q13))*(-2*eps_bx**2 - 2*omeg_b**2 + 1), l4*(2*(-eps_bx*eps_by + eps_bz*omeg_b)*sin(q11)*sin(q12 + q13) + 2*(eps_bx*eps_bz + eps_by*omeg_b)*sin(q12 + q13)*cos(q11) + (-2*eps_bx**2 - 2*omeg_b**2 + 1)*cos(q12 + q13))],
 [0, 1, 0, 2*eps_bx*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*eps_bz*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 4*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), -2*eps_by*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 4*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 2*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_by*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) - 2*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), 2*(-eps_bx*omeg_b + eps_by*eps_bz)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + (2*eps_by**2 + 2*omeg_b**2 - 1)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*(eps_bx*eps_by + eps_bz*omeg_b)*(l3*cos(q12) + l4*cos(q12 + q13)) + 2*(-eps_bx*omeg_b + eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*cos(q11) + (l3*sin(q12) + l4*sin(q12 + q13))*(-2*eps_by**2 - 2*omeg_b**2 + 1)*sin(q11), l4*(-2*(eps_bx*eps_by + eps_bz*omeg_b)*cos(q12 + q13) + 2*(-eps_bx*omeg_b + eps_by*eps_bz)*sin(q12 + q13)*cos(q11) + (-2*eps_by**2 - 2*omeg_b**2 + 1)*sin(q11)*sin(q12 + q13))], 
 [0, 0, 1, 2*eps_bx*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + 2*eps_by*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) - 4*omeg_b*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), -2*eps_bz*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*omeg_b*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), 2*eps_bz*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) + 2*omeg_b*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)), -2*eps_bx*(-l1 + l3*sin(q12) + l4*sin(q12 + q13)) + 2*eps_by*(d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)) - 4*eps_bz*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)), 2*(eps_bx*omeg_b + eps_by*eps_bz)*(-l2*sin(q11) + l3*cos(q11)*cos(q12) + l4*cos(q11)*cos(q12 + q13)) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*(l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)), 2*(-eps_bx*eps_bz + eps_by*omeg_b)*(l3*cos(q12) + l4*cos(q12 + q13)) - 2*(eps_bx*omeg_b + eps_by*eps_bz)*(l3*sin(q12) + l4*sin(q12 + q13))*sin(q11) + (l3*sin(q12) + l4*sin(q12 + q13))*(2*eps_bz**2 + 2*omeg_b**2 - 1)*cos(q11), l4*(2*(-eps_bx*eps_bz + eps_by*omeg_b)*cos(q12 + q13) - 2*(eps_bx*omeg_b + eps_by*eps_bz)*sin(q11)*sin(q12 + q13) + (2*eps_bz**2 + 2*omeg_b**2 - 1)*sin(q12 + q13)*cos(q11))]])
#delta = 0.00001
# JT = 1/delta*np.array([
# 	fk_I1(pb_x+delta,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y+delta,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z+delta,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b+delta,eps_bx,eps_by,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx+delta,eps_by,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by+delta,eps_bz,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz+delta,q11,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11+delta,q12,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12+delta,q13) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13),
# 	fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13+delta) - fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13)])

# JA = JT.transpose()


q_init = np.array([[10],[10],[10],[1],[0],[0],[0],[0.3],[0.2],[0.6]])	
resul1 = analytic_jacob_pos.dot(q_init)
resul2 = geometric_jacob_pos.dot(q_init)
# resul1 = JA.dot(q_init)
resutlado = fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13)
eps_by = 0 
eps_bz = 0
eps_bx = 0
omeg_b = 1
resutlado2 = fk_I1(pb_x,pb_y,pb_z,omeg_b,eps_bx,eps_by,eps_bz,q11,q12,q13)
print("Cinemática directa")
print(resutlado)
print(resutlado2)
# print("Jacobiano Analítico")
# print(np.round(analytic_jacob_pos,2))
# print("Jacobiano Geométrico")
# print(np.round(geometric_jacob_pos,2))
# print("POS _ geom ")
# print(np.round(resul2))
# print("POS _ analitico ")
# print(np.round(resul1))
