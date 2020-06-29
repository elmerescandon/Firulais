# -*- coding: utf-8 -*-

import numpy as np 

'''
Para pata 1 se necesita: 
	- Cinemática directa de la pata 1 con respecto al sis inercial 
	- Jacobiano geométrico de la pata 1 con respecto al sis inercial
'''

def fk1(q):
	l1 = 125;l2=25;l3=105;l4=104;d1=75
	q11 = q[0]
	q12 = q[1] 
	q13 = q[2]
	sin = np.sin
	cos = np.cos
	matriz_t = np.array([[sin(q12 + q13), cos(q12 + q13), 0, l1 - l3*sin(q12) - l4*sin(q12 + q13)], 
	    				 [-sin(q11)*cos(q12 + q13), sin(q11)*sin(q12 + q13), cos(q11), d1 + l2*cos(q11) + l3*sin(q11)*cos(q12) + l4*sin(q11)*cos(q12 + q13)], 
	    				 [cos(q11)*cos(q12 + q13), -sin(q12 + q13)*cos(q11), sin(q11), l2*sin(q11) - l3*cos(q11)*cos(q12) - l4*cos(q11)*cos(q12 + q13)], 
	    				 [0, 0, 0, 1]])
	return matriz_t[0:3,3]


# Configuración inicial
q = np.array[[],[],[],# Posición de b.f
				  [],[],[],[], # Orientación de b.f
				  [],[],[],[],[],[],[],[],[],[],[],[]] # Grados de articulaciones 

# Posición deseada (únicamente cartesiana R3)
p_d = np.array([[-10],[-10],[-100]])
p_act = fk_1 (q_init)
l = 1 # Ganancia de la tarea cinemática
min_error = 0.01

while True:
	# Obtener la cinemática directa y diferencial
	Jg_pos = Jgp1(q_init)
	p_act = fk_1(q)

	# Error 
	e = p_d - p_act 

	# Error 
	q_dot = l*(np.linalg.pinv(Jg_pos).dot(e))
	q = q + q_dot*dt

	if (np.linalg.norm(e) < min_error )
		break

