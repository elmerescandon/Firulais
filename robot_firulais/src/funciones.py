# -*- coding: utf-8 -*-
import numpy as np


# ===============================
#  Funciones Numéricas Primarias
# ===============================


def sdh(d, th, a, alpha):
    """ Matriz de transformación homogénea de Denavit-Hartenberg
    	Entradas: parámetros
    	Salida: Matriz tranformación
    """
    cth = np.cos(th); sth = np.sin(th)
    ca = np.cos(alpha); sa = np.sin(alpha)
    Tdh = np.array([ [cth, -ca*sth,  sa*sth, a*cth],
                     [sth,  ca*cth, -sa*cth, a*sth],
                     [0,        sa,     ca,      d],
                     [0,         0,      0,      1]])
    return Tdh


def rot2quaternion(R):
    """
    Funcion que retorna el vector de quaterion unitario
    a partir de una matriz de rotacion.
    No considera la forma adicional de operar cuando el angulo es 180
    Lo de vuelve de la forma:
    q = (w,ex,ey,ez)
    """
    omega = ((1+R[0, 0]+R[1, 1]+R[2, 2])**0.5)*0.5
    ex = (1/(4*omega))*(R[2, 1]-R[1, 2])
    ey = (1/(4*omega))*(R[0, 2]-R[2, 0])
    ez = (1/(4*omega))*(R[1, 0]-R[0, 1])
    q = np.array([[omega],
                  [ex],
                  [ey],
                  [ez]])
    return q

def sTrasl(x, y, z):
    """ Matriz de transformada homogenea de traslacion
    """
    T = np.array([[1,0,0,x],
                    [0,1,0,y],
                    [0,0,1,z],
                    [0,0,0,1]])
    return T

def sTrotx(ang):
    """ Matriz de transformada homogenea alrededor de X
    """

    Tx = np.array([[1, 0,0,0],
                    [0, np.cos(ang),-np.sin(ang),0],
                    [0, np.sin(ang), np.cos(ang),0],
                    [0, 0, 0, 1]])
    return Tx

def sTroty(ang):
    """ Matriz de transformada homogenea alrededor de Y
    """
    Ty = np.array([[np.cos(ang),0,np.sin(ang),0],
                    [0,1,0,0],
                    [-np.sin(ang),0,np.cos(ang),0],
                    [0,0,0,1]])
    return Ty

def sTrotz(ang):
    """ Matriz de transformada homogenea alrededor de Z
    """
    Tz = np.Matrix([[np.cos(ang),-np.sin(ang),0,0],
                    [np.sin(ang), np.cos(ang),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])
    return Tz

# =================================
#  Funciones Cinemática Directa
# =================================

# Cinemática directa de cada pata considerando
# únicamente la posición 

def fk_pata1_pos(q):
	# Pata 1
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q11 = q[0];q12 = q[1];q13 = q[2];
	T1_01 = sdh(l1, q11, 0, np.pi/2)
	T1_12 = sdh(l2, q12, -l3, 0)
	T1_23 = sdh(0, q13, -l4, 0)
	T1_03 = T1_01.dot(T1_12).dot(T1_23)
	T1_B0 = sTroty(-np.pi/2).dot(sTrotx(np.pi)).dot(sTrasl(0,-d1,0))
	T1_B3 = T1_B0.dot(T1_03)
	return T1_B3


def fk_pata2_pos(q):
	# Pata 2
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q21 = q[0];q22 = q[1];q23 = q[2];
	T2_01 = sdh(l1, np.pi+q21, 0, np.pi/2)
	T2_12 = sdh(l2, q22, l3, 0)
	T2_23 = sdh(0, q23, l4, 0)
	T2_03 = (T2_01.dot(T2_12)).dot(T2_23)
	T2_B0 = sTroty(-np.pi/2).dot(sTrasl(0,d1,0))
	T2_B3 = T2_B0.dot(T2_03)
	return T2_B3

def fk_pata3_pos(q):
	# Pata 3
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q31 = q[0];q32 = q[1];q33 = q[2];
	T3_01 = sdh(l1, np.pi+q31, 0, np.pi/2)
	T3_12 = sdh(l2, q32, l3, 0)
	T3_23 = sdh(0, q33, l4, 0)
	T3_03 = T3_01.dot(T3_12).dot(T3_23)
	T3_B0 = sTroty(-np.pi/2).dot(sTrotx(np.pi)).dot(sTrasl(0,d1,0))
	T3_B3 = T3_B0.dot(T3_03)
	return T3_B3 


def fk_pata4_pos(q):
	# Pata 4
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075	
	q41 = q[0];q42 = q[1];q43 = q[2];
	T4_01 = sdh(l1, q41, 0, np.pi/2)
	T4_12 = sdh(l2, q42, -l3, 0)
	T4_23 = sdh(0, q43, -l4, 0)
	T4_03 = T4_01.dot(T4_12).dot(T4_23)
	T4_B0 = sTroty(-np.pi/2).dot(sTrasl(0,-d1,0))
	T4_B3 = T4_B0.dot(T4_03)
	return T4_B3

# ==============================================
# Cinemática directa de cada pata considerando
# el estaoo de posición y orientación (quaternion)

def fk_pata1(q):
	# Pata 1
	l1 = 0.125;l2=0.025;l3=0.103;l4=0.104;d1=0.075
	q11 = q[0];q12 = q[1];q13 = q[2];
	T1_01 = sdh(l1, q11, 0, np.pi/2)
	T1_12 = sdh(l2, q12, -l3, 0)
	T1_23 = sdh(0, q13, -l4, 0)
	T1_03 = T1_01.dot(T1_12).dot(T1_23)
	R = T1_03[0:3,0:3]
	quat_1 = rot2quaternion(R)
	x  = np.array([[T1_03[0,3]],
				   [T1_03[1,3]],
				   [T1_03[2,3]],
				   [quat_1[0,0]],
				   [quat_1[1,0]],
				   [quat_1[2,0]],
				   [quat_1[3,0]]])
	return x


def fk_pata2(q):
	# Pata 2
	l1 = 0.125;l2=0.025;l3=0.103;l4=0.104;d1=0.075
	q21 = q[0];q22 = q[1];q23 = q[2];
	T2_01 = sdh(l1, np.pi+q21, 0, np.pi/2)
	T2_12 = sdh(l2, q22, l3, 0)
	T2_23 = sdh(0, q23, l4, 0)
	T2_03 = T2_01.dot(T2_12).dot(T2_23)
	R = T2_03[0:3,0:3]
	quat_2 = rot2quaternion(R)
	x  = np.array([[T2_03[0,3]],
				   [T2_03[1,3]],
				   [T2_03[2,3]],
				   [quat_2[0,0]],
				   [quat_2[1,0]],
				   [quat_2[2,0]],
				   [quat_2[3,0]]])
	return x

def fk_pata3(q):
	# Pata 3
	l1 = 0.125;l2=0.025;l3=0.103;l4=0.104;d1=0.075
	q31 = q[0];q32 = q[1];q33 = q[2];
	T3_01 = sdh(l1, np.pi+q31, 0, np.pi/2)
	T3_12 = sdh(l2, q32, l3, 0)
	T3_23 = sdh(0, q33, l4, 0)
	T3_03 = T3_01.dot(T3_12).dot(T3_23)
	R = T3_03[0:3,0:3]
	quat_3 = rot2quaternion(R)
	x  = np.array([[T3_03[0,3]],
				   [T3_03[1,3]],
				   [T3_03[2,3]],
				   [quat_3[0,0]],
				   [quat_3[1,0]],
				   [quat_3[2,0]],
				   [quat_3[3,0]]])
	return x


def fk_pata4(q):
	l1 = 0.125;l2=0.025;l3=0.103;l4=0.104;d1=0.075
	q41 = q[0];q42 = q[1];q43 = q[2];
	T4_01 = sdh(l1, q41, 0, np.pi/2)
	T4_12 = sdh(l2, q42, -l3, 0)
	T4_23 = sdh(0, q43, -l4, 0)
	T4_03 = T4_01.dot(T4_12).dot(T4_23)
	R = T4_03[0:3,0:3]
	quat_4 = rot2quaternion(R)
	x  = np.array([[T4_03[0,3]],
				   [T4_03[1,3]],
				   [T4_03[2,3]],
				   [quat_4[0,0]],
				   [quat_4[1,0]],
				   [quat_4[2,0]],
				   [quat_4[3,0]]])
	return x



