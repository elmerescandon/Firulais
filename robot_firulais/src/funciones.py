# -*- coding: utf-8 -*-
import numpy as np
from copy import copy

# ==================================
#  Funciones Numéricas Primarias
# ==================================


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

def crossproduct(a, b):
    '''
    Funcion numérica que retorna producto
    cruz de los vectores a y b (ambos arrays)
    '''
    x = np.array([[a[1]*b[2] - a[2]*b[1]],
                  [a[2]*b[0] - a[0]*b[2]],
                  [a[0]*b[1] - a[1]*b[0]]])
    return x


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

def fk_pata_pos(q,pata):
    if pata == 1:
    	T = fk_pata1_pos(q)
    elif pata == 2: 
    	T = fk_pata2_pos(q)
    elif pata == 3: 
    	T = fk_pata3_pos(q)
    elif pata == 4: 
    	T = fk_pata4_pos(q)
    return T


# ==============================================
# Cinemática directa de cada pata de cada transformación

def fk_pata1(q):
	# Pata 1
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q11 = q[0];q12 = q[1];q13 = q[2];
	T1_01 = sdh(l1, q11, 0, np.pi/2)
	T1_12 = sdh(l2, q12, -l3, 0)
	T1_23 = sdh(0, q13, -l4, 0)
	T1_03 = T1_01.dot(T1_12).dot(T1_23)
	T1_B0 = sTroty(-np.pi/2).dot(sTrotx(np.pi)).dot(sTrasl(0,-d1,0))
	T1_B1 = T1_B0.dot(T1_01)
	T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
	T1_B3 = T1_B0.dot(T1_03)
	T = [T1_B0,T1_B1,T1_B2,T1_B3]
	return T


def fk_pata2(q):
	# Pata 2
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q21 = q[0];q22 = q[1];q23 = q[2];
	T2_01 = sdh(l1, np.pi+q21, 0, np.pi/2)
	T2_12 = sdh(l2, q22, l3, 0)
	T2_23 = sdh(0, q23, l4, 0)
	T2_03 = (T2_01.dot(T2_12)).dot(T2_23)
	T2_B0 = sTroty(-np.pi/2).dot(sTrasl(0,d1,0))
	T2_B1 = T2_B0.dot(T2_01)
	T2_B2 = T2_B0.dot(T2_01).dot(T2_12)
	T2_B3 = T2_B0.dot(T2_03)
	T = [T2_B1,T2_B2,T2_B3]
	return T

def fk_pata3(q):
	# Pata 3
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075
	#l1 = 125;l2 = 25;l3 = 103;l4 = 104;d1 = 75
	q31 = q[0];q32 = q[1];q33 = q[2];
	T3_01 = sdh(l1, np.pi+q31, 0, np.pi/2)
	T3_12 = sdh(l2, q32, l3, 0)
	T3_23 = sdh(0, q33, l4, 0)
	T3_03 = T3_01.dot(T3_12).dot(T3_23)
	T3_B0 = sTroty(-np.pi/2).dot(sTrotx(np.pi)).dot(sTrasl(0,d1,0))
	T3_B1 = T3_B0.dot(T3_01)
	T3_B2 = T3_B0.dot(T3_01).dot(T3_12)
	T3_B3 = T3_B0.dot(T3_03)
	T = [T3_B1,T3_B2,T3_B3]	
	return T 


def fk_pata4(q):
	# Pata 4
	l1 = 0.125;l2=0.025;l3=0.105;l4=0.104;d1=0.075	
	q41 = q[0];q42 = q[1];q43 = q[2];
	T4_01 = sdh(l1, q41, 0, np.pi/2)
	T4_12 = sdh(l2, q42, -l3, 0)
	T4_23 = sdh(0, q43, -l4, 0)
	T4_03 = T4_01.dot(T4_12).dot(T4_23)
	T4_B0 = sTroty(-np.pi/2).dot(sTrasl(0,-d1,0))
	T4_B1 = T4_B0.dot(T4_01)
	T4_B2 = T4_B0.dot(T4_01).dot(T4_12)
	T4_B3 = T4_B0.dot(T4_03)
	T = [T4_B1,T4_B2,T4_B3]
	return T



# ===============================================
#  Funciones Cinemática Diferencial
# ===============================================

def jacob_g(q,pata):
    '''
    Función que realiza el Jacobiano geométrico dado los parámetros
    Denavit-Hatenberg. Retorna la matriz jacobiana considerando 
    la velocidad lineal y angular del sistema
    '''
    if pata == 1:
    	T_ref = fk_pata1(q)
    elif pata == 2: 
    	T_ref = fk_pata2(q)
    elif pata == 3: 
    	T_ref = fk_pata3(q)
    elif pata == 4: 
    	T_ref = fk_pata4(q)
    d1=0.075	
    p0 = T_ref[0][0:3,3]
    z0 = T_ref[0][0:3,2]
    T_ref = T_ref[1:4]
    J = [];
    for n in range(len(T_ref)):
        if (n == 0):
            z = z0
            p = p0
        else:
            z = T_ref[n-1][0:3, 2]
            p = T_ref[n-1][0:3, 3]

        # Todas son articulaciones de revolución
        Jv = crossproduct(z, T_ref[len(T_ref)-1][0:3, 3] - p)
        Jw = np.array([[z[0]],[z[1]],[z[2]]])
        J1 = np.vstack((Jv, Jw))
        if n == 0: 
        	J = J1
        else:
        	J = np.hstack((J, J1))
    return J


def jacob_a_pos(q,pata,delta=0.0001):
    '''
    Función que realiza el Jacobiano analítico de la posición
    de un pata con respecto a la base del robot
    '''
    # Crear una matriz 3x3
    J = np.zeros((3,3))
    # Transformacion homogenea inicial (usando q)
    T = fk_pata_pos(q,pata)
    # Iteracion para la derivada de cada columna
    for i in xrange(3):
        # Copiar la configuracion articular inicial
        dq = copy(q)
        # Incrementar la articulacion i-esima usando un delta
        dq[i]=dq[i]+delta
        # Transformacion homogenea luego del incremento (q+delta)
        dT=fk_pata_pos(dq,pata)
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J[:,[i]]=(dT[:3,[3]]-T[:3,[3]])/delta
    return J
