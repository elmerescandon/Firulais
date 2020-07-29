# -*- coding: utf-8 -*-
import numpy as np
from copy import copy
from pyquaternion import Quaternion

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
    if omega == 0: 
        quat = Quaternion(matrix=R)
        q = np.array([q[0],q[1],q[2],q[3]])
        return q
    else :
        ex = (1/(4*omega))*(R[2, 1]-R[1, 2])
        ey = (1/(4*omega))*(R[0, 2]-R[2, 0])
        ez = (1/(4*omega))*(R[1, 0]-R[0, 1])
        return np.array([omega,ex,ey,ez])

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
    Tz = np.array([[np.cos(ang),-np.sin(ang),0,0],
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

def quaternionMult(q1, q2):

    qout = np.zeros(4)
    qout[0] = -q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3] + q1[0] * q2[0]
    qout[1] = q1[0] * q2[1] - q1[3] * q2[2] + q1[2] * q2[3] + q1[1] * q2[0]
    qout[2] = q1[3] * q2[1] + q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0]
    qout[3] = -q1[2] * q2[1] + q1[1] * q2[2] + q1[0] * q2[3] + q1[3] * q2[0]
    return qout

def error_quater(qact, qd):
    print(qd);print(qact)
    # Conversión a vector columna
    eps_act = np.array([[qact[1]],[qact[2]],[qact[3]]])
    eps_d = np.array([[qd[1]],[qd[2]],[qd[3]]])
    # Obtención de ew_act y ew_des 
    ew_act = qact[0]
    ew_d = qd[0]
    # Error de quaterniones
    ew_err = ew_d*ew_act + (eps_d.transpose()).dot(eps_act) 
    cross = crossproduct(eps_d,eps_act) 
    cross = np.array([cross[0,0],cross[1,0],cross[2,0]])
    eps_err= -ew_d*eps_act + ew_act*eps_d - cross
    quater_err = np.array([ew_err,[eps_err[0,0]],[eps_err[1,0]],[eps_err[2,0]]])
    return quater_err

def error_quaterv2(qact, qd):
    # Conversión a vector columna
    qact = np.array([qact[0],-qact[1],-qact[2],-qact[3]])
    quater_err =  quaternionMult(qd,qact)
    return np.array([[quater_err[0]],[quater_err[1]],[quater_err[2]],[quater_err[3]]])


# =================================
#  Funciones Cinemática Directa
# =================================

# Cinemática directa de cada pata considerando
# únicamente la posición 

def fk_pata1_pos(q, modo=''):
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(dx, dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, -0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] - np.deg2rad(90-76.111), -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] + np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
    if modo == 'pose':
        quater = rot2quaternion(T1_B3[0:3,0:3])
        position = T1_B3[0:3,3] 
        pose = np.hstack((position,quater))
        return pose
    else:
        return T1_B3


def fk_pata2_pos(q, modo=''):
    # Pata 2
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(-dx, dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, -0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] - np.deg2rad(90-76.111), -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] + np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
    
    if modo=='pose':
        quater = rot2quaternion(T1_B3[0:3,0:3])
        position = T1_B3[0:3,3] 
        pose = np.hstack((position,quater))
        return pose
    else:
        return T1_B3

def fk_pata3_pos(q, modo=''):
    # Pata 3
	dx = 0.116940
	dy = 0.055642
	dz = 0.012717
	T1_B0 = sTrasl(dx, -dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(-np.pi / 2))
	T1_01 = sdh(0, q[0]+np.pi/2, 0.01, np.pi / 2)
	T1_12 = sdh(0, q[1] + np.deg2rad(90-76.111) + np.pi, -0.105627, 0)
	T1_23 = sdh(0.0732435, q[2] - np.deg2rad(90-76.111) + np.pi, 0.1, 0)

	# T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
	T1_B1 = T1_B0.dot(T1_01)
	T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
	T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
	if modo == 'pose':
		quater = rot2quaternion(T1_B3[0:3,0:3])
		position = T1_B3[0:3,3] 
		pose = np.hstack((position,quater))
		return pose
	else:
		return T1_B3


def fk_pata4_pos(q, modo=''):
    # Pata 4
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(-dx, -dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(-np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, 0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] + np.deg2rad(90-76.111) + np.pi, -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] - np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
    if modo == 'pose':
        quater = rot2quaternion(T1_B3[0:3,0:3])
        position = T1_B3[0:3,3] 
        pose = np.hstack((position,quater))
        return pose
    else:
        return T1_B3

def fk_pata_pos(q,pata,modo=''):
    if pata == 1:
    	T = fk_pata1_pos(q,modo)
    elif pata == 2: 
    	T = fk_pata2_pos(q,modo)
    elif pata == 3: 
    	T = fk_pata3_pos(q,modo)
    elif pata == 4: 
    	T = fk_pata4_pos(q,modo)
    return T


# ==============================================
# Cinemática directa de cada pata de cada transformación

def fk_pata1(q):
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(dx, dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, -0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] - np.deg2rad(90-76.111), -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] + np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)


    T = [T1_B0,T1_01,T1_12,T1_23]
    return T


def fk_pata2(q):
	# Pata 2
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(-dx, dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, -0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] - np.deg2rad(90-76.111), -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] + np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)

    T = [T1_B0,T1_01,T1_12,T1_23]
    return T

def fk_pata3(q):
	# Pata 3
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(dx, -dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(-np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, 0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] + np.deg2rad(90-76.111) + np.pi, -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] - np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
    T = [T1_B0,T1_01,T1_12,T1_23]
    return T 


def fk_pata4(q):
	# Pata 4
    dx = 0.116940
    dy = 0.055642
    dz = 0.012717
    T1_B0 = sTrasl(-dx, -dy, dz).dot(sTroty(np.pi / 2)).dot(sTrotz(-np.pi / 2))
    T1_01 = sdh(0, q[0]+np.pi/2, 0.01, np.pi / 2)
    T1_12 = sdh(0, q[1] + np.deg2rad(90-76.111) + np.pi, -0.105627, 0)
    T1_23 = sdh(0.0732435, q[2] - np.deg2rad(90-76.111) + np.pi, 0.1, 0)

    # T1_B0 = sTrasl(bx,by,bz).dot(sTroty(np.pi/2)).dot(sTrotz(np.pi))
    T1_B1 = T1_B0.dot(T1_01)
    T1_B2 = T1_B0.dot(T1_01).dot(T1_12)
    T1_B3 = T1_B0.dot(T1_01).dot(T1_12).dot(T1_23)
    T = [T1_B0,T1_01,T1_12,T1_23]
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
    T = fk_pata_pos(q,  pata)
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

def jacobian_a_pose(q,pata,delta=0.0001):
    """
    Jacobiano analitico para la posicion y orientacion (usando un
    cuaternion). Retorna una matriz de 7x6 y toma como entrada el vector de
    configuracion articular q=[q1, q2, q3, q4, q5, q6]

    """
    # Crear una matriz 3x6
    JT = np.zeros((3,7))
    # Transformacion homogenea inicial (usando q)
    T=fk_pata_pos(q,pata)
    quater = rot2quaternion(T[0:3,0:3])
    position = T[0:3,3] # Posición en el espacio cartesiano
    # Pose tiene la forma = x,y,z,ew,ex,ey,ez
    pose = np.hstack((position,quater))
    # Iteracion para la derivada de cada columna
    for i in xrange(3):
        # Copiar la configuracion articular inicial
        dq = copy(q)
        # Incrementar la articulacion i-esima usando un delta
        dq[i]=dq[i]+delta#incremento de delta en la artic. i 
        # Transformacion homogenea luego del incremento (q+delta)
        dT=fk_pata_pos(dq,pata)
        dquater = rot2quaternion(dT[0:3,0:3])
        dposition = dT[0:3,3] # Posición en el espacio cartesiano
        # Pose tiene la forma = x,y,z,ew,ex,ey,ez
        dpose = np.hstack((dposition,dquater))
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        JT[[i],:]=(dpose-pose)/delta#se va ajustando artic. por artic.
        J = JT.T
    return J

def jacobian_a_posev2(q, pata, delta=0.0001):
    """
    Jacobiano analitico para la posicion y orientacion (usando un
    cuaternion). Retorna una matriz de 7x6 y toma como entrada el vector de
    configuracion articular q=[q1, q2, q3, q4, q5, q6]

    """
    #J = np.zeros((7,6))
    # Transformacion homogenea inicial (usando q)
    q1 = q[0];q2=q[1];q3=q[2];

    # Implementar este Jacobiano aqui
    JT = 1/delta*np.array([fk_pata_pos(np.array([q1+delta,q2,q3]),pata,'pose')-fk_pata_pos(q,pata,'pose'),
                  fk_pata_pos(np.array([q1,q2+delta,q3]),pata,'pose')-fk_pata_pos(q,pata,'pose'),
                  fk_pata_pos(np.array([q1,q2,q3+delta]),pata,'pose')-fk_pata_pos(q,pata,'pose')])
    J = JT.transpose()
    return J

# ===============================================
#  Funciones de control cinemático diferencial
# ===============================================

def control_fkdiff(x,xd,q1,dt,pata,k=1):
    """
     Función que resuelve la tarea del error para obtener
     la velocidad angular de cada articulación
     Entradas: 
        - x (posición actual - Espacio Operacional)
        - xd (posición deseada - Espacio Operacional)
        - q1 (posición actual - espacio articular)
        - k (Ganancia proporicual - Lambda) / Por defecto = 1   
    """ 
    ep = xd[0:3] - x[0:3] 
    Qe = error_quaterv2(x[3:7], xd[3:7])
    e0 = np.array([[Qe[0,0]-1],
           [Qe[1,0]],
           [Qe[2,0]],
           [Qe[3,0]]])
    e = np.array([[ep[0]],[ep[1]],[ep[2]],[e0[0,0]],[e0[1,0]],[e0[2,0]],[e0[3,0]]]) 

    e_dot = k*e
    J = jacobian_a_posev2(q1,pata)
    try: 
        J_mul = np.linalg.pinv(J)
    except np.linalg.LinAlgError:
        J_mul = (J.transpose()).dot(np.linalg.inv(J.dot(J.transpose()) + k*np.eye(7)))

    q_dot = J_mul.dot(e_dot)
    q_dot_vect = np.array([q_dot[0,0],q_dot[1,0],q_dot[2,0]])
    q1 =  q1 + dt*q_dot_vect
    return q1


# Actualización de variables - Estado Inicial
def update_initial_state(q0):
    x_pos = []
    T_pos = []
    pata_values = [[3,6],[9,12],[0,3],[6,9]] 
    for a in range(4):
        tmp = pata_values[a]
        q = q0[tmp[0]:tmp[1]]
        T = fk_pata_pos(q,a+1)
        quater = rot2quaternion(T[0:3,0:3])
        position = T[0:3,3] # Posición en el espacio cartesiano
        x = np.hstack((position,quater))
        x_pos.append(x)
        T_pos.append(T)
    return x_pos,T_pos