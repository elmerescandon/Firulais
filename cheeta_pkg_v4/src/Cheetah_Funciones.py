import numpy as np
from copy import copy
from funciones import *
  

def fkine_cheetah(q):
    """
    Calcular la cinematica directa del robot UR5 dados sus valores articulares. 
    q es un vector numpy de la forma [q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12]
    """
    # Longitudes (en metros)
    l0 = 0.054238 #
    l1 = 0.11694 #
    l2 = 0.09243 #
    l3 = 0.1075 #
    d0 = 0.012455 #
    d1 = 0.002717 #
    d2 = 0.025355 #
    d3 = 0.035 #

    # Matrices DH (completar), emplear la funcion dh con los parametros DH para cada articulacion Pata1
    T10 = Troty(-np.pi/2)*Trotx(np.pi)*Trasl(d0,-l0,0)
    T11 = dh(l1, q[0], l1, np.pi/2) 
    T12 = dh(l2, q[1], l2, 0)
    T13 = dh(l3, q[2], l3, 0)
        # Efector final con respecto a la base
    T1 = T10.dot(T11).dot(T12).dot(T13)

    # Matrices DH (completar), emplear la funcion dh con los parametros DH para cada articulacion Pata2
    T20 = Troty(-np.pi/2)*Trasl(d0,l0,0)
    T21 = dh(l1, np.pi+q[3], d1, np.pi/2)
    T22 = dh(l2, q[4], d2, 0)
    T23 = dh(l3, q[5], d3, 0)
        # Efector final con respecto a la base
    T2 = T20.dot(T21).dot(T22).dot(T23)

    # Matrices DH (completar), emplear la funcion dh con los parametros DH para cada articulacion Pata3
    T30 = Troty(-np.pi/2)*Trotx(np.pi)*Trasl(d0,l0,0)
    T31 = dh(l1, np.pi+q[6], d1, np.pi/2)
    T32 = dh(l2, q[7], d2, 0)
    T33 = dh(l3, q[8], d3, 0)
        # Efector final con respecto a la base
    T3 = T30.dot(T31).dot(T32).dot(T33)

    # Matrices DH (completar), emplear la funcion dh con los parametros DH para cada articulacion Pata4
    T40 = Troty(-np.pi/2)*Trasl(d0,-l0,0)
    T41 = dh(l1, q[9], d1, np.pi/2)
    T42 = dh(l2, q[10], d2, 0)
    T43 = dh(l3, q[11], d3, 0)
        # Efector final con respecto a la base
    T4 = T40.dot(T41).dot(T42).dot(T43)
    return T1,T2,T3,T4

def jacobian_cheetah(q, delta=0.0001):
    """
    Jacobiano analitico para la posicion. Retorna una matriz de 3x12 y toma como
    entrada el vector de configuracion articular q=[q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12]
    """
    # Crear una matriz 3x12
    J1i = np.zeros((3,3))
    J2i = np.zeros((3,3))
    J3i = np.zeros((3,3))
    J4i = np.zeros((3,3))

    # Transformacion homogenea inicial (usando q)
    T103_q, T203_q, T303_q, T403_q = fkine_cheetah(q)
    P103_q = T103_q[0:3,3]
    P203_q = T203_q[0:3,3]
    P303_q = T303_q[0:3,3]
    P403_q = T403_q[0:3,3]
    
    # Iteracion para la derivada de cada columna
    for i in xrange(3):
        # Copiar la configuracion articular inicial
        dq = copy(q)
        # Incrementar la articulacion i-esima usando un delta    

        dq[i] = dq[i] + delta
        dq[i+3] = dq[i+3] + delta
        dq[i+6] = dq[i+6] + delta
        dq[i+9] = dq[i+9] + delta
        # Transformacion homogenea luego del incremento (q+delta)
        T103_dq, T203_dq, T303_dq, T403_dq = fkine_cheetah(dq)
        P103_dq = T103_dq[0:3,3]
        P203_dq = T203_dq[0:3,3]
        P303_dq = T303_dq[0:3,3]
        P403_dq = T403_dq[0:3,3]
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J1i[:,i] = (P103_dq - P103_q)/delta
        J2i[:,i] = (P203_dq - P203_q)/delta
        J3i[:,i] = (P303_dq - P303_q)/delta
        J4i[:,i] = (P403_dq - P403_q)/delta
        
    return J1i, J2i, J3i, J4i


def ikine_cheetah(xdes, q0):
    """
    Calcular la cinematica inversa de UR5 numericamente a partir de la configuracion articular inicial de q0. 
    """
    #Parametros
    epsilon  = 0.0001
    max_iter = 1000
    delta    = 0.00001
    # Copiar la configuracion articular inicial
    q  = copy(q0)
    e = np.zeros((12))
    # Main loop
    for i in range(max_iter):
        #Calcular jacobiano analitico
        J1, J2, J3, J4 = jacobian_cheetah(q, delta)
        #Hallar la posicion del efector final
        T103_q, T203_q, T303_q, T403_q = fkine_cheetah(q)
        f1 = T103_q[0:3,3]
        f2 = T203_q[0:3,3]
        f3 = T303_q[0:3,3]
        f4 = T403_q[0:3,3]
        
        #Error de posicion
        e1 = xdes[0:3]-f1
        e2 = xdes[3:6]-f2
        e3 = xdes[6:9]-f3
        e4 = xdes[9:12]-f4

        q[0:3] = q[0:3] + np.dot(np.linalg.pinv(J1), e1)
        q[3:6] = q[3:6] + np.dot(np.linalg.pinv(J2), e2)
        q[6:9] = q[6:9] + np.dot(np.linalg.pinv(J3), e3)
        q[9:12] = q[9:12] + np.dot(np.linalg.pinv(J4), e4)

        e[0:3] = e1
        e[3:6] = e2
        e[6:9] = e3
        e[9:12] = e4

        #Condicion de termino
        if (np.linalg.norm(e) < epsilon):
            break

    return q