# -*- coding: utf-8 -*-
import numpy as np
# Funciones útiles 

cos=np.cos; sin=np.sin; pi=np.pi
import rbdl


class Robot(object):
    def __init__(self, q0, dq0, ndof, dt):
        self.q = q0    # numpy array (ndof x 1)
        self.dq = dq0  # numpy array (ndof x 1)
        self.M = np.zeros([ndof, ndof])
        self.b = np.zeros(ndof)
        self.dt = dt
        self.robot = rbdl.loadModel('../urdf/cheeta_pkg_v4.urdf')

    def send_command(self, tau):
        rbdl.CompositeRigidBodyAlgorithm(self.robot, self.q, self.M)
        rbdl.NonlinearEffects(self.robot, self.q, self.dq, self.b)
        ddq = np.linalg.inv(self.M).dot(tau-self.b)
        self.q = self.q + self.dt*self.dq
        self.dq = self.dq + self.dt*ddq

    def read_joint_positions(self):
        return self.q

    def read_joint_velocities(self):
        return self.dq



def dh(d, theta, a, alpha):
    """
    Calcular la matriz de transformacion homogenea asociada con los parametros
    de Denavit-Hartenberg.
    Los valores d, theta, a, alpha son escalares.
    """
    # Escriba aqui la matriz de transformacion homogenea en funcion de los valores de d, theta, a, alpha
    T = np.array([[cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta)],
                    [sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
                    [0,        sin(alpha),     cos(alpha),      d],
                    [0,         0,      0,      1]])
    return T
    

def Trasl(x, y, z):
    """ Matriz de transformada homogenea de traslacion
    """
    T = np.array([[1,0,0,x],
                    [0,1,0,y],
                    [0,0,1,z],
                    [0,0,0,1]])
    return T

def Trotx(ang):
    """ Matriz de transformada homogenea alrededor de X
    """
    Tx = np.array([[1, 0,0,0],
                    [0, cos(ang),-sin(ang),0],
                    [0, sin(ang), cos(ang),0],
                    [0, 0, 0, 1]])
    return Tx

def Troty(ang):
    """ Matriz de transformada homogenea alrededor de Y
    """
    Ty = np.array([[cos(ang),0,sin(ang),0],
                    [0,1,0,0],
                    [-sin(ang),0,cos(ang),0],
                    [0,0,0,1]])
    return Ty

def Trotz(ang):
    """ Matriz de transformada homogenea alrededor de Z
    """
    Tz = np.array([[cos(ang),-sin(ang),0,0],
                    [sin(ang), cos(ang),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])
    return Tz