import numpy as np 
from funciones import *


ang = np.linspace(0,6*np.pi,300)
temp = np.mod(ang,2*np.pi)
radio = 0.10 
z = radio*np.sin(temp)
x = radio*np.cos(temp)
print(x)


	