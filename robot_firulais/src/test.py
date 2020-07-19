import numpy as np 
from funciones import *


q = np.array([3,-0.5,0])

g = jacobian_a_pose(q,4,delta=0.0001)
print(g)
	