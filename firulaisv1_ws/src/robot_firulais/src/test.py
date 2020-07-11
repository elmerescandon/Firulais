import numpy as np 
from funciones import *


q = [0.3,0,0]
x1 = fk_pata1_pos(q);
x2 = fk_pata2_pos(q);
x3 = fk_pata3_pos(q);
x4 = fk_pata4_pos(q);
print(np.round(x1,2))
print(np.round(x2,2))
print(np.round(x3,2))
print(np.round(x4,2))