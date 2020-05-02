'''

This is our particle class. It is basically a data structure.

I plan for it to have the following attributes:
    pos - position, 3d numpy array
    vel - velocity, 3d numpy array
    mass - mass, scalar (float)
    
'''

import scipy.constants as scp
import numpy as np

class Particle:
    def __init__(self, mass=scp.proton_mass,pos=np.array([0,0,0]),vel=np.array([0,0,0])):
        self.mass = mass
        self.pos = pos
        self.vel = vel
