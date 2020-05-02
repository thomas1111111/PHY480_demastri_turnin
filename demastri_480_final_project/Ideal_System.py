"""

This is the script to store the System Object in

The System object is made to set our boundaries, store our particles in, and handle the motion of the gas.

"""
import scipy.constants as scp
import numpy as np
import time
from IPython.display import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Ideal_Particle as ideal_ptr




class System:
    
    def __init__(self,dt = 0.0001, total_time = 100, tau = 300*scp.Boltzmann, m = scp.proton_mass, N=100, size = np.array([1,1,1])):
        #This initiates our system
        #The temperature will be in natural units (Energy)
        #Our container will be a rectangular prism
        #The default values set above are supposed to model a hydrogen gas at room temperature.
        
        self.particles = [None]*N #list of particles
        self.timestep = dt        #timestep for calculations
        self.temp = tau           #temperature of the system
        self.mass = m             #mass of our particles
        self.size = size          #size of our box
        self.runtime = total_time #sets a stopping time for our simulation
        self.p_col = 0            #momentum impulsed on the boundaries 
        
            
    def init_particles(self):
        
        #reset our particle locations
        self.particles = [None]*len(self.particles)
        
        #Scattering N particles inside our box:
        for i in range(len(self.particles)):
            
            #setting particle position
            par = ideal_ptr.Particle(mass = self.mass)
            par.pos = np.array([np.random.uniform(0, self.size[0]), np.random.uniform(0, self.size[1]), np.random.uniform(0, self.size[2])])
            
            #setting particle velocity
            #we set the direction by making a random unit vector 
            
            phi = np.random.uniform(0,2*np.pi)
            theta = np.random.uniform(0,np.pi)
            
            #set our velocity from the energy (temperature)
            vel_magnitude = np.sqrt((2*self.temp)/self.mass)
            par.vel = vel_magnitude*np.array([np.cos(theta)*np.sin(phi),np.sin(theta)*np.sin(phi),np.cos(phi)])
            
            self.particles[i] = par
            
            
    def in_box(self, position):
        
        #this function is to check if a position is outside of our box or not.
        
        size = self.size
        
        #If the particle has drifted out of the box, tell us where it is
        #otherwise, tell us it's in the box
        if (position[0] <= 0):
            return "below-x"
        elif (position[1] <= 0):
            return "below-y"
        elif (position[2] <= 0):
            return "below-z"
        elif (position[0] >= size[0]):
            return "above-x"
        elif (position[1] >= size[1]):
            return "above-y"
        elif (position[2] >= size[2]):
            return "above-z"
        else:
            return "in"    
        
    def update(self):
        
        '''
        
        This function uses the current velocity of each particle to update the positions of our particles during one timestep
        
        It also handles scattering the particle off the boundaries, and keeps track of how much momentum is impulsed on the wall.
        
        This information will be important during our calculation of pressure.
        
        '''
        
        pcol = 0 #momentum change from collisions
        
        for par in self.particles:
            
            new_pos = par.pos + self.timestep*par.vel #calculate new position based on trajectory
            check = self.in_box(new_pos)                   #check if it's in the box
            
            par.pos = new_pos
            if (check == "in"): #If it's in the box, then we don't have to do anything
                continue
            else:               #If it's outside the box:
               
                # Turn it around, and keep track of the momentum change
                if (check == "below-x"):    
                    par.pos[0] = 0
                    pcol += 2*par.mass*np.abs(par.vel[0])
                    par.vel[0] = -1 * par.vel[0]
                    continue
                if (check == "below-y"):
                    par.pos[1] = 0
                    pcol += 2*par.mass*np.abs(par.vel[1])
                    par.vel[1] = -1 * par.vel[1]
                    continue
                if (check == "below-z"):
                    par.pos[2] = 0
                    pcol += 2*par.mass*np.abs(par.vel[2])
                    par.vel[2] = -1 * par.vel[2]
                    continue
                if (check == "above-x"):    
                    par.pos[0] = self.size[0]
                    pcol += 2*par.mass*np.abs(par.vel[0])
                    par.vel[0] = -1 * par.vel[0]
                    continue
                if (check == "above-y"):
                    par.pos[1] = self.size[1]
                    pcol += 2*par.mass*np.abs(par.vel[1])
                    par.vel[1] = -1 * par.vel[1]
                    continue
                if (check == "above-z"):
                    par.pos[2] = self.size[2]
                    pcol += 2*par.mass*np.abs(par.vel[2])
                    par.vel[2] = -1 * par.vel[2]
                    continue
        return pcol
       
    def getparticlepos(self):
        
        #Iterate over the particle list,
        #append the particle positions to a list,
        #and return that as an array
        #This is here for if we want to animate the motion later.
        
        particlepos = []
        for i in range(len(self.particles)):
            particlepos.append(self.particles[i].pos)
        return np.array(particlepos)
    
    
    def run(self,plot=False):
        '''
        
        This initiates the particles and updates the system until we hit the maximum allowed time
        
        '''
        
        #If multiple runs of the same system are being done,
        #we want to reset the positions and the information about collisions.
        self.init_particles()
        self.p_col = 0
        
        if plot == False:
            
            t = 0
            
            while t<self.runtime:
                
                momentum = self.update()
                self.p_col += momentum
                t += self.timestep
                
            return
                
        else:
            
            t = 0
            positions = []
            
            while t<self.runtime:
                
                momentum = self.update()
                self.p_col += momentum
                positions.append(self.getparticlepos())
                t += self.timestep 
                
            return positions
        
                              
    def plot(self, positions, dt=0.1):
        
        '''
        
        This function animates the motion of the particles.
        
        It is only meant as a demonstration for small N, to check that things are working as we expect.
        
        It takes in a list of calculated positions, and displays where they were at
        a specified regular time interval.
        
        
        '''
        
        skip = dt/self.timestep
        
        for i in np.arange(0,len(positions),skip):
            i = int(i)
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot(111, projection='3d')
                
            pos = positions[i]
            
            for j in range(len(pos)):
                ax.scatter(pos[j][0],pos[j][1],pos[j][2],color="blue")
            ax.set_xlim([-1,self.size[0]+1])
            ax.set_ylim([-1,self.size[1]+1])
            ax.set_zlim([-1,self.size[2]+1])
            plt.grid() 
            clear_output(wait=True)
            display(fig)
            fig.clear()
        plt.close()

        
    def calc_pressure(self):
        #Grabbed this equation from Kittel and Kroemer's Thermal Physics, page 391.
        
        surface_area = 2*(self.size[0]*self.size[1] + self.size[1]*self.size[2] + self.size[0]*self.size[2])
        
        return (self.p_col/(self.runtime*surface_area))
    
    def pressure_error(self):
        #calculates the relative error of the pressure
        ideal_pressure = (len(self.particles)*self.temp)/(self.size[0]*self.size[1]*self.size[2])
        calculated_p = self.calc_pressure()
        return (np.abs((calculated_p-ideal_pressure)/ideal_pressure))
        
        