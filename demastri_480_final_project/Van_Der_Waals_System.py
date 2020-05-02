"""

This is the script to store the Vad Der Walls System in in

This object is made to be an extension of the Ideal_System object. 

It implements scattering of the particles off one another and
adds inter-particle interactions via the Velocity-Verlet Algorithm.

It is thus much slower than the Ideal_System.

"""
import scipy.constants as scp
import numpy as np
import time
from IPython.display import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Ideal_Particle as ideal_ptr




class VDW_System:
    
    def __init__(self,dt = 0.0001, total_time = 100, tau = 300*scp.Boltzmann, m = scp.proton_mass, 
                 N=100, size = np.array([1,1,1]),ic = scp.proton_mass, b_0 = 0.005):
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
        self.interaction_const = ic #the interaction should be be proportional to 1/(r**2). I'm setting the default value to the mass of 
                                    #the particle
        self.par_radius = b_0     #the particles in a van der walls gas have hard boundaries, this is their radius
        
            
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
        
        
    def seperation_vectors(self):
        
        '''
        
        This finds the seperation distance between particles.
        
        My idea is that the seperation distance between each particle at the current time will 
        be calculated and then all particles will be moved forward.
        
        That way no particle is given "priority" over another particle.
        
        '''
        
        N = len(self.particles)
        spr_keys = []
        spr_dists = []
        
        for i in range(N-1):
            
            #To avoid calculating seperation distances twice,
            #the following block of code creates a list of integers
            #strictly larger than i and less than the total particle number
            
            check_list = range(i,N,1)
            
            #Now iterate over our particles, and calculate
            #the seperation distance (also keep track of which
            #distance we're calculating).
            
            for j in check_list:
                spr_keys.append(str(i)+str(j))
                spr_list = []
                for l in range(3):
                    spr_list.append(self.particles[i].pos[l]-self.particles[j].pos[l])
                spr_vector = np.array(spr_list)
                spr_dists.append(spr_vector)
            
        #Now let's compound these into a dictionary:
        spr_dict = dict(zip(spr_keys,spr_dists))
            
        return spr_dict
    
    def netForce(self):
        
        '''
        
        Applies the interaction between particles and scatters them if they get too close.
        
        '''
        
        spr_dict = self.seperation_vectors()
        
        #We're going to initiate a list of forces,
        #currently just N '0 vectors'
        spr_dict = self.seperation_vectors()
        N = len(self.particles)
        F = [np.array([0,0,0],dtype=np.float64)]*N
        
        for key in spr_dict.keys():
            
            i = int(key[0])
            j = int(key[1])
            
            #Calculating the unit vectors of our seperation vectors, and the seperation distance
            r = spr_dict[key]
            script_r = np.linalg.norm(r)
            
            #If they get too close, turn them around and scatter them off each other
            if script_r<=self.par_radius:
                self.particles[i].vel = -1 * self.particles[i].vel
                self.particles[j].vel = -1 * self.particles[j].vel
                continue
                    
            r_norm_ij = -r/script_r
            r_norm_ji = r/script_r
            
            #Calulating the interaction force from our knowns
            const = self.interaction_const/(script_r**2)
            
            force_ij = const*r_norm_ij
            force_ji = const*r_norm_ji
            
            F[j] = F[j]+force_ij
            F[i] = F[i]+force_ji
            
 
        return F
        
    def update(self):
        
        '''
        
        This function uses the current velocity of each particle to update the positions of our particles during one timestep
        
        It also handles scattering the particle off the boundaries, and keeps track of how much momentum is impulsed on the wall.
        
        This information will be important during our calculation of pressure.
        
        This method is updated in the VDW_System object to implement the velocity verlet algorithm in handling the interactions between
        particles.
        
        '''
        
        F = self.netForce()
        
        pcol = 0 #momentum change from collisions
        fake_system = VDW_System(dt = self.timestep, total_time = self.runtime, tau = self.temp, m = self.mass, 
                 N=len(self.particles), size = self.size,ic = self.interaction_const, b_0 = self.par_radius)
        accel = []
        dt = self.timestep
        
        #To implement the velocity verlet algorithm
        #we need to get the acceleration at the next step as well
        #This is computationally infefficient, but conserves energy.
        for i in range(len(self.particles)):
            
            particle = self.particles[i]
            force = F[i]
           
            a = force/particle.mass
            accel.append(a)
            
            pos_new = particle.pos+particle.vel*dt+0.5*a*dt**2
            
            fake_particle = ideal_ptr.Particle(mass=particle.mass,pos=pos_new)
            fake_system.particles[i] = fake_particle
        
        F_fake = fake_system.netForce()
        
        for j in range(len(self.particles)):
            
            #Now that we've calculated what the acceleration would be in the next step,
            #we implement the velocity verlet algorithm to update the particles' position
            #in our real system.
            
            par = self.particles[j]
            
            force_fake = F_fake[j]
            a_fake = F_fake[j]/fake_system.particles[j].mass
            
            v_new = self.particles[j].vel + (accel[j]+a_fake)/2*dt
            
            par.vel = v_new
            par.pos = fake_system.particles[j].pos
            
            check = self.in_box(par.pos) #check our particle is still in the box
            
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
        
        del(fake_system)
        
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
        
        