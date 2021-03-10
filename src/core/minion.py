from math import *
import random
import util
import numpy as np

vision_range=5
vision_resolution=3
vision_len=vision_range*vision_resolution
num_feature=5
rough_idim=vision_range**2*num_feature
rough_wdim=rough_idim*12
clear_idim=3*vision_len**2
linear_wdim=clear_idim*12

merge_thres=0.04
mut_per_diff=3.0 #minimum number of mutations required for a species division, if every change is monotonic
mut_rate=merge_thres/mut_per_diff


def set_vision(new_vision_range,new_vision_resolution):
    global vision_range
    global vision_resolution
    vision_range=new_vision_range
    vision_resolution=new_vision_resolution

def set_heredity(new_merge_thres,new_mut_per_diff):
    global merge_thres
    global mut_per_diff
    merge_thres=new_merge_thres
    mut_per_diff=new_mut_per_diff


class DNA():
    def __init__(self):
        pass
    def mergeable(self,dna):
        #DNA*DNA->bool
        #use: determine whether they can make a child, in terms of DNA compatibility
        #called in: Minion.mergeable()
        pass
    def merge(self,dna):
        #DNA*DNA->DNA
        #use: determine the DNA of the child, when 'self' is from mother's and 'dna' is from father's.
        #called in: Minion.get_child()
        pass
    def translate(self):
        #DNA-> ~
        #use: determine the Minion attributes based on 'self'
        #called in: Minion.__init__()
        pass

class Brain():
    def __init__(self):
        pass
    def control(self,mi,snapshot,pos):
        #Brain*Minion*np.ndarray(dim 2)*(int,int)->()
        #use: determine mi.move_direc,mi.move_dist, and mi.action based on 'self','snapshot' and 'pos'.
        #called in: world.World.control_all()
        #effect is used in: world.World.act(),play.draw()
        pass


def difference_color(x,y):
    #(int,int,int)*(int,int,int)->int
    #called in: LinearBrain.mergeable()
    return abs(x[0]-y[0])+abs(x[1]-y[1])+abs(x[2]-y[2])
def difference_list(x,y):
    #float[]*float[]->float
    #called in: LinearBrain.mergeable()
    sum=0
    for i in range(len(x)):
        sum+=abs(x[i]-y[i])
    return sum
def nearest_int(x):
    if x-int(x)<0.5:
        return int(x)
    return 1+int(x)
def alter_color(x,a):
    #int*float->int
    #called in: mutate_color()
    if x+a>255:
        return 255
    elif x+a<0:
        return 0
    else:
        return nearest_int(x+a)
def mutate_color(c,sd):
    #(int,int,int)*float->(int,int,int)
    #called in: LinearBrain.merge()
    r,g,b=c
    return (alter_color(r,random.gauss(0,sd)),alter_color(g,random.gauss(0,sd)),alter_color(b,random.gauss(0,sd)))


def mutate_list(l,sd):
    #float[]*float->()
    #called in: LinearBrain.merge()
    for i in range(len(l)):
        l[i]+=np.random.normal(0,sd)

def mutate_positive(a,sd):
    #float*float->float
    #called in: LinearBrain.merge()
    eps=0.001
    b=random.gauss(a,sd)
    if b>=eps:
        return b
    return eps

def mutate_unit(a,sd): #unit inerval
    #float*float->float
    #called in: LinearBrain.merge()
    b=random.gauss(a,sd)
    if b>1.0:
        return 1.0
    elif b<0.0:
        return 0.0
    else:
        return b

def softmax(l,c):
    #np.ndarray(dim 1)*float->np.ndarray(dim 1)
    #called in: LinearBrain.control()
    m=np.exp(c*l)
    return m/m.sum()

def sigmoid(x,c):
    #float*float->float
    #called in: LinearBrain.control()
    return 1/(1+exp(-c*x))



#---------------------------------------------------------------

#--------------------------------------------------------------

class LinearDNA(DNA): #DNA subclass
    def __init__(self,color,maxsize,uptake,maxage,weights,init=False):
        super(LinearDNA,self).__init__()
        self.color=color
        self.maxsize=maxsize
        self.uptake=uptake
        self.maxage=maxage
        self.brain=LinearBrain(weights=np.array(weights),init=init)
        self.input_type="clear"
    def translate(self):
        #DNA.translate() implementation
        return self.color,self.maxsize,self.uptake,self.maxage,self.brain
    def mergeable(self,dna):
        #DNA.mergeable() implementation
        num_traits=5
        const=0.11523  #erfinv(1-0.5**(1/num_traits))
        #color trait--------------------------------------
        color_dist=difference_color(self.color,dna.color)
        color_range=1+max(self.color[0],dna.color[0])+max(self.color[1],dna.color[1])+max(self.color[2],dna.color[2])
        #maxsize trait-------------------------------------
        maxsize_dist=abs(self.maxsize-dna.maxsize)
        maxsize_range=max(self.maxsize,dna.maxsize)
        #uptake trait
        uptake_dist=abs(self.uptake-dna.uptake)
        uptake_range=max(self.uptake,dna.uptake)
        #maxage trait-------------------------------------
        maxage_dist=abs(self.maxage-dna.maxage)
        maxage_range=max(self.maxage,dna.maxage)
        #weights trait----------------------------------
        weights_dist=difference_list(self.brain.weights,dna.brain.weights)
        weights_range=len(self.brain.weights)
        return util.with_prob((1-erf((const*color_dist)/(merge_thres*color_range)))*(1-erf((const*maxsize_dist)/(merge_thres*maxsize_range)))*(1-erf((const*uptake_dist)/(merge_thres*uptake_range)))*(1-erf((const*maxage_dist)/(merge_thres*maxage_range))))*(1-erf((const*weights_dist)/(merge_thres*weights_range)))
    def merge(self,dna):
        #DNA.merge() implementation
        #color trait--------------------------------------
        r=random.choice([self.color[0],dna.color[0]])
        g=random.choice([self.color[1],dna.color[1]])
        b=random.choice([self.color[2],dna.color[2]])
        new_color=mutate_color((r,g,b),255*mut_rate)
        #maxsize trait-------------------------------------
        new_maxsize=random.choice([self.maxsize,dna.maxsize])
        new_maxsize=mutate_positive(new_maxsize,mut_rate*new_maxsize)
        #uptake trait-------------------------------------
        new_uptake=random.choice([self.uptake,dna.uptake])
        new_uptake=mutate_unit(new_uptake,mut_rate*new_uptake)
        #maxage trait-------------------------------------
        new_maxage=random.choice([self.maxage,dna.maxage])
        new_maxage=mutate_positive(new_maxage,mut_rate*new_maxage)
        #weights trait----------------------------------
        new_weights=np.zeros(linear_wdim)
        for i in range(12):
            new_weights[clear_idim*i:clear_idim*(i+1)]=random.choice([self.brain.weights[clear_idim*i:clear_idim*(i+1)],dna.brain.weights[clear_idim*i:clear_idim*(i+1)]])
        mutate_list(new_weights,mut_rate)
        return LinearDNA(new_color,new_maxsize,new_uptake,new_maxage,new_weights,init=False)




class LinearBrain(Brain): #Brain subclass
    def __init__(self,weights,init=False):
        super(LinearBrain,self).__init__()
        if init:
            self.weights=np.array([random.uniform(-1,1) for _ in range(linear_wdim)])
        else:
            self.weights=np.array(weights)
    def control(self,mi,inp):
        #Brain.control() implementation
        if mi.frozen:
            return
        hidden=np.matmul(inp,self.weights.reshape(clear_idim,-1))
        c=1/sqrt(clear_idim*5/16) #variance of uniform(0,1)*uniform(-1,1) is 5/16 => sd of sum of (clear_idim) independence these variables is sqrt((clear_idim)*5/16)
        direc_probs=softmax(hidden[:4],c)
        action_probs=softmax(hidden[8:12],c)
        mi.move_direc=util.rand_with_probs(direc_probs)
        mi.move_dist=sigmoid(hidden[4+mi.move_direc],c)
        mi.action=util.rand_with_probs(action_probs)




class Minion():
    def __init__(self,dna):
        self.dna=dna #DNA
        self.id=None #int
        self.dead=False
        self.color,self.maxsize,self.uptake,self.maxage,self.brain=dna.translate()
        self.alen=1
        self.mass=9
        self.energy=5.0*self.mass
        self.adjust_energy_from_mass()
        self.avg_consum_rate=0.025
        self.basal_metabolic_rate=self.avg_consum_rate/4 #this*self.mass per moment
        self.move_consum_rate=self.avg_consum_rate/40
        self.stretch_consum_rate=self.avg_consum_rate/4
        self.sex_consum_rate=self.avg_consum_rate/4
        self.birth_consum_rate=0.3
        self.counter=0
        self.partner=None #None or Minion
        self.age=0
        self.action=0
        self.move_dist=0.0 #don't alter this default value, it is used for freezing
        self.move_direc=0
        self.frozen=False
        self.cum_dist=0
    def freeze(self):
        #Minion->()
        #called in: world.World.__init__()
        #effect is used in: world.World.act()
        self.move_dist=0.0
        self.move_direc=0
        self.frozen=True
    def defreeze(self):
        #Minion->()
        #effect is used in: world.World.act()
        self.frozen=False
    def increase_age(self):
        #Minion->bool->()
        #use: self.age++ and return whether the maximum age is attained
        #called in: world.World.kill_elderly()
        self.age+=1
        if self.age>self.maxage:
            return False
        return True
    def adjust_energy_from_mass(self):
        self.max_energy=10.0*self.mass
        if self.energy>self.max_energy:
            self.energy=self.max_energy
    def take_mass(self,amount):
        #Minion*float->()
        #use: increase mass and adjust size accordingly
        #called in: world.World.take_mass()
        self.mass+=amount
        self.adjust_energy_from_mass()
        self.alen=int((sqrt(self.mass)-1)/2)
    def loss_mass(self,amount):
        #Minion*float->()
        #use: reduce mass and adjust size accordingly
        #called in: world.World.loss_mass()
        self.mass-=amount
        self.adjust_energy_from_mass()
        self.alen=int((sqrt(self.mass)-1)/2)
    def mergeable(self,mi):
        #Minion*Minion->bool
        #use: determine whether the act of reproduction was successful
        #called in: world.World.childbirth()
        return self.alen>=2 and self.dna.mergeable(mi.dna)
    def get_child(self,mi):
        #Minion*Minion->Minion
        #use: construct and return the child out of mother and father
        #called in: world.World.childbirth()
        return Minion(self.dna.merge(mi.dna))

    def get_input(self,snapshot,pos):
        #Minion*np.ndarray(dim 2)*(int,int)->np.ndarray(dim 3)
        #use: make subjective input data based on 'snapshot' from world.World.snapshot
        #called in: world.World.control_all()
        if self.dna.input_type=="clear":
            xsize, ysize = len(snapshot), len(snapshot[0])
            x_min = pos[0] - (1 + 2 * self.alen) * ((vision_range - 1) // 2) - self.alen
            x_max = pos[0] + (1 + 2 * self.alen) * ((vision_range - 1) // 2) + self.alen
            y_min = pos[1] - (1 + 2 * self.alen) * ((vision_range - 1) // 2) - self.alen
            y_max = pos[1] + (1 + 2 * self.alen) * ((vision_range - 1) // 2) + self.alen
            x_indices = np.linspace(x_min, x_max, vision_len, dtype=int) % xsize
            y_indices = np.linspace(y_min, y_max, vision_len, dtype=int) % ysize
            sub = snapshot[np.ix_(x_indices, y_indices)].reshape(-1)
            r_inp = (sub // 256 ** 2) / 256
            g_inp = ((sub % 256 ** 2) // 256) / 256
            b_inp = (sub % 256) / 256
            inp = np.concatenate([r_inp, g_inp, b_inp])
            return inp