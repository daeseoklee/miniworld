from math import *
import numpy as np
import copy
from fast_random_py import uniform,randint,randbool,gaussian
from fast_random_py import multiple_uniform,multiple_gaussian
from fast_random_py import randint_with_probs,seed
from random import sample,shuffle
seed()
num_out=12
vision_range=5
vision_resolution=3
vision_len=vision_range*vision_resolution
num_feature=5
rough_idim=vision_range**2*num_feature
rough_wdim=rough_idim*num_out
clear_idim=3*vision_len**2
linear_wdim=clear_idim*num_out

merge_thres=0.04
mut_per_diff=3.0 #minimum number of mutations required for a species division, if every change is monotonic
mut_rate=merge_thres/mut_per_diff

xsize=200
ysize=150
min_maxsize=3
max_maxsize=8
min_uptake=0.1
max_uptake=0.7
min_maxage=1500
max_maxage=1500

def set_size(x,y):
    global xsize
    global ysize
    xsize=x
    ysize=y

def set_traits_range(ms1,ms2,ut1,ut2,ma1,ma2):
    global min_maxsize,max_maxsize,min_uptake,max_uptake,min_maxage,max_maxage
    min_maxsize=ms1
    max_maxsize=ms2
    min_uptake=ut1
    max_uptake=ut2
    min_maxage=ma1
    max_maxage=ma2

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

def nearest_int(x):
    if x-int(x)<0.5:
        return int(x)
    return 1+int(x)

class ColorTrait():
    def __init__(self,r,g,b):
        self.r=r
        self.g=g
        self.b=b
    def opposite(self):
        return ColorTrait(255-self.r,255-self.g,255-self.b)
    def mutate(self,rate):
        dr=256*rate*gaussian()
        dg=256*rate*gaussian()
        db=256*rate*gaussian()
        self.r=nearest_int(self.r+dr)
        self.g=nearest_int(self.g+dg)
        self.b=nearest_int(self.b+db)
        if self.r>255:
            self.r=255
        elif self.r<0:
            self.r=0
        if self.g>255:
            self.g=255
        elif self.g<0:
            self.g=0
        if self.b>255:
            self.b=255
        elif self.b<0:
            self.b=0
    def mixed(t1,t2):
        if randbool():
            r=t1.r
        else:
            r=t2.r
        if randbool():
            g=t1.g
        else:
            g=t2.g
        if randbool():
            b=t1.b
        else:
            b=t2.b
        return ColorTrait(r,g,b)
    def normalized_difference(t1,t2):
        difference = abs(t1.r - t2.r) + abs(t1.g - t2.g) + abs(t1.b - t2.b)
        normalizer = 1 + max(t1.r, t2.r) + max(t1.g, t2.g) + max(t1.b, t2.b)
        return difference / normalizer
class PositiveTrait():
    def __init__(self,a):
        self.eps=0.001
        if a>self.eps:
            self.a=a
        else:
            self.a=self.eps
    def mutate(self,rate):
        self.a+=rate*self.a*gaussian()
        if self.a<self.eps:
            self.a=self.eps
    def mixed(t1,t2):
        if randbool():
            a=t1.a
        else:
            a=t2.a
        return PositiveTrait(a)
    def normalized_difference(t1,t2):
        difference=abs(t1.a-t2.a)
        normalizer=max(t1.a,t2.a)
        return difference/normalizer
class UnitTrait():
    def __init__(self,a):
        self.eps=0.001
        if a>self.eps:
            self.a=a
        elif a>1:
            self.a=1
        else:
            self.a=self.eps
    def mutate(self,rate):
        self.a+=rate*gaussian()
        if self.a<self.eps:
            self.a=self.eps
        elif self.a>1:
            self.a=1
    def mixed(t1,t2):
        if randbool():
            a=t1.a
        else:
            a=t2.a
        return UnitTrait(a)
    def normalized_difference(t1,t2):
        difference=abs(t1.a-t2.a)
        normalizer=max(t1.a,t2.a)
        return difference/normalizer
class FloatListTrait():
    def __init__(self,l,group_sizes):
        self.l=l
        self.group_sizes=group_sizes
    def mutate(self,rate):
        n=len(self.l)
        dl=multiple_gaussian(n)
        for i in range(n):
            self.l[i]+=rate*dl[i]
    def mixed(t1,t2):
        n=len(t1.l)
        m=len(t1.group_sizes)
        ini_index=0
        l=np.empty(n,dtype=float)
        for j in range(m):
            fin_index=ini_index+t1.group_sizes[j]
            if randbool():
                for i in range(ini_index,fin_index):
                    l[i]=t1.l[i]
            else:
                for i in range(ini_index,fin_index):
                    l[i]=t2.l[i]
            ini_index=fin_index
        return FloatListTrait(l,t1.group_sizes)
    def normalized_difference(t1,t2):
        n=len(t1.l)
        difference=0
        for i in range(n):
            difference+=abs(t1.l[i]-t2.l[i])
        normalizer=n
        return difference/normalizer

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

def randLinearDNA():
    r=randint(0,255)
    g=randint(0,255)
    b=randint(0,255)
    maxsize=randint(min_maxsize,max_maxsize)
    uptake=min_uptake+(max_uptake-min_uptake)*uniform()
    maxage=randint(min_maxage,max_maxage)
    weights=np.array(list(multiple_uniform(linear_wdim)),dtype=float)
    for i in range(linear_wdim):
        weights[i]=1-2*weights[i]
    group_sizes=np.empty(num_out,dtype=int)
    for i in range(num_out):
        group_sizes[i]=clear_idim
    colorTrait=ColorTrait(r,g,b)
    maxsizeTrait=PositiveTrait(maxsize)
    uptakeTrait=UnitTrait(uptake)
    maxageTrait=PositiveTrait(maxage)
    weightsTrait=FloatListTrait(weights,group_sizes)
    return LinearDNA(colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)


class LinearDNA(DNA): #DNA subclass
    def __init__(self,colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait):
        super(LinearDNA,self).__init__()
        self.colorTrait=colorTrait
        self.maxsizeTrait=maxsizeTrait
        self.uptakeTrait=uptakeTrait
        self.maxageTrait=maxageTrait
        self.weightsTrait=weightsTrait
        self.input_type="clear"
    def get_traits(self):
        return (self.colorTrait,self.maxsizeTrait,self.uptakeTrait,self.maxageTrait,self.weightsTrait)
    def translate(self):
        #DNA.translate() implementation
        return (self.colorTrait.r,self.colorTrait.g,self.colorTrait.b),\
               self.maxsizeTrait.a,self.uptakeTrait.a,self.maxageTrait.a
    def translate_weights(self):
        return self.weightsTrait.l
    def mergeable(self,dna):
        #DNA.mergeable() implementation
        num_traits=5
        const=0.11523  #erfinv(1-0.5**(1/num_traits))
        #--------------------------------------
        d1=ColorTrait.normalized_difference(self.colorTrait,dna.colorTrait)
        d2=PositiveTrait.normalized_difference(self.maxsizeTrait,dna.maxsizeTrait)
        d3=UnitTrait.normalized_difference(self.uptakeTrait,dna.uptakeTrait)
        d4=PositiveTrait.normalized_difference(self.maxageTrait,dna.maxageTrait)
        d5=FloatListTrait.normalized_difference(self.weightsTrait,dna.weightsTrait)
        p=1
        p*=1-erf(const*d1/merge_thres)
        p*=1-erf(const*d2/merge_thres)
        p*=1-erf(const*d3/merge_thres)
        p*=1-erf(const*d4/merge_thres)
        p*=1-erf(const*d5/merge_thres)
        return uniform()<p
    def merge(self,dna):
        #DNA.merge() implementation
        colorTrait=ColorTrait.mixed(self.colorTrait,dna.colorTrait)
        colorTrait.mutate(mut_rate)
        maxsizeTrait=PositiveTrait.mixed(self.maxsizeTrait,dna.maxsizeTrait)
        maxsizeTrait.mutate(mut_rate)
        uptakeTrait=UnitTrait.mixed(self.uptakeTrait,dna.uptakeTrait)
        uptakeTrait.mutate(mut_rate)
        maxageTrait=PositiveTrait.mixed(self.maxageTrait,dna.maxageTrait)
        maxageTrait.mutate(mut_rate)
        weightsTrait=FloatListTrait.mixed(self.weightsTrait,dna.weightsTrait)
        weightsTrait.mutate(mut_rate)
        return LinearDNA(colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)




class LinearBrain(Brain): #Brain subclass
    def __init__(self,weights,init=False):
        super(LinearBrain,self).__init__()
        if init:
            weights=1-2*np.random.rand(linear_wdim)
        self.weights=weights
    def control(self,mi,inp):
        #Brain.control() implementation
        if mi.frozen:
            return
        hidden=np.matmul(inp,self.weights.reshape(clear_idim,-1))
        c=1/sqrt(clear_idim*5/16) #variance of uniform(0,1)*uniform(-1,1) is 5/16 => sd of sum of (clear_idim) independence these variables is sqrt((clear_idim)*5/16)
        direc_probs=softmax(hidden[:4],c)
        action_probs=softmax(hidden[8:12],c)
        mi.move_direc=randint_with_probs(0,3,direc_probs)
        mi.move_dist=sigmoid(hidden[4+mi.move_direc],c)
        mi.action=randint_with_probs(0,3,action_probs)


def construct_minion(dna,alen,pos,do_freeze):
    mi=Minion(dna)
    mi.take_mass((1+2*alen)**2-9)
    if pos==(-1,-1):
        pos=(randint(0,xsize-1),randint(0,ysize-1))
    mi.pos=pos
    if do_freeze:
        mi.freeze()
    return mi

class Minion():
    def __init__(self,dna):
        self.dna=dna #DNA
        self.id=None #int
        self.pos=(0,0)
        self.dead=False
        self.color,self.maxsize,self.uptake,self.maxage=dna.translate()
        self.brain=LinearBrain(self.dna.translate_weights(),init=False)
        self.idim=clear_idim
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






def dist(x,y,T):
    #int*int*int->int
    #use: distance on a periodic line of period 'T'
    #called in: World.huntable()
    if x>y:
        return min(x-y,T+y-x)
    else:
        return min(y-x,T+x-y)
class World():
    def __init__(self,xsize,ysize,total_mass,mis,no_age=False,no_birth=False,no_eat=False,\
                 no_energy=False, no_excrete=False,\
                 no_hunt=False,halluc=False,record_pedigree=False):
        self.snapshot=np.array([[0 for _ in range(ysize)] for _ in range(xsize)])
        self.xsize=xsize #int
        self.ysize=ysize #int
        self.moment=0
        self.new_id=0
        self.mins=[[0 for _ in range(ysize)] for _ in range(xsize)]
        nut=total_mass
        for mi in mis:
            nut-=mi.mass
        for _ in range(nut):
            i=randint(0,xsize-1)
            j=randint(0,ysize-1)
            self.mins[i][j]+=1
        self.pedigree=[] #(int,int)[]
        self.occupy_map=[[[] for _ in range(self.ysize)] for _ in range(self.xsize)] #minion.Minion[][][]
        self.mis=[] #minion.Minion[]
        for mi in mis:
            self.register(mi)
        #additional attributes for test
        self.no_age=no_age
        self.no_birth=no_birth #bool
        self.no_eat=no_eat #bool
        self.no_energy=no_energy #bool
        self.no_excrete=no_excrete #bool
        self.no_hunt=no_hunt #bool
        self.messiness=0
        self.halluc=halluc #bool
        self.hidden_mass=0
        #measurement option
        self.record_pedigree=record_pedigree
    def get_mins(self):
        return self.mins
    def get_moment(self):
        return self.moment
    def get_xsize(self):
        return self.xsize
    def get_ysize(self):
        return self.ysize
    #---------------------------------------------

    def register(self,mi):
        #World*minion.Minion*(int,int)->()
        #use: register 'mi' in the world and locate it at 'pos'
        #called in: World.__init__(), World.childbirth()
        self.mis.append(mi) #don't change this order, because it's relevant in some tests
        mi.id=self.new_id
        self.new_id+=1
        for pos in self.bodyposs(mi):
            self.occupy_map[pos[0]][pos[1]].append(mi)
    def unregister_deads(self):
        #World->()
        #use: remove all the data associated with every mi in self.mis for which mi.dead==True. called after each execution of 1 moment.
        #called in: World.one_step()
        i=0
        while i<len(self.mis):
            mi=self.mis[i]
            if mi.dead:
                self.mis.pop(i)
                del mi
                continue
            i+=1
    #position----------------------------------
    def on(self,pos):
        #World*(int,int)->minion.Minion[]
        #use: list of minions whose body covers the position
        return self.occupy_map[pos[0]][pos[1]]
    def bodyposs(self,mi):
        #World*minion.Minion->generator[(int,int)]
        #use: yield positions that are covered by the body of 'mi'
        #called in: World.mk_corpse(),World.kill(),World.move(),World.eat() etc
        x,y=mi.pos
        for i in range(-mi.alen,1+mi.alen):
            for j in range(-mi.alen,1+mi.alen):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)
    def skewbodyposs(self,mi):
        #World*minion.Minion->generator[(int,int)]
        #use: Similar to bodyposs, but only yield a proper subset which is sufficient for detecting all the other minions whose body is contained in that of 'mi'; more efficient for those cases
        #called in: World.try_hunt()
        x,y=mi.pos
        for i in range(-mi.alen,1+mi.alen,3):
            for j in range(-mi.alen,1+mi.alen,3):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)

    def genipos(self,mi):
        #World*minion.Minion->(int,int)
        #use: determine the position of the reproductive organ
        #called in: World.stretch()
        x,y=mi.pos
        if mi.move_direc==0:
            return ((x+2*mi.alen)%self.xsize,y)
        elif mi.move_direc==1:
            return (x,(y+2*mi.alen)%self.ysize)
        elif mi.move_direc==2:
            return ((x-2*mi.alen)%self.xsize,y)
        elif mi.move_direc==3:
            return (x,(y-2*mi.alen)%self.ysize)

    def plusminusposs(self,alen_before,mi):
        #World*int*minion.Minion->generator[(int,int)]
        #use: helper function for self.occupy_map adjustion. Assuming the size of 'mi' has changed from alen_before to mi.alen, determine the positions to add/substract.
        #called in: World.take_mass(),World.loss_mass()
        x,y=mi.pos
        if mi.alen>alen_before:
            shorter,longer=alen_before,mi.alen
        else:
            shorter,longer=mi.alen,alen_before
        for i in range(-longer,1+longer):
            for j in range(-longer,-shorter):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)
            for j in range(1+shorter,1+longer):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)
        for j in range(-shorter,1+shorter):
            for i in range(-longer,-shorter):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)
            for i in range(1+shorter,1+longer):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)

    #------------------------------------
    def take_mass(self,mi,amount):
        #World*minion.Minion*float->()
        #use: call mi.take_mass() and adjust self.occupy_map accordingly
        #called in: World.digest()
        alen_before=mi.alen
        mi.take_mass(amount)
        if mi.alen>alen_before:
            for pos in self.plusminusposs(alen_before,mi):
                self.occupy_map[pos[0]][pos[1]].append(mi)
    def loss_mass(self,mi,amount):
        #World*minion.Minion*float->()
        #use: call mi.loss_mass() and adjust self.occupy_map accordingly
        #called in: World.childbirth()
        """
        print("lost")
        print("alen:",mi.alen,"mass:",mi.mass,"amount:",amount,"pos:",mi.pos)
        for pos in self.bodyposs(mi):
            if not mi in self.occupy_map[pos[0]][pos[1]]:
                raise Exception("here1")
        print("no problem here1")
        """
        alen_before=mi.alen
        mi.loss_mass(amount)
        if mi.alen<alen_before:
            for pos in self.plusminusposs(alen_before,mi):
                self.occupy_map[pos[0]][pos[1]].remove(mi)
    def take_energy(self,mi,amount):
        #World*minion.Minion*float->()
        #use: energy uptake
        #called in: World.digest()
        if amount<=mi.max_energy-mi.energy:
            mi.energy+=amount
        else:
            mi.energy=mi.max_energy
    def loss_energy(self,mi,amount,fromm=None):
        #World*minion.Minion*float(*string)->bool
        #use: energy consumption and return whether 'mi' is still alive
        #called in: World.childbirth(), World.stretch(), World.move(), World.basal_metabolism()
        if self.no_energy:
            return True
        mi.energy-=amount
        if mi.energy<=0:
            self.kill(mi)
            #print("death by energy")
            #print("from:",fromm)
            #print("body mass:",mi.mass)
            return False
        return True
    def childbirth(self,girl,boy):
        #World*minion.Minion*minion.Minion->()
        #use: childbirth and consequences
        #called in: World.stretch()
        if self.no_birth:
            return
        if girl.mergeable(boy):
            #print("birth!!!")
            child=girl.get_child(boy)
            pos=girl.pos
            if girl.move_direc==0:
               x,y=(pos[0]-girl.alen-3)%self.xsize,pos[1]
            elif girl.move_direc==1:
                x,y=pos[0],(pos[1]-girl.alen-3)%self.ysize
            elif girl.move_direc==2:
                x,y=(pos[0]+girl.alen+3)%self.xsize,pos[1]
            elif girl.move_direc==3:
                x,y=pos[0],(pos[1]+girl.alen+3)%self.ysize
            self.loss_mass(girl,child.mass)
            alive=self.loss_energy(girl,girl.birth_consum_rate*girl.mass,fromm="birth")
            if alive and girl.alen<1:
                self.kill(girl)
                #print("death by childbirth")
            child.pos=(x,y)
            self.register(child)
            if self.record_pedigree:
                self.pedigree.append((girl.id,boy.id,child.id))

    def huntable(self,pred,pray):
        #World*minion.Minion*minion.Minion->bool
        #use: determine whether 'pred' can hunt 'pray'
        #called in: World.try_hunt()
        return (not self.no_hunt) and max(dist(pred.pos[0],pray.pos[0],self.xsize),dist(pred.pos[1],pray.pos[1],self.ysize))<=pred.alen-pray.alen and pred!=pray
    def mk_corpse(self,mi):
        #World*minion.Minion->()
        #use: leave minerals on the occupying area of 'mi'
        #called in: World.kill()
        area=(1+2*mi.alen)**2
        q=mi.mass/area
        for pos in self.bodyposs(mi):
            self.mins[pos[0]][pos[1]]+=q

    def kill(self,mi,corpse=True):
        #World*minion.Minion(*bool)->()
        #use: Make 'mi' in a to-be-removed state, by setting mi.dead True, making corpse, and removing 'mi' from 'self' attributes. The remaining jobs will be done by World.unregister_deads(). The corpse is not made if corpse==False.
        #called in: World.lose_energy(),World.childbirth(), World.hunt(),World.kill_elderly()
        if mi.dead:
            #print("doublekill")
            return
        if corpse:
            self.mk_corpse(mi)
        for pos in self.bodyposs(mi):
            self.occupy_map[pos[0]][pos[1]].remove(mi)
        mi.dead=True
    def kill_elderly(self):
        #World->()
        #use: increase ages and kill those that attained the maximum age
        #called in: World.one_step()
        if self.no_age:
            return
        for mi in self.mis:
            if not mi.increase_age():
                self.kill(mi)
                #print("death by age")
    #control and action---------------------------------------
    def control_all(self):
        #World->()
        #use: for each 'mi' in self.mis(), determine mi.move_direc, mi.move_dist, and mi.action based on their relation to the world.
        #called in: World.one_step()
        if not self.halluc:
            for mi in self.mis:
                inp=mi.get_input(self.snapshot,mi.pos)
                mi.brain.control(mi,inp)
        else:
            inps=[mi.get_input(self.snapshot,mi.pos) for mi in self.mis]
            shuffle(inps)
            for i in range(len(self.mis)):
                mi=self.mis[i]
                inp=inps[i]
                mi.brain.control(mi,inp)

    def move(self,mi):
        #World*minion.Minion->()
        #use: perform movement, based on mi.move_direc and mi.move_dist
        #called in: World.act()
        move_direc=mi.move_direc
        move_dist = int(2 * mi.move_dist * sqrt(mi.mass))
        #print("dist:",move_dist)
        if move_dist==0:
            return
        x,y=mi.pos
        a,b=0,0
        if move_direc==0:
            a+=move_dist
        elif move_direc==1:
            b+=move_dist
        elif move_direc==2:
            a-=move_dist
        elif move_direc==3:
            b-=move_dist
        for pos in self.bodyposs(mi):
            self.occupy_map[pos[0]][pos[1]].remove(mi)
        mi.pos=((x+a)%self.xsize,(y+b)%self.ysize)
        for pos in self.bodyposs(mi):
            self.occupy_map[pos[0]][pos[1]].append(mi)
        self.loss_energy(mi,mi.move_consum_rate*mi.mass*move_dist,fromm="move")
        mi.cum_dist+=move_dist
    def stretch(self,mi):
        #World*minion.Minion->()
        #use: mating trial, by stretching its reproductive organ
        #called in: World.act()
        alive=self.loss_energy(mi,mi.stretch_consum_rate*mi.mass,fromm="stretch")
        if not alive:
            return
        genipos=self.genipos(mi)
        l=copy.copy(self.occupy_map[genipos[0]][genipos[1]])
        for girl in l:
            self.childbirth(girl,mi)
            self.messiness+=1
            alive=self.loss_energy(mi,mi.sex_consum_rate*mi.mass,fromm="sex")
            if not alive:
                return

    def excrete(self,mi,amount):
        #World*minion.Minion*float->()
        #use: excretion
        #called in: World.digest()
        if self.no_excrete:
            self.hidden_mass+=amount
            return
        a,b,dist=0,0,4*mi.alen
        area=(1+2*mi.alen)**2
        q,r=amount//area,amount%area
        if mi.move_direc==0:
            a-=dist
        elif mi.move_direc==1:
            b-=dist
        elif mi.move_direc==2:
            a+=dist
        elif mi.move_direc==3:
            b+=dist
        poss=self.bodyposs(mi)
        if q==0:
            for _ in range(int(r)):
                pos=next(poss)
                self.mins[(pos[0]+a)%self.xsize][(pos[1]+b)%self.ysize]+=1
            pos=next(poss)
            self.mins[(pos[0]+a)%self.xsize][(pos[1]+b)%self.ysize]+=r-int(r)
        else:
            for i,pos in enumerate(poss):
                if i<=r-1:
                    self.mins[(pos[0]+a)%self.xsize][(pos[1]+b)%self.ysize]+=q+1
                elif r-1<i<r:
                    self.mins[(pos[0]+a)%self.xsize][(pos[1]+b)%self.ysize]+=q+r-int(r)
                else:
                    self.mins[(pos[0]+a)%self.xsize][(pos[1]+b)%self.ysize]+=q


    def digest(self,mi,amount):
        #World*minion.Minion*float->()
        #use: digestion
        #called in: World.hunt(), World.eat()
        take=min(mi.uptake*amount,(1+2*mi.maxsize)**2-mi.mass)
        out=amount-take
        self.take_mass(mi,take)
        self.excrete(mi,out)
        self.take_energy(mi,out)

    def hunt(self,pred,pray):
        #World*minion.Minion*minion.Minion->()
        #use: predation
        #called in: World.try_hunt()
        self.digest(pred,pray.mass)
        self.kill(pray,corpse=False)
        #print("death by subsumption")
    def try_hunt(self,mi):
        #World*minion.Minion->()
        #use: predation trial
        #called in: World.act()
        to_try=[]
        for pos in self.skewbodyposs(mi):
            for pray in self.on(pos):
                if not pray in to_try:
                    to_try.append(pray)
        for pray in to_try:
            if self.huntable(mi,pray):
                self.hunt(mi,pray)
    def eat(self,mi):
        #World*minion.Minion->()
        #use: eat minerals
        #called in: World.act()
        if self.no_eat:
            return
        total=0
        for pos in self.bodyposs(mi):
            total+=self.mins[pos[0]][pos[1]]
            self.mins[pos[0]][pos[1]]=0
        self.digest(mi,total)
    def act(self,mi):
        #World*minion.Minion->()
        #use: perform action based on mi.action
        #called in: World.act_all(0
        if mi.frozen:
            return
        #print("action:",mi.action)
        if mi.action==0:
            self.move(mi)
        elif mi.action==1:
            self.stretch(mi)
        elif mi.action==2:
            self.try_hunt(mi)
        elif mi.action==3:
            self.eat(mi)
    def basal_metabolism(self,mi):
        #World*minion.Minion->()
        #use: basal metabolism
        #called in: World.basal_metabolism_all()
        self.loss_energy(mi,mi.basal_metabolic_rate*mi.mass,fromm="basal")

    #outmost-----------------------------------------

    def act_all(self):
        #World->()
        #use: for each 'mi' in self.mis peform mi.action, but for newly appended ones during the execution.
        #called in: World.one_step()
        for i in range(len(self.mis)): #don't use for mi in self.mis
            mi=self.mis[i]
            if not mi.dead:
                self.act(mi)
    def basal_metabolism_all(self):
        #World->()
        #use: for each 'mi' in self.mis perform mi.bassal_metabolism()
        #called in: World.one_step()
        for i in range(len(self.mis)):
            mi=self.mis[i]
            if not mi.dead:
                self.basal_metabolism(mi)
    def render(self):
        #World->()
        #use: save the capture image of the world in self.snapshot
        #called in: World.one_step()
        for i in range(self.xsize):
            for j in range(self.ysize):
                self.snapshot[i][j]=0 #black
        for mi in self.mis:
            c=mi.color
            for pos in self.bodyposs(mi):
                self.snapshot[pos[0]][pos[1]]=256**2*c[0]+256*c[1]+c[2]

        for i in range(self.xsize):
            for j in range(self.ysize):
                if self.mins[i][j]>0:
                    self.snapshot[i][j]=256**3-1 #white
    #public API-------------------------------------------
    #main
    def one_step(self):
        #World->()
        #use: do all the things that should happen in per-moment basis
        #called in: play.py
        #print("mass:",self.total_mass())
        self.control_all()
        self.act_all()
        self.basal_metabolism_all()
        self.moment+=1
        self.kill_elderly()
        self.unregister_deads()
        self.render()
    #basic
    def get_xsize(self):
        return self.xsize
    def get_ysize(self):
        return self.ysize
    def total_mass(self):
        #World->float
        #use: total mass in the world
        total=0
        for mi in self.mis:
            total+=mi.mass
        for i in range(self.xsize):
            for j in range(self.ysize):
                total+=self.mins[i][j]
        return total+self.hidden_mass
    def total_min(self):
        #World->float
        #use: total number of minerals in the world
        total=0
        for i in range(self.xsize):
            for j in range(self.ysize):
                total+=self.mins[i][j]
        return total
    #end conditions
    def limit(self,n):
        return self.moment>=n
    def exhausted(self):
        #World->bool
        #use: whether there's no mineral left
        #called in: play.py
        for i in range(self.xsize):
            for j in range(self.ysize):
                if self.mins[i][j]!=0:
                    return False
        return True
    def extincted(self):
        return self.mis==[]
    #measurements
    def get_pedigree(self):
        return self.pedigree
    def sample_dna(self,n):
        if n>=len(self.mis):
            subpopulation=sample(self.mis,n)
        else:
            subpopulation=self.mis
        dnas=[]
        for mi in subpopulation:
            dnas.append(mi.dna)
        return dnas
    def get_mis(self):
        return self.mis
    def get_population(self):
        return len(self.mis)
    def get_colors(self):
        for mi in self.mis:
            yield mi.color
    def get_avg_r(self):
        if self.mis==[]:
            return 0
        sum=0
        for mi in self.mis:
            sum+=mi.color[0]
        return sum/len(self.mis)
    def get_avg_g(self):
        if self.mis==[]:
            return 0
        sum=0
        for mi in self.mis:
            sum+=mi.color[1]
        return sum/len(self.mis)
    def get_avg_b(self):
        if self.mis==[]:
            return 0
        sum=0
        for mi in self.mis:
            sum+=mi.color[2]
        return sum/len(self.mis)

    #for use in visualization
    def body_rect(self,k):
        #called in: play.py
        for mi in self.mis:
            x,y=mi.pos
            yield mi.color,(k*(x-mi.alen),k*(y-mi.alen),k*(1+2*mi.alen),k*(1+2*mi.alen))
    def geni_rect(self,k):
        #called in: play.py
        for mi in self.mis:
            if mi.action==1:
                x,y=mi.pos
                if mi.move_direc==0:
                    yield mi.color,(k*x,k*y,k*(1+2*mi.alen),k)
                elif mi.move_direc==1:
                    yield mi.color,(k*x,k*y,k,k*(1+2*mi.alen))
                elif mi.move_direc==2:
                    yield mi.color,(k*x+k-1,k*y,-k*(1+2*mi.alen),k)
                elif mi.move_direc==3:
                    yield mi.color,(k*x,k*y+k-1,k,-k*(1+2*mi.alen))