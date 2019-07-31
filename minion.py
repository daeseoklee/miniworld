from math import *
import random
vision=5
num_feature=5
merge_thres=0.3
mut_per_diff=1.6 #minimum number of mutations required for a species division, if every change is monotonic
mut_rate=merge_thres/mut_per_diff
idim=vision**2*num_feature
wdim=9*idim
scarcity=7.5
def translate(dna):
    return dna.get_c(),dna.get_brain()

class DNA():
    def __init__(self):
        pass
    def get_c(self):
        pass
    def get_brain(self):
        pass
    def mergeable(self,dna):
        pass
    def merge(self,dna):
        pass

def alter_c(x,a):
    if x+a>255:
        return 255
    elif x+a<0:
        return 0
    else:
        return int(x+a)
def difference_c(x,y):
    return abs(x[0]-y[0])+abs(x[1]-y[1])+abs(x[2]-y[2])
def difference_l(x,y):
    sum=0
    for i in range(len(x)):
        sum+=abs(x[i]-y[i])
    return sum

def mutate_c(c,sd):
    r,g,b=c
    return (alter_c(r,random.gauss(0,sd)),alter_c(g,random.gauss(0,sd)),alter_c(b,random.gauss(0,sd)))

def mutate_l(l,sd):
    for i in range(len(l)):
        l[i]+=random.gauss(0,sd)

class LinearDNA(DNA):
    def __init__(self,c,weights,init=False):
        super(LinearDNA,self).__init__()
        self.c=c
        self.brain=LinearBrain(weights=weights,init=init)
    def get_c(self):
        return self.c
    def get_brain(self):
        return self.brain
    def mergeable(self,dna):
        c_dist=difference_c(self.c,dna.c)
        w_dist=difference_l(self.brain.weights,dna.brain.weights)
        return random.uniform(0,1)<(1-erf(0.3*(c_dist/255)/(3*merge_thres)))*(1-erf(0.3*w_dist/(wdim*merge_thres)))
    def merge(self,dna):
        #return FoolDNA(self.c)
        r=random.choice([self.c[0],dna.c[0]])
        g=random.choice([self.c[1],dna.c[1]])
        b=random.choice([self.c[2],dna.c[2]])
        weights=[None for _ in range(wdim)]
        for i in range(9):
            weights[idim*i:idim*(i+1)]=random.choice([self.brain.weights[idim*i:idim*(i+1)],dna.brain.weights[idim*i:idim*(i+1)]])
        mutate_l(weights,mut_rate)
        return LinearDNA(mutate_c((r,g,b),255*mut_rate),weights,init=False)



class Brain():
    def __init__(self):
        pass
    def control(self,mi,snapshot,pos):
        pass

def in_range(x,y,a,b):
    return (0<=a<x) and (0<=b<y)
def softmax(l,c):
    sum=0
    for x in l:
        sum+=exp(c*x)
    return list(map(lambda x:exp(c*x)/sum,l))
def sigmoid(x,c):
    return 1/(1+exp(-c*x))

class LinearBrain(Brain):
    def __init__(self,weights,init=False):
        super(LinearBrain,self).__init__()
        if init:
            self.weights=[random.uniform(-1,1) for _ in range(wdim)]
        else:
            self.weights=weights
    def control(self,mi,snapshot,pos):
        xsize,ysize=len(snapshot),len(snapshot[0])
        input=[0.0 for _ in range(idim)]
        for i in range(-2,3):
            for j in range(-2,3):
                m,r,g,b,num=0,0,0,0,0
                for k in range(-mi.alen,1+mi.alen):
                    for l in range(-mi.alen,1+mi.alen):
                        x=pos[0]+i*(1+2*mi.alen)+k
                        y=pos[1]+j*(1+2*mi.alen)+l
                        c=snapshot[x%xsize][y%ysize]
                        if c==(255,255,255):
                            m+=1
                        elif c!=(0,0,0):
                            num+=1
                            r+=c[0]
                            g+=c[1]
                            b+=c[2]
                m/=(1+2*mi.alen)**2
                if num==0:
                    r,g,b=0.5,0.5,0.5
                else:
                    r,g,b=r/(255*num),g/(255*num),b/(255*num)
                num/=(1+2*mi.alen)**2
                input[num_feature * (5 * (2 + i) + (2 * j))] = scarcity*m
                input[num_feature*(5 * (2 + i) + (2 * j))+1] = r
                input[num_feature * (5 * (2 + i) + (2 * j)) + 2] = g
                input[num_feature * (5 * (2 + i) + (2 * j)) + 3] = b
                input[num_feature * (5 * (2 + i) + (2 * j)) + 4] = scarcity*num
        hidden=[0.0 for _ in range(9)]
        for i in range(9):
            for j in range(idim):
                hidden[i]+=self.weights[idim*i+j]*input[j]
        c=1/sqrt(idim*5/16) #variance of uniform(0,1)*uniform(-1,1) is 5/16 => sd of sum of (idim) independence these variables is sqrt((idim)*5/16)
        p=softmax(hidden[:4],c)
        u=random.uniform(0,1)
        if u<p[0]:
            mi.move_direc=0
        elif u<p[0]+p[1]:
            mi.move_direc=1
        elif u<p[0]+p[1]+p[2]:
            mi.move_direc=2
        else:
            mi.move_direc=3
        if random.uniform(0,1)<sigmoid(hidden[8],c):
            mi.move_dist=-1
        else:
            mi.move_dist=sigmoid(hidden[4+mi.move_direc],c)



class Minion():
    def __init__(self,dna):
        self.dna=dna
        self.id=None
        self.c,self.brain=translate(dna)
        self.alen=1
        self.mass=9
        self.stretched=False
        self.counter=0
        self.partner=None
        self.age=0
        self.move_dist=0
        self.move_direc=0.0
    def increase_age(self):
        self.age+=1
        if self.age==500:
            return False
        return True
    def count(self,partner):
        """
        if self.partner==partner:
            self.counter+=1
            if self.counter==2:
                self.counter=0
                self.partner=None
                return True
        else:
            self.partner=partner
            self.counter=1
            return False
        """
        return True

    def stretch(self):
        self.stretched=True
    def hold(self):
        self.stretched=False
    def take_mass(self,amount):
        self.mass+=amount
        self.alen=int((sqrt(self.mass)-1)/2)
    def loss_mass(self,amount):
        self.mass-=amount
        self.alen=int((sqrt(self.mass)-1)/2)
    def fatal(self):
        return self.alen<1
    def mergeable(self,mi):
        return self.alen>=2 and self.dna.mergeable(mi.dna)
    def get_child(self,mi):
        return Minion(self.dna.merge(mi.dna))
