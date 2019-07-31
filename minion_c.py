from math import *
import random

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

def alter(x,a):
    if x+a>255:
        return 255
    elif x+a<0:
        return 0
    else:
        return x+a
def difference(x,y):
    return abs(x[0]-y[0])+abs(x[1]-y[1])+abs(x[2]-y[2])

def mutate_c(c):
    r,g,b=c
    sd=10
    return (alter(r,random.gauss(0,sd)),alter(g,random.gauss(0,sd)),alter(b,random.gauss(0,sd)))

class FoolDNA(DNA):
    def __init__(self,c):
        super(FoolDNA,self).__init__()
        self.c=c
        self.brain=FoolBrain()
    def get_c(self):
        return self.c
    def get_brain(self):
        return self.brain
    def mergeable(self,dna):
        #return self.c==dna.c
        return random.uniform(0,1)>erf(difference(self.c,dna.c)/600)
    def merge(self,dna):
        #return FoolDNA(self.c)
        r=random.choice([self.c[0],dna.c[0]])
        g=random.choice([self.c[1],dna.c[1]])
        b=random.choice([self.c[2],dna.c[2]])
        #r=int((self.c[0]+dna.c[0])/2)
        #g = int((self.c[1] + dna.c[1]) / 2)
        #b = int((self.c[2] + dna.c[2]) / 2)
        return FoolDNA(mutate_c((r,g,b)))


class Brain():
    def __init__(self):
        pass
    def control(self,input):
        pass

class FoolBrain(Brain):
    def __init__(self):
        pass
    def control(self,mi,snapshot): #snapshot : whole map
        if random.uniform(0,1)<0.4:
            mi.move_dist=-1
        else:
            mi.move_dist = random.uniform(0, 1)
        mi.move_direc=random.randint(0,3)

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
        if self.age==200:
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
