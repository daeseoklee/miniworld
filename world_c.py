import minion
import random
from math import *

class World():
    def __init__(self,xsize,ysize,total_mass,mis):
        self.snapshot=[[(0,0,0) for _ in range(ysize)] for _ in range(xsize)]
        self.xsize=xsize
        self.ysize=ysize
        self.moment=0
        self.new_id=0
        self.total_mass=total_mass
        self.mins=dict()
        for i in range(xsize):
            for j in range(ysize):
                self.mins[(i,j)]=0
        nut=total_mass
        for mi in mis:
            nut-=mi.mass
        for _ in range(nut):
            i=random.randint(0,xsize-1)
            j=random.randint(0,ysize-1)
            self.mins[(i,j)]+=1
        self.mis=[]
        self.pedigree=[]
        self.poss=dict()
        for mi in mis:
            xpos=random.randint(0,xsize-1)
            ypos=random.randint(0,ysize-1)
            self.register(mi,(xpos,ypos))
    def in_range(self,pos):
        return (pos[0] in range(self.xsize)) and (pos[1] in range(self.ysize))
    def in_body(self,mi,pos):
        x,y=self.poss[mi]
        return (x-mi.alen<=pos[0]<=x+mi.alen) and (y-mi.alen<=pos[1]<=y+mi.alen)
    def register(self,mi,pos):
        self.mis.append(mi)
        mi.id=self.new_id
        self.new_id+=1
        self.poss[mi]=pos

    def sexy(self,boy,girl):
        x,y=self.poss[boy]
        if boy.move_direc==0:
            return self.in_body(girl,(x+2*boy.alen,y))
        elif boy.move_direc==1:
            return self.in_body(girl,(x,y+2*boy.alen))
        elif boy.move_direc==2:
            return self.in_body(girl,(x-2*boy.alen,y))
        elif boy.move_direc==3:
            return self.in_body(girl,(x,y-2*boy.alen))
    def childbirth(self,girl:minion.Minion,boy:minion.Minion):
        if girl.mergeable(boy):
            child:minion.Minion=girl.get_child(boy)
            pos=self.poss[girl]
            if girl.move_direc==0:
               x,y=max(0,pos[0]-girl.alen-3),pos[1]
            elif girl.move_direc==1:
                x,y=pos[0],max(0,pos[1]-girl.alen-3)
            elif girl.move_direc==2:
                x,y=min(self.xsize-1,pos[0]+girl.alen+3),pos[1]
            elif girl.move_direc==3:
                x,y=pos[0],min(self.ysize-1,pos[1]+girl.alen+3)
            girl.loss_mass(child.mass)
            if girl.fatal():
                self.kill(girl)
            self.register(child,(x,y))
            self.pedigree.append((girl.id,boy.id,child.id))

    def subsumable(self,pred:minion.Minion,pray:minion.Minion):
        return pred.alen<7 and max(abs(self.poss[pred][0]-self.poss[pray][0]),abs(self.poss[pred][1]-self.poss[pray][1]))<pred.alen-pray.alen
    def subsume(self,pred:minion.Minion,pray:minion.Minion):
        pred.take_mass(pray.mass)
        self.kill(pray,distribute=False)
    def kill(self,mi,distribute=True):
        self.mis.remove(mi)
        self.poss.__delitem__(mi)
        if distribute:
            for _ in range(mi.mass):
                i=random.randint(0,self.xsize-1)
                j=random.randint(0,self.ysize-1)
                self.mins[(i,j)]+=1
        del mi
    def kill_elderly(self):
        for mi in self.mis:
            if not mi.increase_age():
                self.kill(mi)
    def control_all(self):
        for mi in self.mis:
            mi.brain.control(mi,self.snapshot)
    def eat(self,mi:minion.Minion):
        x,y=self.poss[mi]
        sum=0
        for i in range(-mi.alen,mi.alen+1):
            for j in range(-mi.alen,mi.alen+1):
                if self.in_range((x+i,y+j)):
                    sum+=self.mins[(x+i,y+j)]
                    self.mins[(x+i,y+j)]=0
        mi.take_mass(sum)
    def try_move(self,mi:minion.Minion):
        x,y=self.poss[mi]
        a,b=0,0
        if mi.move_dist==-1:
            for girl in self.mis:
                if self.sexy(mi,girl):
                    girl:minion.Minion
                    if girl.count(mi):
                        self.childbirth(girl,mi)
            return
        move_direc=mi.move_direc
        move_dist = int(1.5 * mi.move_dist * sqrt(mi.mass))
        if move_direc==0:
            a+=min(move_dist,self.xsize-1-x)
        elif move_direc==1:
            b+=min(move_dist,self.ysize-1-y)
        elif move_direc==2:
            a-=min(move_dist,x)
        elif move_direc==3:
            b-=min(move_dist,y)
        self.poss[mi]=(x+a,y+b)


    def move_all(self):
        i=0
        while i<len(self.mis): #self.mis is dynamic. If a new child is added, it is added at the tail.
            mi=self.mis[i]
            self.try_move(mi)
            self.eat(mi)
            died=False
            for prey in self.mis:
                if prey!=mi and self.subsumable(mi,prey):
                    if self.mis.index(prey)<i:
                        i-=1
                    self.subsume(mi,prey)
            for pred in self.mis:
                if pred!=mi and self.subsumable(pred,mi):
                    died=True
                    self.subsume(pred,mi)
                    break
            if not died:
                i+=1
    def render(self):
        for i in range(self.xsize):
            for j in range(self.ysize):
                self.snapshot[i][j]=None
        for mi in self.mis:
            x,y=self.poss[mi]
            for i in range(-mi.alen,mi.alen+1):
                for j in range(-mi.alen,mi.alen+1):
                    if self.in_range((x+i,y+j)):
                        self.snapshot[x+i][y+j]=mi.c
            if mi.move_dist==-1:
                if mi.move_direc==0:
                    a,b=1,0
                elif mi.move_direc==1:
                    a,b=0,1
                elif mi.move_direc==2:
                    a,b=-1,0
                elif mi.move_direc==3:
                    a,b=0,-1
                for i in range(mi.alen):
                    if self.in_range((x+a*(mi.alen+i),y+b*(mi.alen+j))):
                        self.snapshot[x+a*(mi.alen+i)][y+b*(mi.alen+j)]=mi.c
        for i in range(self.xsize):
            for j in range(self.ysize):
                if self.snapshot[i][j]==None:
                    if self.mins[(i,j)]==0:
                        self.snapshot[i][j]=(0,0,0)
                    else:
                        self.snapshot[i][j]=(255,255,255)
    def one_step(self):
        self.control_all()
        self.move_all()
        self.moment+=1
        self.kill_elderly()
        self.render()