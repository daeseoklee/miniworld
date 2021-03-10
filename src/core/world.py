#import minion
import random
from math import *
import numpy as np
import copy



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
                 no_hunt=False,halluc=False,alen_list=None,pos_list=None,freeze_list=None,record_pedigree=False):
        self.snapshot=np.array([[0 for _ in range(ysize)] for _ in range(xsize)])
        self.xsize=xsize #int
        self.ysize=ysize #int
        self.moment=0
        self.new_id=0
        self.mins=[[0 for _ in range(ysize)] for _ in range(xsize)]
        nut=total_mass
        for i,mi in enumerate(mis):
            if alen_list!=None:
                mi.take_mass((1+2*alen_list[i])**2-mi.mass)
            nut-=mi.mass
        for _ in range(nut):
            i=random.randint(0,xsize-1)
            j=random.randint(0,ysize-1)
            self.mins[i][j]+=1
        self.mis=[] #minion.Minion[]
        self.pedigree=[] #(int,int)[]
        self.pos_dict=dict() #dict[minion.Minion](int,int)
        self.occupy_map=[[[] for _ in range(self.ysize)] for _ in range(self.xsize)] #minion.Minion[][][]
        for i,mi in enumerate(mis):
            if pos_list==None:
                xpos=random.randint(0,xsize-1)
                ypos=random.randint(0,ysize-1)
            else:
                xpos,ypos=pos_list[i]
            self.register(mi,(xpos,ypos))
            if freeze_list!=None and freeze_list[i]:
                mi.freeze()

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

    def register(self,mi,pos):
        #World*minion.Minion*(int,int)->()
        #use: register 'mi' in the world and locate it at 'pos'
        #called in: World.__init__(), World.childbirth()
        self.mis.append(mi) #don't change this order, because it's relevant in some tests
        mi.id=self.new_id
        self.new_id+=1
        self.pos_dict[mi]=pos
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
        x,y=self.pos_dict[mi]
        for i in range(-mi.alen,1+mi.alen):
            for j in range(-mi.alen,1+mi.alen):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)
    def skewbodyposs(self,mi):
        #World*minion.Minion->generator[(int,int)]
        #use: Similar to bodyposs, but only yield a proper subset which is sufficient for detecting all the other minions whose body is contained in that of 'mi'; more efficient for those cases
        #called in: World.try_hunt()
        x,y=self.pos_dict[mi]
        for i in range(-mi.alen,1+mi.alen,3):
            for j in range(-mi.alen,1+mi.alen,3):
                yield ((x+i)%self.xsize,(y+j)%self.ysize)

    def genipos(self,mi):
        #World*minion.Minion->(int,int)
        #use: determine the position of the reproductive organ
        #called in: World.stretch()
        x,y=self.pos_dict[mi]
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
        x,y=self.pos_dict[mi]
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
        print("alen:",mi.alen,"mass:",mi.mass,"amount:",amount,"pos:",self.pos_dict[mi])
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
            pos=self.pos_dict[girl]
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
            self.register(child,(x,y))
            if self.record_pedigree:
                self.pedigree.append((girl.id,boy.id,child.id))

    def huntable(self,pred,pray):
        #World*minion.Minion*minion.Minion->bool
        #use: determine whether 'pred' can hunt 'pray'
        #called in: World.try_hunt()
        return (not self.no_hunt) and max(dist(self.pos_dict[pred][0],self.pos_dict[pray][0],self.xsize),dist(self.pos_dict[pred][1],self.pos_dict[pray][1],self.ysize))<=pred.alen-pray.alen and pred!=pray
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
        self.pos_dict.__delitem__(mi)
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
                inp=mi.get_input(self.snapshot,self.pos_dict[mi])
                mi.brain.control(mi,inp)
        else:
            inps=[mi.get_input(self.snapshot,self.pos_dict[mi]) for mi in self.mis]
            random.shuffle(inps)
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
        x,y=self.pos_dict[mi]
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
        self.pos_dict[mi]=((x+a)%self.xsize,(y+b)%self.ysize)
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
                    self.snapshot[i][j]=256**3-1
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
    #measurements
    def get_pedigree(self):
        return self.pedigree
    def sample_dna(self,n):
        if n>=len(self.mis):
            subpopulation=random.sample(self.mis,n)
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
            x,y=self.pos_dict[mi]
            yield mi.color,(k*(x-mi.alen),k*(y-mi.alen),k*(1+2*mi.alen),k*(1+2*mi.alen))
    def geni_rect(self,k):
        #called in: play.py
        for mi in self.mis:
            if mi.action==1:
                x,y=self.pos_dict[mi]
                if mi.move_direc==0:
                    yield mi.color,(k*x,k*y,k*(1+2*mi.alen),k)
                elif mi.move_direc==1:
                    yield mi.color,(k*x,k*y,k,k*(1+2*mi.alen))
                elif mi.move_direc==2:
                    yield mi.color,(k*x+k-1,k*y,-k*(1+2*mi.alen),k)
                elif mi.move_direc==3:
                    yield mi.color,(k*x,k*y+k-1,k,-k*(1+2*mi.alen))



