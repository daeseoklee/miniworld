# distutils: language=c++
# cython: profile=True

import numpy as np
from libc.math cimport exp,sqrt,erf,log
from fast_random cimport uniform,randint,randbool,gaussian
from fast_random cimport multiple_uniform,multiple_gaussian
from fast_random cimport randint_with_probs,seed
from random import sample
from cython.view cimport array as cvarray

seed()


cdef:
    int num_out=12
    int vision_range=5
    int vision_resolution=9
    int vision_len=vision_range*vision_resolution
    int num_feature=5
    int rough_idim=(vision_range)**2*num_feature
    int rough_wdim=rough_idim*num_out
    int clear_idim=3*vision_len**2
    int linear_wdim=clear_idim*num_out
    double merge_thres=0.04
    double mut_per_diff=3.0 #minimum number of mutations required for a species division, if every change is monotonic
    double mut_rate=merge_thres/mut_per_diff
cdef:
    int xsize=200
    int ysize=150
    int min_maxsize=3
    int max_maxsize=8
    double min_uptake=0.1
    double max_uptake=0.7
    int min_maxage=1500
    int max_maxage=1500
cdef:
    double avg_consum_rate=0.01
    double consum_exp=0.7

def set_consum(a,e):
    global avg_consum_rate
    global consum_exp
    avg_consum_rate=a
    consum_exp=e

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

cdef class ColorTrait():
    def __cinit__(self,int r,int g,int b):
        self.r=r
        self.g=g
        self.b=b
    def get_c(self):
        return (self.r,self.g,self.b)
    def opposite(self):
        return ColorTrait(255-self.r,255-self.g,255-self.b)
    cdef void mutate(self,double rate):
        cdef double dr=256*rate*gaussian()
        cdef int intdr=<int>dr
        cdef double dg=256*rate*gaussian()
        cdef int intdg=<int>dg
        cdef double db=256*rate*gaussian()
        cdef int intdb=<int>db
        if dr>0:
            if dr-<double>intdr<0.5:
                self.r+=intdr
            else:
                self.r+=intdr+1
        else:
            if <double>intdr-dr<0.5:
                self.r+=intdr
            else:
                self.r+=intdr-1
        if self.r>255:
            self.r=255
        elif self.r<0:
            self.r=0

        if dg>0:
            if dg-<double>intdg<0.5:
                self.g+=intdg
            else:
                self.g+=intdg+1
        else:
            if <double>intdg-dg<0.5:
                self.g+=intdg
            else:
                self.g+=intdg-1
        if self.g>255:
            self.g=255
        elif self.g<0:
            self.g=0

        if db>0:
            if db-<double>intdb<0.5:
                self.b+=intdb
            else:
                self.b+=intdb+1
        else:
            if <double>intdb-db<0.5:
                self.b+=intdb
            else:
                self.b+=intdb-1
        if self.b>255:
            self.b=255
        elif self.b<0:
            self.b=0
    cdef ColorTrait mixed(ColorTrait t1,ColorTrait t2):
        cdef int r,g,b
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
    def get_normalized_difference(t1,t2):
        return ColorTrait.normalized_difference(t1,t2)
    cdef double normalized_difference(ColorTrait t1,ColorTrait t2):
        cdef double difference, normalizer
        difference=abs(t1.r-t2.r)+abs(t1.g-t2.g)+abs(t1.b-t2.b)
        normalizer=1+max(t1.r,t2.r)+max(t1.g,t2.g)+max(t1.b,t2.b)
        return difference/normalizer


cdef class PositiveTrait():
    def __cinit__(self,double a):
        self.eps=0.001
        if a>self.eps:
            self.a=a
        else:
            self.a=self.eps
    def get_a(self):
        return self.a
    cdef void mutate(self,double rate):
        self.a+=rate*self.a*gaussian()
        if self.a<self.eps:
            self.a=self.eps
    cdef PositiveTrait mixed(PositiveTrait t1, PositiveTrait t2):
        cdef double a
        if randbool():
            a=t1.a
        else:
            a=t2.a
        return PositiveTrait(a)
    cdef double normalized_difference(PositiveTrait t1,PositiveTrait t2):
        cdef double difference,normalizer
        difference=abs(t1.a-t2.a)
        normalizer=max(t1.a,t2.a)
        return difference/normalizer


cdef class UnitTrait():
    def __cinit__(self,double a):
        self.eps=0.001
        if a>self.eps:
            self.a=a
        elif a>1:
            self.a=1
        else:
            self.a=self.eps
    def get_a(self):
        return self.a
    cdef void mutate(self,double rate):
        self.a+=rate*gaussian()
        if self.a<self.eps:
            self.a=self.eps
        elif self.a>1:
            self.a=1
    cdef UnitTrait mixed(UnitTrait t1,UnitTrait t2):
        cdef double a
        if randbool():
            a=t1.a
        else:
            a=t2.a
        return UnitTrait(a)
    cdef double normalized_difference(UnitTrait t1,UnitTrait t2):
        cdef double difference,normalizer
        difference=abs(t1.a-t2.a)
        normalizer=max(t1.a,t2.a)
        return difference/normalizer
cdef class FloatListTrait():
    def __cinit__(self,double[:] l,size_t[:] group_sizes):
        self.l=l
        self.group_sizes=group_sizes
    def yield_l(self):
        for i in range(self.l.shape[0]):
            yield self.l[i]
    def get_group_sizes(self):
        n=self.group_sizes.shape[0]
        l=[self.group_sizes[i] for i in range(n)]
        return l
    cdef void mutate(self,double rate):
        #implemented in low level for efficiency
        cdef size_t i,n,m,remainder
        cdef double[:] dl
        cdef double w,x1,x2,r1,r2
        n=self.l.shape[0]
        m=n/2
        remainder=n-2*m
        for i in range(m):
            w = 2.0
            while (w >= 1.0):
                x1 = 2.0 * uniform() - 1.0
                x2 = 2.0 * uniform() - 1.0
                w = x1 * x1 + x2 * x2
            w = sqrt((-2.0 * log(w)) / w)
            self.l[2*i]+=rate*x1*w
            self.l[2*i+1]+=rate*x2*2
        if remainder==1:
            w = 2.0
            while (w >= 1.0):
                x1 = 2.0 * uniform() - 1.0
                x2 = 2.0 * uniform() - 1.0
                w = x1 * x1 + x2 * x2
            w = sqrt((-2.0 * log(w)) / w)
            self.l[n-1]+=rate*x1*w

    cdef FloatListTrait mixed(FloatListTrait t1,FloatListTrait t2):
        cdef double[:] l
        cdef size_t n,i,m,j,ini_index,fin_index
        n=t1.l.shape[0]
        m=t1.group_sizes.shape[0]
        ini_index=0
        l=np.empty(n,dtype=np.float64,order='C')
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

    cdef double normalized_difference(FloatListTrait t1,FloatListTrait t2):
        cdef double difference,normalizer
        cdef size_t n,i
        difference=0
        n=t1.l.shape[0]
        for i in range(n):
            difference+=abs(t1.l[i]-t2.l[i])
        normalizer=n
        return difference/normalizer


cdef void take_softmax(double[:] l,double c):
    #np.ndarray(dim 1)*float->np.ndarray(dim 1)
    #called in: LinearBrain.control()
    cdef double sum
    cdef size_t n,i
    sum=0
    n=l.shape[0]
    for i in range(n):
        l[i]=exp(c*l[i])
        sum+=l[i]
    for i in range(n):
        l[i]/=sum

cdef double sigmoid(double x,double c):
    return 1/(1+exp(-c*x))

#---------------------------------------------------------------

cpdef LinearDNA randLinearDNA_with((int,int,int) c,double mst_a,double utt_a,double mat_a):
    cdef int r,g,b
    cdef double maxsize,uptake,maxage
    cdef double[:] weights
    cdef size_t[:] group_sizes
    cdef ColorTrait colorTrait
    cdef PositiveTrait maxsizeTrait
    cdef UnitTrait uptakeTrait
    cdef PositiveTrait maxageTrait
    cdef FloatListTrait weightsTrait
    cdef size_t i
    r=c[0]
    g=c[1]
    b=c[2]
    weights=multiple_uniform(linear_wdim)
    for i in range(<size_t>linear_wdim):
        weights[i]=1-2*weights[i]
    group_sizes=np.empty(num_out,dtype=np.uintp)
    for i in range(<size_t>num_out):
        group_sizes[i]=clear_idim
    colorTrait=ColorTrait(r,g,b)
    maxsizeTrait=PositiveTrait(mst_a)
    uptakeTrait=UnitTrait(utt_a)
    maxageTrait=PositiveTrait(mat_a)
    weightsTrait=FloatListTrait(weights,group_sizes)
    return LinearDNA(colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)
cpdef LinearDNA randLinearDNA():
    cdef int r,g,b
    cdef double maxsize,uptake,maxage
    cdef double[:] weights
    cdef size_t[:] group_sizes
    cdef ColorTrait colorTrait
    cdef PositiveTrait maxsizeTrait
    cdef UnitTrait uptakeTrait
    cdef PositiveTrait maxageTrait
    cdef FloatListTrait weightsTrait
    cdef size_t i
    r=randint(0,255)
    g=randint(0,255)
    b=randint(0,255)
    maxsize=randint(min_maxsize,max_maxsize)
    uptake=min_uptake+(max_uptake-min_uptake)*uniform()
    maxage=randint(min_maxage,max_maxage)
    weights=multiple_uniform(linear_wdim)
    for i in range(<size_t>linear_wdim):
        weights[i]=1-2*weights[i]
    group_sizes=np.empty(num_out,dtype=np.uintp)
    for i in range(<size_t>num_out):
        group_sizes[i]=clear_idim
    colorTrait=ColorTrait(r,g,b)
    maxsizeTrait=PositiveTrait(maxsize)
    uptakeTrait=UnitTrait(uptake)
    maxageTrait=PositiveTrait(maxage)
    weightsTrait=FloatListTrait(weights,group_sizes)
    return LinearDNA(colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)

def write_list_of_dnas_file(l,filename):
    f=open(filename,'w')
    f.write(str(len(l))+"\n")
    f.write(str(l[0][0].get_wdim())+"\n")
    group_sizes=l[0][0].get_group_sizes()
    f.write(str(len(group_sizes))+"\n")
    for a in group_sizes:
        f.write(str(a)+"\n")
    j=0
    for dnas in l:
        j+=1
        f.write("#dnas"+str(j)+"\n")
        i=0
        for dna in dnas:
            i+=1
            f.write("#dna"+str(i)+"\n")
            ct,mst,utt,mat,wt=dna.get_traits()
            f.write("<color>\n")
            r,g,b=ct.get_c()
            f.write(str(r)+"\n")
            f.write(str(g)+"\n")
            f.write(str(b)+"\n")
            f.write("<maxsize>\n")
            a=mst.get_a()
            f.write(str(a)+"\n")
            f.write("<uptake>\n")
            a=utt.get_a()
            f.write(str(a)+"\n")
            f.write("<maxage>\n")
            a=mat.get_a()
            f.write(str(a)+"\n")
            f.write("<weights>\n")
            for a in wt.yield_l():
                f.write(str(a)+"\n")
            f.write("#end_dna\n")
        f.write("#end_dnas\n")
    f.write("#end")
    f.close()
def read_list_of_dnas_file(filename):
    f=open(filename,'r')
    lines=f.readlines()
    f.close()
    num=int(lines[0])
    wdim=int(lines[1])
    group_sizes_len=int(lines[2])
    group_sizes=np.empty(group_sizes_len,dtype=np.uintp)
    for i in range(group_sizes_len):
        group_sizes[i]=int(lines[3+i])
    at=3+group_sizes_len
    dnas_list=[]
    for _ in range(num):
        dnas=[]
        while True: #parsing one dna per loop
            while len(lines[at])<4 or lines[at][:4]!="#dna":
                at+=1
            while len(lines[at])<1 or lines[at][0]!="<":
                at+=1
            assert lines[at]=="<color>\n"
            r=int(lines[at+1])
            g=int(lines[at+2])
            b=int(lines[at+3])
            ct=ColorTrait(r,g,b)
            at+=3
            while len(lines[at])<1 or lines[at][0]!="<":
                at+=1
            assert lines[at]=="<maxsize>\n"
            a=float(lines[at+1])
            mst=PositiveTrait(a)
            at+=1
            while len(lines[at])<1 or lines[at][0]!="<":
                at+=1
            assert lines[at]=="<uptake>\n"
            a=float(lines[at+1])
            utt=UnitTrait(a)
            at+=1
            while len(lines[at])<1 or lines[at][0]!="<":
                at+=1
            assert lines[at]=="<maxage>\n"
            a=float(lines[at+1])
            mat=PositiveTrait(a)
            at+=1
            while len(lines[at])<1 or lines[at][0]!="<":
                at+=1
            assert lines[at]=="<weights>\n"
            at+=1
            l=cvarray(shape=(wdim,), itemsize=sizeof(double), format="d")
            i=0
            while lines[at]!="#end_dna\n":
                l[i]=float(lines[at])
                at+=1
                i+=1
            at+=1
            wt=FloatListTrait(l,group_sizes)
            dna=LinearDNA(ct,mst,utt,mat,wt)
            dnas.append(dna)
            if lines[at]=="#end_dnas\n":
                break
        dnas_list.append(dnas)
    return dnas_list








cdef class LinearDNA():
    def __cinit__(self,ColorTrait colorTrait,PositiveTrait maxsizeTrait,UnitTrait uptakeTrait,PositiveTrait maxageTrait,FloatListTrait weightsTrait):
        self.colorTrait=colorTrait
        self.maxsizeTrait=maxsizeTrait
        self.uptakeTrait=uptakeTrait
        self.maxageTrait=maxageTrait
        self.weightsTrait=weightsTrait
    def get_wdim(self):
        return self.weightsTrait.l.shape[0]
    def get_group_sizes(self):
        return self.weightsTrait.get_group_sizes()
    def get_traits(self):
        return (self.colorTrait,self.maxsizeTrait,self.uptakeTrait,self.maxageTrait,self.weightsTrait)
    def get_colorTrait(self):
        return self.colorTrait
    cdef bint mergeable(self,LinearDNA dna):
        if self.colorTrait.r!=dna.colorTrait.r:
            return False
        if self.colorTrait.g!=dna.colorTrait.g:
            return False
        if self.colorTrait.b!=dna.colorTrait.b:
            return False
        cdef double d
        d=FloatListTrait.normalized_difference(self.weightsTrait,dna.weightsTrait)
        return uniform()<1-erf(d/merge_thres)

    cdef LinearDNA merge(self,LinearDNA dna):
        cdef ColorTrait colorTrait
        colorTrait=self.colorTrait
        cdef PositiveTrait maxsizeTrait
        maxsizeTrait=self.maxsizeTrait
        cdef UnitTrait uptakeTrait
        uptakeTrait=self.uptakeTrait
        cdef PositiveTrait maxageTrait
        maxageTrait=self.maxageTrait
        cdef FloatListTrait weightsTrait
        weightsTrait=FloatListTrait.mixed(self.weightsTrait,dna.weightsTrait)
        weightsTrait.mutate(mut_rate)
        return LinearDNA(colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)
    cdef ((int,int,int),double,double,double) translate(self):
        return (self.colorTrait.r,self.colorTrait.g,self.colorTrait.b),\
               self.maxsizeTrait.a,self.uptakeTrait.a,self.maxageTrait.a
    cdef double[:] translate_weights(self):
        return self.weightsTrait.l

#--------------------------------------------------------------

cdef class Brain():
    cdef void control(self,Minion mi,double[:] inp):
        pass


cdef double[:] apply_linear(double[:] v,double[:] A):
    cdef size_t n,m,i,j
    cdef double[:] result
    cdef double sum
    n=v.shape[0]
    m=A.shape[0]
    m/=n
    result = cvarray(shape=(m,), itemsize=sizeof(double), format="d")
    for j in range(m):
        sum=0
        for i in range(n):
            sum+=v[i]*A[j*n+i]
        result[j]=sum
    return result



cdef class LinearBrain(Brain):
    def __cinit__(self,double[:] weights,bint init=False):
        if init:
            weights=1-2*np.random.rand(linear_wdim)
        self.weights=weights
    cdef void control(self,Minion mi,double[:] inp):
        if mi.frozen:
            return
        cdef double[:] out,direc_probs,action_probs
        cdef double c=1/sqrt((<double>clear_idim)*5/16)
        out=apply_linear(inp,self.weights)
        direc_probs=out[:4]
        take_softmax(direc_probs,c)
        action_probs=out[8:12]
        take_softmax(action_probs,c)
        mi.move_direc=randint_with_probs(0,3,direc_probs)
        mi.move_dist=sigmoid(out[4+<size_t>mi.move_direc],c)
        mi.action=randint_with_probs(0,3,action_probs)





#--------------------------------------------------------------

cpdef Minion construct_minion(LinearDNA dna,int alen,(int,int) pos,bint do_freeze):
    cdef Minion mi
    mi=Minion(dna)
    mi.take_mass(<double>((1+2*alen)**2-9))
    if pos==(-1,-1):
        pos[0]=randint(0,xsize-1)
        pos[1]=randint(0,ysize-1)
    mi.pos[0]=<size_t>pos[0]
    mi.pos[1]=<size_t>pos[1]
    if do_freeze:
        mi.freeze()
    return mi



cdef class Minion():
    def __cinit__(self,LinearDNA dna):
        self.dna=dna
        self.color,self.maxsize,self.uptake,self.maxage=self.dna.translate()
        self.brain=LinearBrain(self.dna.translate_weights(),init=False)
        self.idim=clear_idim
        self.id=-1
        self.pos=(0,0)
        self.age=0
        self.dead=False
        self.alen=1
        self.mass=9.0
        self.action=0
        self.move_direc=0
        self.move_dist=0.0
        self.cum_dist=0
        self.frozen=False
        self.node=None


        self.energy=5.0*self.mass
        self.adjust_energy_from_mass()
        self.avg_consum_rate=avg_consum_rate
        self.basal_metabolic_rate=self.avg_consum_rate/4 #this*self.mass per moment
        self.move_consum_rate=self.avg_consum_rate/40
        self.stretch_consum_rate=self.avg_consum_rate/4
        self.sex_consum_rate=self.avg_consum_rate/4
        self.birth_consum_rate=0.3
    def get_dna(self):
        return self.dna
    def get_color(self):
        return self.color
    def get_mass(self):
        return self.mass
    def get_alen(self):
        return self.alen
    def get_pos(self):
        return self.pos
    def get_action(self):
        return self.action
    def get_move_direc(self):
        return self.move_direc
    def get_frozen(self):
        return self.frozen
    def get_cum_dist(self):
        return self.cum_dist
    def freeze(self):
        self.frozen=True
    def defreeze(self):
        self.frozen=False
    def set_alen(self,n):
        self.alen=n
    def set_pos(self,a,b):
        self.pos=(a,b)

    cdef bint increase_age(self):
        self.age+=1
        if <double>self.age>self.maxage:
            return False
        return True

    cdef double energy_with_constant(self,double const):
        return const*self.mass*(self.mass/9)**consum_exp

    cdef void take_energy(self,double amount):
        if amount<=self.max_energy-self.energy:
            self.energy+=amount
        else:
            self.energy=self.max_energy
    cdef bint loss_energy(self,double amount):
        self.energy-=amount
        if self.energy<=0:
            return False
        return True

    cdef void adjust_energy_from_mass(self):
        self.max_energy=10*self.mass
        if self.energy>self.max_energy:
            self.energy=self.max_energy
    cdef void take_mass(self,double amount):
        #print("amount",amount)
        self.mass+=amount
        self.adjust_energy_from_mass()
        self.alen=<int>((sqrt(self.mass)-1)/2)
    cdef void loss_mass(self,double amount):
        self.mass-=amount
        self.adjust_energy_from_mass()
        self.alen=<int>((sqrt(self.mass)-1)/2)
    cdef bint mergeable(self,Minion mi):
        return self.alen>=2 and self.dna.mergeable(mi.dna)
    cdef Minion get_child(self,Minion mi):
        return Minion(self.dna.merge(mi.dna))

    cdef double[:] get_input(self,int[:,:] snapshot):
        cdef int xsize,ysize
        cdef double x_min,x_max,y_min,y_max
        cdef size_t i,j, x,y,vlen
        cdef int n
        cdef double r,g,b
        cdef double[:] inp
        xsize=snapshot.shape[0]
        ysize=snapshot.shape[1]
        x_min = self.pos[0] - (1 + 2 * self.alen) * ((vision_range - 1) / 2) - self.alen
        x_max = self.pos[0] + (1 + 2 * self.alen) * ((vision_range - 1) / 2) + self.alen
        y_min = self.pos[1] - (1 + 2 * self.alen) * ((vision_range - 1) / 2) - self.alen
        y_max = self.pos[1] + (1 + 2 * self.alen) * ((vision_range - 1) / 2) + self.alen
        inp=cvarray(shape=(clear_idim,), itemsize=sizeof(double), format="d")
        vlen=<size_t>vision_len
        for i in range(vlen):
            x=<size_t>((<int>(x_min+i*(x_max-x_min)/(vision_len-1)))%xsize)
            for j in range(vlen):
                y=<size_t>((<int>(y_min+j*(y_max-y_min)/(vision_len-1)))%ysize)
                n=snapshot[x,y]
                r=(<double>(n/256**2))/256
                g=(<double>((n%256**2)/256))/256
                b=(<double>(n%256))/256
                inp[i*vision_len+j]=r
                inp[vision_len**2+i*vision_len+j]=g
                inp[2*vision_len**2+i*vision_len+j]=b
        #print(inp[7*vision_len+7],inp[vision_len**2+7*vision_len+7],inp[2*vision_len**2+7*vision_len+7])
        return inp


#--------------------------------------------------------------
cdef class MinionDLLNode():
    def __cinit__(self,MinionDLL dll,Minion mi,MinionDLLNode prev,MinionDLLNode next):
        self.mi=mi
        self.prev=prev
        self.next=next
        if next==None:
            self.is_tail=True
        else:
            self.is_tail=False
        self.dll=dll
    def get_mi(self):
        return self.mi
def dll_from_list(mi_list):
    dll=MinionDLL()
    for mi in mi_list:
        dll.push_py(mi)
    return dll
cdef class MinionDLL():
    def __cinit__(self):
        self.len=0
        self.head=None
        self.tail=None
        self.current=None
    def is_empty(self):
        return self.len==0
    def push_py(self,mi):
        self.push(mi)
    def remove_py(self,mi):
        self.remove_by_search(mi)
    def contains_py(self,mi):
        return self.contains(mi)
    def get_nth_py(self,n):
        if self.len<=n:
            raise Exception("MinionDLL out of range!!")
        node=self.head
        if n==0:
            return node.get_mi()
        for _ in range(n-1):
            node=self.next(node)
        return node.get_mi()
    def __iter__(self):
        self.current=self.head
        return self
    def __next__(self):
        if self.current==None:
            raise StopIteration
        mi=self.current.get_mi()
        self.current=self.next(self.current)
        return mi
    cdef MinionDLLNode prev(self,MinionDLLNode node):
        return node.prev
    cdef MinionDLLNode next(self,MinionDLLNode node):
        return node.next
    cdef MinionDLLNode push(self,Minion mi):
        cdef MinionDLLNode node,current_tail
        if self.len==0:
            node=MinionDLLNode(self,mi,None,None)
            self.head=node
            self.tail=node
            self.len=1
            return
        current_tail=self.tail
        node=MinionDLLNode(self,mi,current_tail,None)
        current_tail.next=node
        current_tail.is_tail=False
        self.tail=node
        self.len+=1
        return node
    cdef MinionDLLNode remove_and_get_next(self,MinionDLLNode node):
        cdef MinionDLLNode prev,next
        if self.head==node:
            if self.len==1:
                self.head=None
                self.tail=None
                self.len=0
                return None
            next=self.next(node)
            next.prev=None
            self.head=next
            self.len-=1
            return next
        if self.tail==node:
            prev=self.prev(node)
            prev.next=None
            prev.is_tail=True
            self.tail=prev
            self.len-=1
            return None
        prev=self.prev(node)
        next=self.next(node)
        prev.next=next
        next.prev=prev
        self.len-=1
        return next
    cdef void remove_by_link(self,Minion mi):
        self.remove_and_get_next(mi.node)
    cdef void remove_by_search(self,Minion mi):
        if self.len==0:
            return
        cdef MinionDLLNode node
        node=self.head
        while (node.mi is not mi) and (not node.is_tail):
            node=self.next(node)
        if node.mi is mi:
            self.remove_and_get_next(node)
    cdef bint contains(self,Minion mi):
        if self.len==0:
            return False
        cdef MinionDLLNode node
        node=self.head
        while not node.is_tail:
            if node.mi==mi:
                return True
            node=self.next(node)
        if node.mi==mi:
            return True
        return False


#previous world.py--------------------------------------------------------------



cdef int dist(int x,int y,int T):
    if x>y:
        return min(x-y,T+y-x)
    return min(y-x,T+x-y)

cdef class World():
    def __cinit__(self,int xsize,int ysize,double total_mass,MinionDLL mis,\
                  bint no_age=False,bint no_birth=False,bint no_eat=False,\
                  bint no_energy=False,bint no_excrete=False,bint no_hunt=False,\
                  bint halluc=False,bint record_pedigree=False):
        seed()
        self.snapshot=np.array([[-1 for _ in range(ysize)] for _ in range(xsize)],dtype=np.int32,order='C')
        self.xsize=xsize
        self.ysize=ysize
        self.moment=0
        self.new_id=0
        self.mins=np.array([[0. for _ in range(ysize)] for _ in range(xsize)],dtype=np.float64,order='C')
        cdef double initial_total_min=total_mass
        cdef size_t x,y
        for mi in mis:
            initial_total_min-=mi.get_mass()
        for _ in range(<size_t>initial_total_min):
            x=<size_t>randint(0,xsize-1)
            y=<size_t>randint(0,ysize-1)
            self.mins[x,y]+=1

        self.mis=MinionDLL()
        #how two assign an empty vector?
        self.occupy_map=np.array([[MinionDLL() for _ in range(ysize)] for _ in range(xsize)])
        for mi in mis:
            self._register(mi)
        self.no_age=no_age
        self.no_birth=no_birth
        self.no_eat=no_eat
        self.no_energy=no_energy
        self.no_excrete=no_excrete
        self.no_hunt=no_hunt
        self.halluc=halluc
        self.messiness=0
        self.hidden_mass=0.
        self.record_pedigree=record_pedigree
        self.pedigree=[] #python object


    def _register(self,mi):
        self.register(mi)

    cdef void register(self,Minion mi):
        cdef MinionDLLNode node
        node=self.mis.push(mi)
        mi.node=node
        mi.id=self.new_id
        self.new_id+=1
        cdef int xpos,ypos
        cdef size_t i,j,x,y
        cdef MinionDLL dll
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        for i in range(1+2*<size_t>mi.alen):
            for j in range(1+2*<size_t>mi.alen):
                x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                dll=self.occupy_map[x,y]
                dll.push(mi)
    cdef void unregister_deads(self):
        cdef MinionDLLNode node
        cdef Minion mi
        if self.mis.len==0:
            return
        node=self.mis.head
        while not node.is_tail:
            if node.mi.dead:
                node=self.mis.remove_and_get_next(node)
            else:
                node=self.mis.next(node)
        if node.mi.dead:
            self.mis.remove_and_get_next(node)


    cdef void take_mass(self,Minion mi,double amount):
        cdef int alen_before,xpos,ypos
        cdef size_t i,j,x,y
        cdef MinionDLL dll
        alen_before=mi.alen
        mi.take_mass(amount)
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        if mi.alen>alen_before:
            for i in range(1+2*mi.alen):
                x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                for j in range(mi.alen-alen_before):
                    y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                    dll=self.occupy_map[x,y]
                    dll.push(mi)
                for j in range(1+mi.alen+alen_before,1+2*mi.alen):
                    y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                    dll=self.occupy_map[x,y]
                    dll.push(mi)
            for j in range(mi.alen-alen_before,1+mi.alen+alen_before):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                for i in range(mi.alen-alen_before):
                    x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                    dll=self.occupy_map[x,y]
                    dll.push(mi)
                for i in range(1+mi.alen+alen_before,1+2*mi.alen):
                    x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                    dll=self.occupy_map[x,y]
                    dll.push(mi)
    cdef void loss_mass(self,Minion mi,double amount):
        cdef int alen_before,xpos,ypos
        cdef size_t i,j,x,y
        cdef MinionDLL dll
        alen_before=mi.alen
        mi.loss_mass(amount)
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        if alen_before>mi.alen:
            for i in range(1+2*alen_before):
                x=<size_t>((xpos-alen_before+<int>i)%self.xsize)
                for j in range(alen_before-mi.alen):
                    y=<size_t>((ypos-alen_before+<int>j)%self.ysize)
                    dll=self.occupy_map[x,y]
                    dll.remove_by_search(mi)
                for j in range(1+alen_before+mi.alen,1+2*alen_before):
                    y=<size_t>((ypos-alen_before+<int>j)%self.ysize)
                    dll=self.occupy_map[x,y]
                    dll.remove_by_search(mi)
            for j in range(alen_before-mi.alen,1+alen_before+mi.alen):
                y=<size_t>((ypos-alen_before+<int>j)%self.ysize)
                for i in range(alen_before-mi.alen):
                    x=<size_t>((xpos-alen_before+<int>i)%self.xsize)
                    dll=self.occupy_map[x,y]
                    dll.remove_by_search(mi)
                for i in range(1+alen_before+mi.alen,1+2*alen_before):
                    x=<size_t>((xpos-alen_before+<int>i)%self.xsize)
                    dll=self.occupy_map[x,y]
                    dll.remove_by_search(mi)
    cdef void take_energy(self,Minion mi,double amount):
        if self.no_energy:
            return
        mi.take_energy(amount)
    cdef bint loss_energy(self,Minion mi,double amount):
        if self.no_energy:
            return True
        cdef bint alive
        alive=mi.loss_energy(amount)
        if not alive:
            #self.death_by("lack of energy")
            self.kill(mi,corpse=True)
        return alive

    cdef void childbirth(self,Minion girl,Minion boy):
        if self.no_birth:
            return
        cdef Minion child
        cdef size_t x,y
        cdef int xpos,ypos
        cdef bint alive
        xpos=<int>girl.pos[0]
        ypos=<int>girl.pos[1]
        if girl.mergeable(boy):
            child=girl.get_child(boy)
            if girl.move_direc==0:
                x=<size_t>((xpos-girl.alen-3)%self.xsize)
                y=<size_t>ypos
            elif girl.move_direc==1:
                x=<size_t>xpos
                y=<size_t>((ypos-girl.alen-3)%self.ysize)
            elif girl.move_direc==2:
                x=<size_t>((xpos+girl.alen+3)%self.xsize)
                y=<size_t>ypos
            elif girl.move_direc==3:
                x=<size_t>xpos
                y=<size_t>((ypos+girl.alen+3)%self.ysize)
            child.pos=(x,y)
            self.loss_mass(girl,child.mass)
            alive=self.loss_energy(girl,girl.energy_with_constant(girl.birth_consum_rate))
            if alive and girl.alen<1:
                self.kill(girl,corpse=True)
            self.register(child)
            if self.record_pedigree:
                self.pedigree.push_back((girl.id,boy.id,child.id))
    cdef bint huntable(self,Minion pred,Minion pray):
        return (not self.no_hunt) and max(dist(pred.pos[0],pray.pos[0],self.xsize),dist(pred.pos[1],pray.pos[1],self.ysize))<=pred.alen-pray.alen and pred!=pray

    cdef void mk_corpse(self,Minion mi):
        cdef int area,xpos,ypos,k
        cdef double q,r
        cdef size_t i,j,x,y
        area=(1+2*mi.alen)**2
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        q=<double>((<int>mi.mass)/(<int>area))
        r=mi.mass%area
        for i in range(1+2*<size_t>mi.alen):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(1+2*<size_t>mi.alen):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                k=(<int>i)*(1+2*mi.alen)+<int>j
                if (<double>k)<=r-1:
                    self.mins[x,y]+=q+1
                elif r-1<(<double>k)<r:
                    self.mins[x,y]+=q+r-<double><int>r
                else:
                    self.mins[x,y]+=q

    cdef void kill(self,Minion mi,bint corpse):
        if mi.dead:
            return
        if corpse:
            self.mk_corpse(mi)
        cdef int xpos,ypos
        cdef size_t i,j,x,y
        cdef MinionDLL dll
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        for i in range(1+2*<size_t>mi.alen):
            for j in range(1+2*<size_t>mi.alen):
                x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                dll=self.occupy_map[x,y]
                dll.remove_by_search(mi)
        mi.dead=True

    cdef void kill_elderly(self):
        if self.no_age or self.mis.len==0:
            return
        cdef MinionDLLNode node,next_node
        cdef Minion mi
        node=self.mis.head
        while not node.is_tail:
            next_node=self.mis.next(node)
            mi=node.mi
            if not mi.increase_age():
                #self.death_by("age")
                self.kill(mi,corpse=True)
            node=next_node
        mi=node.mi
        if not mi.increase_age():
            #self.death_by("age")
            self.kill(mi,corpse=True)

    cdef void hunt(self,Minion pred,Minion pray):
        self.digest(pred,pray.mass)
        #self.death_by("hunt")
        self.kill(pray,corpse=False)


    cdef void move(self,Minion mi):
        cdef int move_direc,move_dist,a,b,xpos,ypos
        cdef size_t x,y,i,j
        cdef MinionDLL dll
        if mi.move_dist==0.0:
            return
        move_direc=mi.move_direc
        move_dist=<int>(2*mi.move_dist*mi.mass**0.5)
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        a=0
        b=0
        if move_direc==0:
            a+=move_dist
        elif move_direc==1:
            b+=move_dist
        elif move_direc==2:
            a-=move_dist
        elif move_direc==3:
            b-=move_dist
        for i in range(1+2*<size_t>mi.alen):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(1+2*<size_t>mi.alen):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                dll=self.occupy_map[x,y]
                dll.remove_by_search(mi)
        mi.pos[0]=((xpos+a)%self.xsize)
        mi.pos[1]=((ypos+b)%self.ysize)
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        for i in range(1+2*<size_t>mi.alen):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(1+2*<size_t>mi.alen):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                dll=self.occupy_map[x,y]
                dll.push(mi)
        self.loss_energy(mi,<double>move_dist*mi.energy_with_constant(mi.move_consum_rate))
        mi.cum_dist+=move_dist
    cdef void stretch(self,Minion mi):
        cdef int xpos,ypos
        cdef size_t x,y
        cdef bint alive
        cdef MinionDLL dll

        alive=self.loss_energy(mi,mi.energy_with_constant(mi.stretch_consum_rate))
        if not alive:
            return

        xpos=mi.pos[0]
        ypos=mi.pos[1]
        if mi.move_direc==0:
            x=<size_t>((xpos+2*mi.alen)%self.xsize)
            y=<size_t>ypos
        elif mi.move_direc==1:
            x=<size_t>xpos
            y=<size_t>((ypos+2*mi.alen)%self.ysize)
        elif mi.move_direc==2:
            x=<size_t>((xpos-2*mi.alen)%self.xsize)
            y=<size_t>ypos
        elif mi.move_direc==3:
            x=<size_t>xpos
            y=<size_t>((ypos-2*mi.alen)%self.ysize)

        #be careful for newly born child being reached!
        #(that doesn't happen as long as child's occupying area and
        #mother's occupying area are guaranteed to be disjoint
        dll=self.occupy_map[x,y]
        if dll.len==0:
            return
        cdef MinionDLLNode node,next_node
        cdef Minion other
        node=dll.head
        while not node.is_tail:
            next_node=dll.next(node)
            other=node.mi
            self.childbirth(other,mi)
            self.messiness+=1
            alive=self.loss_energy(mi,mi.energy_with_constant(mi.sex_consum_rate))
            if not alive:
                return
            node=next_node
        other=node.mi
        self.childbirth(other,mi)
        self.messiness+=1

    cdef void excrete(self,Minion mi,double amount):
        if self.no_excrete:
            self.hidden_mass+=amount
            return
        cdef int dist,xpos,ypos
        cdef size_t i,j,x,y
        cdef double area,q,r
        cdef int k
        dist=4*mi.alen
        if mi.move_direc==0:
            xpos=mi.pos[0]-dist
            ypos=mi.pos[1]
        elif mi.move_direc==1:
            xpos=mi.pos[0]
            ypos=mi.pos[1]-dist
        elif mi.move_direc==2:
            xpos=mi.pos[0]+dist
            ypos=mi.pos[1]
        elif mi.move_direc==3:
            xpos=mi.pos[0]
            ypos=mi.pos[1]+dist
        area=<double>(1+2*mi.alen)**2
        q=<double>((<int>amount)/(<int>area))
        r=amount%area
        for i in range(1+2*<size_t>mi.alen):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(1+2*<size_t>mi.alen):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                k=(<int>i)*(1+2*mi.alen)+<int>j
                if (<double>k)<=r-1:
                    self.mins[x,y]+=q+1
                elif r-1<(<double>k)<r:
                    self.mins[x,y]+=q+r-<double><int>r
                else:
                    self.mins[x,y]+=q
    cdef void digest(self,Minion mi,double amount):
        cdef double take,out
        take=min(mi.uptake*amount,(1+2*mi.maxsize)**2-mi.mass)
        out=amount-take
        self.take_mass(mi,take)
        self.excrete(mi,out)
        self.take_energy(mi,out)
    cdef void try_hunt(self,Minion mi):
        cdef int xpos,ypos
        cdef size_t i,j,x,y
        cdef Minion other
        cdef MinionDLL occupying,to_try
        cdef MinionDLLNode node
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        to_try=MinionDLL()
        for i in range(0,1+2*<size_t>mi.alen,3):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(0,1+2*<size_t>mi.alen,3):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                occupying=self.occupy_map[x,y]
                if occupying.len==0:
                    continue
                node=occupying.head
                while not node.is_tail:
                    other=node.mi
                    if self.huntable(mi,other) and (not to_try.contains(other)):
                        to_try.push(other)
                    node=occupying.next(node)
                other=node.mi
                if self.huntable(mi,other) and (not to_try.contains(other)):
                    to_try.push(other)
        if to_try.len==0:
            return
        node=to_try.head
        while not node.is_tail:
            other=node.mi
            self.hunt(mi,other)
            node=to_try.next(node)
        other=node.mi
        self.hunt(mi,other)
    cdef void eat(self,Minion mi):
        cdef double total
        cdef int xpos,ypos
        cdef size_t i,j,x,y
        if self.no_eat:
            return
        total=0
        xpos=mi.pos[0]
        ypos=mi.pos[1]
        for i in range(1+2*<size_t>mi.alen):
            x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
            for j in range(1+2*<size_t>mi.alen):
                y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                total+=self.mins[x,y]
                self.mins[x,y]=0
        self.digest(mi,total)
    cdef void act(self,Minion mi):
        if mi.frozen:
            return
        #print("mi.action",mi.action)
        if mi.action==0:
            self.move(mi)
            #print("mi.move_dist",mi.move_dist)
        elif mi.action==1:
            self.stretch(mi)
        elif mi.action==2:
            self.try_hunt(mi)
        elif mi.action==3:
            self.eat(mi)
    cdef void basal_metabolism(self,Minion mi):
        self.loss_energy(mi,mi.basal_metabolic_rate*mi.mass)



    cdef void control_all(self):
        if self.mis.len==0:
            return
        cdef MinionDLLNode node
        cdef Minion mi
        cdef double[:] inp
        cdef double[:,:] inps
        cdef size_t i
        if not self.halluc:
            node=self.mis.head
            while not node.is_tail:
                mi=node.mi
                inp=mi.get_input(self.snapshot)
                #self.snapshot_error_detector("right before control")
                #if node.prev==None:
                 #   print(list(inp))
                #print(inp[7*vision_len+7],inp[vision_len**2+7*vision_len+7],inp[2*vision_len**2+7*vision_len+7])
                mi.brain.control(mi,inp)
                node=self.mis.next(node)
            mi=node.mi
            inp=mi.get_input(self.snapshot)
            mi.brain.control(mi,inp)
            return
        v=np.empty((self.mis.len,self.mis.head.mi.idim),dtype=np.float64)
        inps=v
        node=self.mis.head
        i=0
        while not node.is_tail:
            mi=node.mi
            inps[i,:]=mi.get_input(self.snapshot)
            node=self.mis.next(node)
            i+=1
        mi=node.mi
        inps[i,:]=mi.get_input(self.snapshot)
        np.random.shuffle(v)
        node=self.mis.head
        i=0
        while not node.is_tail:
            mi=node.mi
            inp=inps[i,:]
            mi.brain.control(mi,inp)
            node=self.mis.next(node)
            i+=1
        mi=node.mi
        inp=inps[i,:]
        mi.brain.control(mi,inp)
    cdef void act_all(self):
        cdef MinionDLLNode node
        cdef Minion mi
        if self.mis.len==0:
            return
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            if not mi.dead:
                self.act(mi)
            node=self.mis.next(node)
        mi=node.mi
        if not mi.dead:
            self.act(mi)
    cdef void basal_metabolism_all(self):
        cdef MinionDLLNode node
        cdef Minion mi
        if self.mis.len==0:
            return
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            if not mi.dead:
                self.basal_metabolism(mi)
            node=self.mis.next(node)
        mi=node.mi
        if not mi.dead:
            self.basal_metabolism(mi)
    cdef void render(self):
        cdef size_t xsize,ysize,x,y,i,j
        cdef int xpos,ypos
        cdef MinionDLLNode node
        cdef Minion mi
        cdef (int,int,int) c
        xsize=<size_t>self.xsize
        ysize=<size_t>self.ysize
        for x in range(xsize):
            for y in range(ysize):
                self.snapshot[x,y]=0
        if self.mis.len==0:
            pass
        else:
            node=self.mis.head
            while not node.is_tail:
                mi=node.mi
                c=mi.color
                xpos=mi.pos[0]
                ypos=mi.pos[1]
                for i in range(1+2*<size_t>mi.alen):
                    x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                    for j in range(1+2*<size_t>mi.alen):
                        y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                        self.snapshot[x,y]=256**2*c[0]+256*c[1]+c[2]
                node=self.mis.next(node)
            mi=node.mi
            c=mi.color
            xpos=mi.pos[0]
            ypos=mi.pos[1]
            for i in range(1+2*<size_t>mi.alen):
                x=<size_t>((xpos-mi.alen+<int>i)%self.xsize)
                for j in range(1+2*<size_t>mi.alen):
                    y=<size_t>((ypos-mi.alen+<int>j)%self.ysize)
                    self.snapshot[x,y]=256**2*c[0]+256*c[1]+c[2]
        for x in range(xsize):
            for y in range(ysize):
                if self.mins[x,y]>0:
                    self.snapshot[x,y]=256**3-1

        #for i in range(self.xsize):
         #   print(list(self.snapshot[i,:]))
    #public API helpers
    cdef double _total_minion_mass(self):
        cdef double total
        cdef MinionDLLNode node
        cdef Minion mi
        total=0
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            if not mi.dead:
                total+=mi.mass
            node=self.mis.next(node)
        mi=node.mi
        if not mi.dead:
            total+=mi.mass
        return total
    cdef double _total_mass(self):
        cdef double total
        cdef size_t xsize,ysize,x,y
        cdef MinionDLLNode node
        cdef Minion mi
        total=self.hidden_mass

        xsize=<size_t>self.xsize
        ysize=<size_t>self.ysize
        for x in range(xsize):
            for y in range(ysize):
                total+=self.mins[x,y]
        if self.mis.len==0:
            return total

        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            if not mi.dead:
                total+=mi.mass
            node=self.mis.next(node)
        mi=node.mi
        if not mi.dead:
            total+=mi.mass
        return total
    cdef double _total_min(self):
        cdef double total
        cdef size_t x,y
        total=0
        for x in range(<size_t>self.xsize):
            for y in range(<size_t>self.ysize):
                total+=self.mins[x,y]
        return total
    cdef bint _exhausted(self):
        cdef size_t x,y
        for x in range(<size_t>self.xsize):
            for y in range(<size_t>self.ysize):
                if self.mins[x,y]>0:
                    return False
        return True
    cdef double _get_avg_r(self):
        if self.mis.len==0:
            return 0
        cdef double sum
        cdef MinionDLLNode node
        cdef Minion mi
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            sum+=<double>mi.color[0]
            node=self.mis.next(node)
        mi=node.mi
        sum+=<double>mi.color[0]
        return sum/(<double>self.mis.len)
    cdef double _get_avg_g(self):
        if self.mis.len==0:
            return 0
        cdef double sum
        cdef MinionDLLNode node
        cdef Minion mi
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            sum+=<double>mi.color[1]
            node=self.mis.next(node)
        mi=node.mi
        sum+=<double>mi.color[1]
        return sum/(<double>self.mis.len)
    cdef double _get_avg_b(self):
        if self.mis.len==0:
            return 0
        cdef double sum
        cdef MinionDLLNode node
        cdef Minion mi
        node=self.mis.head
        while not node.is_tail:
            mi=node.mi
            sum+=<double>mi.color[2]
            node=self.mis.next(node)
        mi=node.mi
        sum+=<double>mi.color[2]
        return sum/(<double>self.mis.len)

    #public API-------------------------------------------
    #main
    def occupy_check(self):
        for x in range(self.xsize):
            for y in range(self.ysize):
                dll=self.occupy_map[x,y]
                for mi in dll:
                    xpos,ypos=mi.get_pos()
                    if max(dist(xpos,x,self.xsize),dist(ypos,y,self.ysize))>mi.get_alen():
                        return 1,False,mi.get_pos(),mi.get_alen(),x,y
        for mi in self.mis:
            xpos,ypos=mi.get_pos()
            for a in range(-mi.get_alen(),1+mi.get_alen()):
                for b in range(-mi.get_alen(),1+mi.get_alen()):
                    x=(xpos+a)%self.xsize
                    y=(ypos+b)%self.ysize
                    dll=self.occupy_map[x,y]
                    if not dll.contains_py(mi):
                        return 2,False,mi.get_pos(),mi.get_alen(),x,y
        return True
    def total_mass_error_detector(self,prompt):
        if self.total_mass()>15000.2:
            print("error during "+prompt+"!!!")
            raise Exception()
    def mins_negative_error_detector(self,prompt):
        for i in range(self.get_xsize()):
            for j in range(self.get_ysize()):
                if self.get_mins()[i,j]<0:
                    print("mins negative in "+prompt+" at ",i,j)
                    return
    def snapshot_error_detector(self,prompt):
        for i in range(self.get_xsize()):
            for j in range(self.get_ysize()):
                if self.snapshot[i,j]==1:
                    print("snapshot not assigned in "+prompt)
                    return
        print("no problem in "+prompt)
    def max_mass(self):
        m=-1
        for mi in self.mis:
            if mi.get_mass()>m:
                m=mi.get_mass()
        return m
    def death_by(self,prompt):
        print("death by "+prompt)
    def status(self):
        red_num=0
        blue_num=0
        for mi in self.mis:
            if mi.get_dna().get_colorTrait().get_c()==(255,0,0):
                red_num+=1
            if mi.get_dna().get_colorTrait().get_c()==(0,0,255):
                blue_num+=1
        return (red_num,blue_num)
    def append_red_dnas(self,l):
        for mi in self.mis:
            if mi.get_dna().get_colorTrait().get_c()==(255,0,0):
                l.append(mi.get_dna())
    def append_blue_dnas(self,l):
        for mi in self.mis:
            if mi.get_dna().get_colorTrait().get_c()==(0,0,255):
                l.append(mi.get_dna())
    def one_step(self):
        self.render()
        self.control_all()
        self.act_all()
        self.basal_metabolism_all()
        self.moment+=1
        self.kill_elderly()
        self.unregister_deads()
        if self.moment%100==0:
            seed()
    #basic
    def get_xsize(self):
        return self.xsize
    def get_ysize(self):
        return self.ysize
    def get_moment(self):
        return self.moment
    def get_new_id(self):
        return self.new_id
    def get_mins(self):
        return self.mins
    def total_mass(self):
        return self._total_mass()
    def total_minion_mass(self):
        return self._total_minion_mass()
    def total_min(self):
        return self._total_min()
    #end conditions
    def limit(self,n):
        return self.moment>=n
    def exhausted(self):
        return self._exhausted()
    def extincted(self):
        return self.mis.is_empty()
    #measurements
    def get_pedigree(self):
        return list(self.pedigree)
    def sample_dna(self,n):
        if n>=self.mis.len:
            subpopulation=sample(list(self.mis),n)
        else:
            subpopulation=list(self.mis)
        dnas=[]
        for mi in subpopulation:
            dnas.append(mi.get_dna())
        return dnas


    def get_mis(self):
        return list(self.mis)
    def get_nth_mi(self,n):
        return self.mis.get_nth_py(n)
    def get_population(self):
        return self.mis.len
    def get_colors(self):
        for mi in self.mis:
            yield mi.color
    def get_avg_r(self):
        return self._get_avg_r()
    def get_avg_g(self):
        return self._get_avg_g()
    def get_avg_b(self):
        return self._get_avg_b()
    def get_messiness(self):
        return self.messiness

    #for use in visualization
    def body_rect(self,k):
        #called in: play.py
        for mi in self.mis:
            x,y=mi.get_pos()
            alen=mi.get_alen()
            yield mi.get_color(),(k*(x-alen),k*(y-alen),k*(1+2*alen),k*(1+2*alen))
    def geni_rect(self,k):
        #called in: play.py
        for mi in self.mis:
            if mi.get_action()==1:
                x,y=mi.get_pos()
                alen=mi.get_alen()
                color=mi.get_color()
                move_direc=mi.get_move_direc()
                if move_direc==0:
                    yield color,(k*x,k*y,k*(1+2*alen),k)
                elif move_direc==1:
                    yield color,(k*x,k*y,k,k*(1+2*alen))
                elif move_direc==2:
                    yield color,(k*x+k-1,k*y,-k*(1+2*alen),k)
                elif move_direc==3:
                    yield color,(k*x,k*y+k-1,k,-k*(1+2*alen))
