import exputil_cy as exputil
import world_cy as world
import random
import util
import pickle
import os
from fast_random_py import randint



basedir=os.getcwd()
resultdir=basedir+"/results/"
xsize, ysize, k = 300, 200, 5
total_mass = 20000
test_xsize,test_ysize,test_k=30,20,15
test_total_mass=200
minmaxsize,maxmaxsize=5,5
minuptake,maxuptake=0.4,0.4
minmaxage,maxmaxage=1600,1600
num_dna = 3
ininum = 150


merge_thres = 0.04
mut_per_diff = 3.0

exputil.set_size(xsize, ysize)
exputil.set_traits_range(ms1=minmaxsize, ms2=maxmaxsize, ut1=minuptake, ut2=maxuptake, ma1=minmaxage, ma2=maxmaxage)
exputil.set_heredity(merge_thres, mut_per_diff)



#>>>initial run>>>
def init_world():
    return exputil.construct_world(xsize, ysize, total_mass, num_dna, ininum)

def do_init(limit,period,runnum,screen_off=False,print_moment=True,print_moment_period=1000):
    """

    :param limit: moment limit
    :param period:
    :param runnum: run id
    :return:
    """
    w=init_world()
    times=[]
    pc_list=[]
    samples=[]
    dof=lambda w:pc_and_sample_dofun(w,pc_list,times,samples,period)
    endf=lambda w:limit_endfun(w,limit)
    returnf=lambda w:pc_and_sample_returnfun(pc_list,times,samples,runnum)
    exputil.do(w,k,dof,endf,returnf,screen_off=screen_off,print_moment=print_moment,print_moment_period=print_moment_period,print_script=str(runnum)+"th run")
def pc_and_sample_dofun(w,pc_list,times,samples,period):
    pc_list.append((w.get_population(),w.get_avg_r(),w.get_avg_g(),w.get_avg_b()))
    if w.get_moment()%period==0:
        times.append(w.get_moment())
        samples.append(list(map(lambda mi:mi.get_dna(),random.sample(list(w.get_mis()),min(20,w.get_population())))))
    w.one_step()
def pc_and_sample_returnfun(pc_list,times,samples,runnum):
    util.write_results(resultdir,runnum,0,"pc",pc_list,list(range(0,len(pc_list))))
    filename=resultdir+str(runnum)+"_sample"
    world.write_list_of_dnas_file(samples, filename)
    filename=resultdir+str(runnum)+"_time"
    f=open(resultdir+str(runnum)+"_time",'wb')
    pickle.dump(times,f)
    f.close()
#<<<initial run<<<

#>>>subsequent tests>>>
def testtype(testname):
    if testname in ["hunt(o)","hunt(e)","rape(o)","rape(e)"]:
        return 2
    return 1
def test_dofun(w,testname):
    w.one_step()
def test_endfun(w,testname,limit=None,thres=None):
    if testname in ["food(n)","food(h)"]: #(n) stands for 'normal' and (h) stands for 'hallucinated'
        return w.total_min()<0.1*w.total_mass()
    elif testname in ["hunt(o)","hunt(e)"]:
        return w.get_population()<2
    elif testname in ["rape(o)","rape(e)"]:
        return w.get_messiness()>=thres
    elif testname in ["birth(n)","birth(h)"]:
        return w.get_moment()>=limit or w.get_population()==0
    elif testname in ["messy(n)","messy(h)"]:
        return w.get_moment()>=limit


def test_returnfun(w,testname,index=None,ininum=40):
    if testname in ["food(n)","food(h)"]:
        return w.get_moment()
    elif testname in ["hunt(o)","hunt(e)"]:
        return w.get_nth_mi(index).get_cum_dist()
    elif testname in ["rape(o)","rape(e)"]: #o:homo, e:hetero
        return w.get_nth_mi(index).get_cum_dist()
    elif testname in ["birth(n)","birth(h)"]:
        return w.get_new_id()-ininum
    elif testname in ["messy(n)","messy(h)"]:
        return w.get_messiness()
def furthest_color_dna_pair(dnas):
    first=None
    second=None
    m=0
    for dna1 in dnas:
        for dna2 in dnas:
            if world.ColorTrait.get_normalized_difference(dna1.get_colorTrait(),dna2.get_colorTrait())>m:
                first=dna1
                second=dna2
                m=world.ColorTrait.get_normalized_difference(dna1.get_colorTrait(),dna2.get_colorTrait())
    return (first,second)

def test_world(testname,dnas):
    if testname=="food(n)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)], \
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=False, no_energy=True, no_excrete=True, no_hunt=True, \
                              halluc=False, record_pedigree=False)


    elif testname=="food(h)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)],
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=False, no_energy=True, no_excrete=True, no_hunt=True, \
                              halluc=True, record_pedigree=False)
    elif testname=="hunt(o)":
        dna1,_=furthest_color_dna_pair(dnas)
        mis = world.dll_from_list( \
            [world.construct_minion(dna1,alen=4,pos=(0,0), do_freeze=False), \
             world.construct_minion(dna1, alen=1, pos=(test_xsize//2,test_ysize//2),do_freeze=True)])
        return world.World(test_xsize, test_ysize, test_total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=True, no_energy=True, no_excrete=True, no_hunt=False, \
                              halluc=False, record_pedigree=False)

    elif testname=="hunt(e)":
        dna1,dna2=furthest_color_dna_pair(dnas)
        mis = world.dll_from_list( \
            [world.construct_minion(dna1,alen=4,pos=(0,0), do_freeze=False), \
             world.construct_minion(dna2, alen=1, pos=(test_xsize//2,test_ysize//2),do_freeze=True)])
        return world.World(test_xsize, test_ysize, test_total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=True, no_energy=True, no_excrete=True, no_hunt=False, \
                              halluc=False, record_pedigree=False)

    elif testname=="rape(o)":
        dna1,_=furthest_color_dna_pair(dnas)
        mis = world.dll_from_list( \
            [world.construct_minion(dna1,alen=1,pos=(0,0), do_freeze=False), \
             world.construct_minion(dna1, alen=4, pos=(test_xsize//2,test_ysize//2),do_freeze=True)])
        return world.World(test_xsize, test_ysize, test_total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=True, no_energy=True, no_excrete=True, no_hunt=True, \
                              halluc=False, record_pedigree=False)
    elif testname=="rape(e)":
        dna1,dna2=furthest_color_dna_pair(dnas)
        mis = world.dll_from_list( \
            [world.construct_minion(dna1,alen=1,pos=(0,0), do_freeze=False), \
             world.construct_minion(dna2, alen=4, pos=(test_xsize//2,test_ysize//2),do_freeze=True)])
        return world.World(test_xsize, test_ysize, test_total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=True, no_energy=True, no_excrete=True, no_hunt=True, \
                              halluc=False, record_pedigree=False)
    elif testname=="birth(n)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)], \
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=False, no_birth=False, \
                              no_eat=False, no_energy=False, no_excrete=False, no_hunt=False, \
                              halluc=False, record_pedigree=False)
    elif testname=="birth(h)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)], \
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=False, no_birth=False, \
                              no_eat=False, no_energy=False, no_excrete=False, no_hunt=False, \
                              halluc=True, record_pedigree=False)
    elif testname=="messy(n)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)], \
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=False, no_energy=True, no_excrete=False, no_hunt=True, \
                              halluc=False, record_pedigree=False)
    elif testname=="messy(h)":
        mis = world.dll_from_list( \
            [world.construct_minion(dnas[i % len(dnas)], \
            alen=1,pos=(randint(0,xsize-1),randint(0,ysize-1)), do_freeze=False)\
             for i in range(ininum)])
        return world.World(xsize, ysize, total_mass, mis, no_age=True, no_birth=True, \
                              no_eat=False, no_energy=True, no_excrete=False, no_hunt=True, \
                              halluc=True, record_pedigree=False)

def test(testname,runnum,iter_per_sample,limit=None,thres=None,index=None,screen_off=True):
    print("doing",runnum,"-",testname)
    filename=resultdir+str(runnum)+"_time"
    f=open(filename,'rb')
    times=pickle.load(f)
    print("times:",times)
    f.close()
    filename=resultdir+str(runnum)+"_sample"
    samples=world.read_list_of_dnas_file(filename)
    assert len(times)==len(samples)
    n=len(times)
    for iteration in range(iter_per_sample):
        results=[]
        for i,(t,dnas) in enumerate(zip(times,samples)): #for each moment
            w=test_world(testname,dnas)
            dof=lambda w:test_dofun(w,testname)
            endf=lambda w:test_endfun(w,testname,limit=limit,thres=thres)
            returnf=lambda w:test_returnfun(w,testname,index=index,ininum=ininum)
            results.append(exputil.do(w,test_k if testtype(testname)==2 else k,dof,endf,returnf,screen_off=screen_off,print_moment=False))
            print("done",runnum,t,"-",iteration)
        util.write_results(resultdir,runnum,iteration,testname,results,times)
    print("done",runnum,"-",testname)
#<<<subsequent tests<<<



def limit_endfun(w,limit):
    return w.get_population()==0 or w.get_moment()>limit


