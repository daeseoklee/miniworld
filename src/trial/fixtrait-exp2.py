import exputil_cy
import world_cy as world
import random
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("n",type=int)
args=parser.parse_args()

merge_thres = 0.04
mut_per_diff = 3.0

xsize, ysize, k = 400, 300, 4
limit = 100000
total_mass = 40000

ininum_red,ininum_blue=120,240
red_maxsize,red_uptake,red_maxage=6,0.6,1600
blue_maxsize,blue_uptake,blue_maxage=4,0.4,1600
avg_consum_rate,consum_exp=0.01,0.5

exputil_cy.set_size(xsize, ysize)
exputil_cy.set_heredity(merge_thres, mut_per_diff)
world.set_consum(avg_consum_rate,consum_exp)

thres=0.1

filename="./results/result"+str(args.n)+".txt"


def dofun(w,red_keep,blue_keep,tmp_red_dnas,tmp_blue_dnas):
    w.one_step()
    r,b=w.status()
    print(r,b)
    with open(filename,'a') as f:
        f.write(str(w.get_moment())+" - "+str(r)+" , "+str(b)+"\n")
    if not red_keep[0]:
        if r<=ininum_red*thres:
            tmp_red_dnas.clear()
            w.append_red_dnas(tmp_red_dnas)
            red_keep[0]=True
    else:
        if r>ininum_red*thres:
            red_keep[0]=False

    if not blue_keep[0]:
        if b<=ininum_blue*thres:
            tmp_blue_dnas.clear()
            w.append_blue_dnas(tmp_blue_dnas)
            blue_keep[0]=True
    else:
        if b>ininum_blue*thres:
            blue_keep[0]=False

def endf(w):
    r,b=w.status()
    return (r==0 or b==0) or w.get_moment()==limit
def returnf(w):
    r,b=w.status()
    if r==0 and b>0:
        return 'blue'
    elif b==0 and r>0:
        return 'red'
    return 'none'




for iternum in range(200):
    if iternum!=0:
        print(win)
    with open(filename,'a') as f:
        f.write("#"+str(iternum)+"th iteration-----------------------------------\n")
    print(str(iternum)+"th iteration-----------------------------------")
    if iternum==0:
        tmp_red_dnas, tmp_blue_dnas = [], []
        red_dnas, blue_dnas = [], []
        red_dnas=[world.randLinearDNA_with(c=(255,0,0),mst_a=red_maxsize,utt_a=red_uptake,mat_a=red_maxage)]
        blue_dnas=[world.randLinearDNA_with(c=(0,0,255),mst_a=blue_maxsize,utt_a=blue_uptake,mat_a=blue_maxage)]
    elif win=="red":
        red_dnas.clear()
        w.append_red_dnas(red_dnas)
        random.shuffle(red_dnas)
        blue_dnas.clear()
        for dna in tmp_blue_dnas:
            blue_dnas.append(dna)
    elif win=="blue":
        blue_dnas.clear()
        w.append_blue_dnas(blue_dnas)
        random.shuffle(blue_dnas)
        red_dnas.clear()
        for dna in tmp_red_dnas:
            red_dnas.append(dna)
    elif win=="none":
        red_dnas.clear()
        w.append_red_dnas(red_dnas)
        random.shuffle(red_dnas)
        blue_dnas.clear()
        w.append_blue_dnas(blue_dnas)
        random.shuffle(blue_dnas)
    mis=[]
    for i in range(ininum_red):
        mis.append(world.construct_minion(red_dnas[i%len(red_dnas)],alen=1,pos=(-1,-1),do_freeze=False))
    for i in range(ininum_blue):
        mis.append(world.construct_minion(blue_dnas[i%len(blue_dnas)],alen=1,pos=(-1,-1),do_freeze=False))
    mis=world.dll_from_list(mis)
    w=world.World(xsize,ysize,total_mass,mis,no_age=False,no_birth=False,no_eat=False,\
              no_energy=False,no_excrete=False,no_hunt=False,halluc=False,record_pedigree=False)
    red_keep,blue_keep=[False],[False]
    tmp_red_dnas.clear()
    w.append_red_dnas(tmp_red_dnas)
    tmp_blue_dnas.clear()
    w.append_blue_dnas(tmp_blue_dnas)
    dof=lambda w:dofun(w,red_keep,blue_keep,tmp_red_dnas,tmp_blue_dnas)
    win=exputil_cy.do(w, k, dof=dof, endf=endf,
                  returnf=returnf, screen_off=True, print_moment=True, print_moment_period=1, screen_period=50)


