import exputil_cy
import world_cy as world
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("n",type=int)
args=parser.parse_args()

merge_thres = 0.04
mut_per_diff = 3.0

xsize, ysize, k = 400, 300, 4
total_mass = 29000+2000*args.n

ininum_red,ininum_blue=120,240
red_maxsize,red_uptake,red_maxage=6,0.6,1600
blue_maxsize,blue_uptake,blue_maxage=4,0.4,1600

resultdir="./results/"
n,m=9,7
num_iter=2
criterion=100000

def endf(w):
    r,b=w.status()
    return (r==0 or b==0) or w.limit(criterion)
def dofun(w,rs,bs):
    r,b=w.status()
    rs.append(r)
    bs.append(b)
    w.one_step()
def returnf(w):
    r,b=w.status()
    return w.get_moment(),w.get_population(),r,b


for i in range(n):
    consum_exp=0.1*i
    for j in range(m):
        for iternum in range(num_iter):
            avg_consum_rate=0.025*(1.5**(j-3))/(5.5**consum_exp)
            exputil_cy.set_size(xsize, ysize)
            exputil_cy.set_heredity(merge_thres, mut_per_diff)
            world.set_consum(avg_consum_rate,consum_exp)
            red_dna = world.randLinearDNA_with(c=(255, 0, 0), mst_a=red_maxsize, utt_a=red_uptake, mat_a=red_maxage)
            blue_dna = world.randLinearDNA_with(c=(0, 0, 255), mst_a=blue_maxsize, utt_a=blue_uptake, mat_a=blue_maxage)
            mis = []
            for _ in range(ininum_red):
                mis.append(world.construct_minion(red_dna, alen=1, pos=(-1, -1), do_freeze=False))
            for _ in range(ininum_blue):
                mis.append(world.construct_minion(blue_dna, alen=1, pos=(-1, -1), do_freeze=False))
            print(len(mis))
            mis = world.dll_from_list(mis)
            w = world.World(xsize, ysize, total_mass, mis, no_age=False, no_birth=False, no_eat=False, \
                            no_energy=False, no_excrete=False, no_hunt=False, halluc=False, record_pedigree=False)
            rs, bs = [], []
            dof=lambda w:dofun(w,rs,bs)
            moment,population,r,b=exputil_cy.do(w, k, dof, endf, returnf, screen_off=True, print_moment=True, print_moment_period=1,
                          screen_period=50)
            with open(resultdir+"result"+str(args.n)+".txt",'a') as f:
                f.write(str(i)+" , "+str(j)+" - "+str(moment)+" , "+str(r)+" , "+str(b)+" , "+str(population-r-b)+"\n")
            with open(resultdir+"population_"+str(args.n)+"_"+str(i)+"_"+str(j)+"_"+str(iternum)+".txt",'w') as f:
                for r,b in zip(rs,bs):
                    f.write(str(r)+" , "+str(b)+"\n")
                f.close()



