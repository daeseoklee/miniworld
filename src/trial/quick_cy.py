import exputil_cy
import world_cy as world

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

red_dna=world.randLinearDNA_with(c=(255,0,0),mst_a=red_maxsize,utt_a=red_uptake,mat_a=red_maxage)
blue_dna=world.randLinearDNA_with(c=(0,0,255),mst_a=blue_maxsize,utt_a=blue_uptake,mat_a=blue_maxage)
mis=[]
for _ in range(ininum_red):
    mis.append(world.construct_minion(red_dna,alen=1,pos=(-1,-1),do_freeze=False))
for _ in range(ininum_blue):
    mis.append(world.construct_minion(blue_dna,alen=1,pos=(-1,-1),do_freeze=False))
print(len(mis))
mis=world.dll_from_list(mis)

w=world.World(xsize,ysize,total_mass,mis,no_age=False,no_birth=False,no_eat=False,\
              no_energy=False,no_excrete=False,no_hunt=False,halluc=False,record_pedigree=False)


def dof(w):
    w.one_step()
    print(w.get_population())
    ct,mst,utt,mat,wt=w.get_nth_mi(0).get_dna().get_traits()
    print(ct.get_c(),mst.get_a(),utt.get_a(),mat.get_a())
exputil_cy.do(w, k, dof=dof, endf=lambda w: w.extincted() or w.limit(limit),
              returnf=lambda w: None, screen_off=True, print_moment=True, print_moment_period=1,screen_period=50)
