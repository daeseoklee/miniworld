import sys
sys.path.append("./core")
#import exputil_py
import exputil_py

xsize,ysize,k=200,150,5
num_dna=2
ininum=120
limit=40000
period=500
total_mass=10000

merge_thres=0.04
mut_per_diff=3.0

w=exputil_py.construct_world(xsize,ysize,total_mass,num_dna,ininum)
exputil_py.do(w,k,dof=lambda w:w.one_step(),endf=lambda w:w.extincted() or w.limit(limit),returnf=lambda w:None,screen_off=False,print_moment=True,print_moment_period=1)
