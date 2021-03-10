import sys
sys.path.append("./core")
import exputil_cy

xsize,ysize,k=200,150,5
num_dna=2
ininum=120
limit=40000
period=500
total_mass=10000

merge_thres=0.04
mut_per_diff=3.0

exputil_cy.set_size(xsize,ysize)
exputil_cy.set_traits_range(ms1=3,ms2=8,ut1=0.1,ut2=0.7,ma1=1500,ma2=1500)
exputil_cy.set_heredity(merge_thres,mut_per_diff)

w=exputil_cy.construct_world(xsize,ysize,total_mass,num_dna,ininum)
exputil_cy.do(w,k,dof=lambda w:w.one_step(),endf=lambda w:w.extincted() or w.limit(limit),returnf=lambda w:None,screen_off=False,print_moment=True,print_moment_period=1)
