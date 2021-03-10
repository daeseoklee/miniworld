import sys
sys.path.append("./core")
import exputil_cy

xsize, ysize, k = 300, 200, 5
num_dna = 3
ininum = 120
limit = 1100
total_mass = 20000

merge_thres = 0.04
mut_per_diff = 3.0

exputil_cy.set_size(xsize, ysize)
exputil_cy.set_traits_range(ms1=3, ms2=8, ut1=0.1, ut2=0.7, ma1=1500, ma2=1500)
exputil_cy.set_heredity(merge_thres, mut_per_diff)

w = exputil_cy.construct_world(xsize, ysize, total_mass, num_dna, ininum,record_pedigree=True)


def get_num_child(pedigree,id,end_moment):
    s={id}
    count=0
    for mid,fid,cid,moment in pedigree:
        if moment>end_moment:
            break
        if mid in s:
            count+=1
            s.add(cid)
            continue
        if fid in s:
            count+=1
            s.add(cid)
    return count

def returnf(w):
    pedigree=w.get_pedigree()
    print(pedigree)
    return
    for id in range(360):
        print(get_num_child(pedigree,id,1000))
    """
    id1=pedigree[0][0]
    id2=pedigree[0][1]
    id3=pedigree[1][0]
    id4=pedigree[2][1]
    print("first dad until 1000:",get_num_child(pedigree,id1,1000))
    print("first dad until 1500:", get_num_child(pedigree, id1, 1500))
    print("first mom until 1000:", get_num_child(pedigree, id2, 1000))
    print("first mom until 1500:", get_num_child(pedigree, id2, 1500))
    print("second dad until 1000:", get_num_child(pedigree, id3, 1000))
    print("second dad until 1500:", get_num_child(pedigree, id3, 1500))
    print("second mom until 1000:", get_num_child(pedigree, id4, 1000))
    print("second mom until 1500:", get_num_child(pedigree, id4, 1500))
    """


exputil_cy.do(w, k, dof=lambda w: w.one_step(), endf=lambda w: w.extincted() or w.limit(limit),
              returnf=returnf, screen_off=True, print_moment=True, print_moment_period=10)

