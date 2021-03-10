import sys
sys.path.append("./core")
import exp1


ininum=0
num_run=10
limit=300000
thres=5
population_period=30000
iter_per_sample=1
test_limit=400


def alternative():
    for runnum in range(0,1):
        for testname in ["food2(n)","food2(h)"]:
            exp1.test(testname,runnum,iter_per_sample,limit=test_limit,thres=thres,index=0,screen_off=True)



def main():
    for runnum in range(ininum,ininum+num_run):
        exp1.do_init(limit,population_period,runnum,screen_off=True,print_moment=True,print_moment_period=100)
    for runnum in range(ininum,ininum+num_run):
        for testname in ["rape(o)","rape(e)","food(n)","food(h)","hunt(o)","hunt(e)","messy(n)","messy(h)","birth(n)","birth(h)"]:
            exp1.test(testname,runnum,iter_per_sample,limit=test_limit,thres=thres,index=0,screen_off=True)
alternative()

