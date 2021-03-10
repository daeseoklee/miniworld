*intro 
This repository contains miniworld 0.4.1(which will be soon documented) implementation,
codes that generated the data cited in "NATURAL EVOLUTION OF BEHAVIORS IN A VIRTUAL
ECOSYSTEM AND ITS VALIDATION", and the data.
An implementation of the core system of miniworld 0.4.1, which is called "the virtual ecosystem" in the paper,
is in the directory 'src/core'.
There are both python and cython implementations. 
The former makes use of multiple cpu cores(maybe due to numpy), while the latter does not.
Normalizing the speed by the cpu usage, the latter was about 10x more efficient in a test. 
Most of the data-generating codes are in 'src'. 
The data are in 'data'. 
If you adjust the *_visual.py file in each subidrectory of 'data' a little bit approriately
and run it, you will get the data appearing in the paper.   

*quick test 
To ensure that you can replicate the simulation and experiments with your own computer, 
you can do a quick test as follows. 
<If you use linux> 
Run 
"python quick_py.py" 
"python quick_cy.py" 
<If you use another os> 
Run 
"python quick_py.py" 
For making use of the cython implementation, 
you might need to 
1. run src/core/setup.py with your own cython 
and c compiler 
2. Adjust the directory names in .py files
before running "python quick_cy.py" 

*dependency 
python 3.x.x
<python library> 
cython
numpy 1.16.(>=4)
pygame 1.9.(>=6)
matplotlib 3.1.x
