from matplotlib import pyplot as plt

fileformat="./results/performance{}.txt"

resultd=dict()

for i in range(1,31):
    filename=fileformat.format(i)
    with open(filename,'r') as f:
        for line in f.readlines():
            if not i in resultd:
                resultd[i]=[float(line)]
            else:
                resultd[i].append(float(line))
l=[]
for j in range(33):
    m=0
    for i in range(1,31):
        if resultd[i][j]>m:
            m=resultd[i][j]
    l.append(m)
plt.plot(l)
plt.show()