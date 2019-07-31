import minion
import world
import pygame

def main():
    xsize,ysize=100,70
    k=5
    total_mass=2000
    num=100
    dnas=[minion.FoolDNA((255,0,0)),minion.FoolDNA((0,255,0)),minion.FoolDNA((0,0,255))]
    mis=[minion.Minion(dnas[i%3]) for i in range(num)]
    w=world.World(xsize,ysize,total_mass,mis)

    screen:pygame.Surface=pygame.display.set_mode((xsize*k,ysize*k))
    screen.fill((0,0,0))
    running=1
    moment=0
    while running:
        event=pygame.event.poll()
        if event.type==pygame.QUIT:
          running=0
        else:
            step(w)
            moment+=1
            print("moment: ",moment)
            for i in range(xsize):
                for j in range(ysize):
                    for a in range(k):
                        for b in range(k):
                            screen.set_at((k*i+a,k*j+b),w.snapshot[i][j])
            pygame.display.flip()
def step(w):
    w.one_step()
    print(w.pedigree)

main()

