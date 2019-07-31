import minion
import world
import pygame

def main():
    xsize,ysize=80,50
    k=10
    total_mass=800
    num=40
    dnas=[minion.LinearDNA((255,0,0),None,init=True),minion.LinearDNA((0,255,0),None,init=True),minion.LinearDNA((0,0,255),None,init=True)]
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
            screen.fill((0,0,0))
            for i in range(xsize):
                for j in range(ysize):
                    if w.mins[(i,j)]>0:
                        pygame.draw.rect(screen,(255,255,255),(k*i,k*j,k,k))
            for mi in w.mis:
                #pygame.draw.rect(screen,mi.c,w.body_rect(mi,k))
                draw_rects(screen, mi.c, w.body_rect(mi, k))
                if mi.move_dist<-0.5:
                    #pygame.draw.rect(screen,mi.c,w.fork_rect(mi,k))
                    draw_rects(screen, mi.c, w.fork_rect(mi, k))
            pygame.display.flip()

def draw_rects(surface,c,rect):
    x,y,w,h=rect
    xsize,ysize=surface.get_width(),surface.get_height()
    if w<0:
        draw_rects(surface,c,(x+w,y,-w,h))
    elif h<0:
        draw_rects(surface,c,(x,y+h,w,-h))
    elif x<0:
        draw_rects(surface,c,(x+xsize,y,-x,h))
        draw_rects(surface,c,(0,y,x+w,h))
    elif y<0:
        draw_rects(surface,c,(x,y+ysize,w,-y))
        draw_rects(surface,c,(x,0,w,y+h))
    elif x+w-1>=xsize:
        draw_rects(surface,c,(x,y,xsize-x,h))
        draw_rects(surface,c,(0,y,w-xsize+x,h))
    elif y+h-1>=ysize:
        draw_rects(surface,c,(x,y,w,ysize-y))
        draw_rects(surface,c,(x,0,w,h-ysize+y))
    else:
        pygame.draw.rect(surface,c,rect)


def step(w):
    w.one_step()
    #print(w.pedigree)

main()

