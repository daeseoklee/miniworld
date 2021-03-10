import world_cy
import pygame
from fast_random import seed_py

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

def draw(w,screen,k):
    screen.fill((0, 0, 0))
    mins=w.get_mins()
    for color,rect in w.body_rect(k):
        # pygame.draw.rect(screen,color,rect)
        draw_rects(screen, color, rect)
    for color,rect in w.geni_rect(k):
        # pygame.draw.rect(screen,color,rect)
        draw_rects(screen, color,rect)
    for i in range(w.get_xsize()):
        for j in range(w.get_ysize()):
            if mins[i][j] > 0:
                pygame.draw.rect(screen, (255, 255, 255), (k * i, k * j, k, k))
    pygame.display.flip()

#public API-------------------------------------------------------
def set_size(x,y):
    world_cy.set_size(x,y)

def set_traits_range(ms1,ms2,ut1,ut2,ma1,ma2):
    world_cy.set_traits_range(ms1,ms2,ut1,ut2,ma1,ma2)
def set_heredity(new_merge_thres,new_mut_per_diff):
    world_cy.set_heredity(new_merge_thres, new_mut_per_diff)
def set_vision(new_vision_range,new_vision_resolution):
    world_cy.set_vision(new_vision_range, new_vision_resolution)

def get_heredity():
    return world_cy.merge_thres, world_cy.mut_per_diff




def dna_with_opposite_color(dna):
    colorTrait,maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait=dna.get_traits()
    return world_cy.LinearDNA(colorTrait.opposite(),maxsizeTrait,uptakeTrait,maxageTrait,weightsTrait)




def construct_world(xsize,ysize,total_mass,num_dna,ininum,\
             no_age=False,no_birth=False,no_eat=False, no_energy=False, no_excrete=False,\
             no_hunt=False,halluc=False,alen_list=None,pos_list=None,freeze_list=None,\
             record_pedigree=False):
    seed_py()
    dnas=[world_cy.randLinearDNA() for _ in range(num_dna)]
    mis=world_cy.dll_from_list(\
        [world_cy.construct_minion(dnas[i%num_dna],\
        alen=1 if alen_list is None else alen_list[i],\
        pos=(-1,-1) if pos_list is None else pos_list[i],\
        do_freeze=False if freeze_list is None else freeze_list[i])\
         for i in range(ininum)])
    return world_cy.World(xsize,ysize,total_mass,mis,no_age=no_age,no_birth=no_birth,\
                       no_eat=no_eat,no_energy=no_energy,no_excrete=no_excrete,no_hunt=no_hunt,\
                       halluc=halluc,record_pedigree=record_pedigree)

def do(w,k,dof,endf,returnf,screen_off=False,print_moment=True,print_moment_period=1,print_script=""):
    if screen_off:
        while not endf(w):
            dof(w)
            if print_moment and w.get_moment() % print_moment_period == 0:
                print(print_script,"-","moment: ", w.get_moment())
                print("population: ",w.get_population())
        return returnf(w)
    else:
        screen=pygame.display.set_mode((w.get_xsize()*k,w.get_ysize()*k))
        screen.fill((0,0,0))
        while True:
            for _ in range(5):
                event=pygame.event.poll()
                if event.type==pygame.QUIT:
                    pygame.quit()
                    return returnf(w)
            else:
                if not endf(w):
                    dof(w)
                else:
                    pygame.quit()
                    return returnf(w)
                if print_moment and w.get_moment() % print_moment_period == 0:
                    print(print_script,"-","moment: ", w.get_moment())
                    print("population: ", w.get_population())
                draw(w,screen,k)



