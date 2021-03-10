# distutils: language=c++

cdef class ColorTrait():
    cdef int r,g,b
    cdef void mutate(self,double rate)
    cdef ColorTrait mixed(ColorTrait t1, ColorTrait t2)
    cdef double normalized_difference(ColorTrait t1,ColorTrait t2)


cdef class PositiveTrait():
    cdef double a
    cdef double eps
    cdef void mutate(self,double rate)
    cdef PositiveTrait mixed(PositiveTrait t1,PositiveTrait t2)
    cdef double normalized_difference(PositiveTrait t1,PositiveTrait t2)

cdef class UnitTrait():
    cdef double a
    cdef double eps
    cdef void mutate(self,double rate)
    cdef UnitTrait mixed(UnitTrait t1, UnitTrait t2)
    cdef double normalized_difference(UnitTrait t1,UnitTrait t2)

cdef class FloatListTrait():
    cdef double[:] l
    cdef size_t[:] group_sizes
    cdef void mutate(self,double rate)
    cdef FloatListTrait mixed(FloatListTrait t1,FloatListTrait t2)
    cdef double normalized_difference(FloatListTrait t1,FloatListTrait t2)

cpdef LinearDNA randLinearDNA_with((int,int,int) c,double mst_a,double utt_a,double mat_a)

cpdef LinearDNA randLinearDNA()

cdef class LinearDNA():
    cdef:
        ColorTrait colorTrait
        PositiveTrait maxsizeTrait
        UnitTrait uptakeTrait
        PositiveTrait maxageTrait
        FloatListTrait weightsTrait

        bint mergeable(self,LinearDNA dna)
        LinearDNA merge(self,LinearDNA dna)
        ((int,int,int),double,double,double) translate(self)
        double[:] translate_weights(self)


cdef class Brain():
    cdef double[:] weights
    cdef void control(self,Minion mi,double[:] inp)

cdef double[:] apply_linear(double[:] v,double[:] A)

cdef class LinearBrain(Brain):
    pass


cpdef Minion construct_minion(LinearDNA dna,int alen,(int,int) pos,bint do_freeze)

cdef class Minion():
    cdef:
        (int,int,int) color
        double maxsize
        double uptake
        double maxage
        Brain brain
        int idim

        LinearDNA dna
        int id
        (int,int) pos
        int age
        bint dead
        int alen
        double mass
        int action
        int move_direc
        double move_dist
        int cum_dist

        bint frozen

        MinionDLLNode node

        double max_energy
        double energy
        double avg_consum_rate
        double basal_metabolic_rate
        double move_consum_rate
        double stretch_consum_rate
        double sex_consum_rate
        double birth_consum_rate

        bint increase_age(self)
        double energy_with_constant(self,double const)
        void take_energy(self,double amount)
        bint loss_energy(self,double amount)
        void adjust_energy_from_mass(self)
        void take_mass(self,double amount)
        void loss_mass(self,double amount)
        bint mergeable(self,Minion mi)
        Minion get_child(self,Minion mi)
        double[:] get_input(self,int[:,:] snapshot)


cdef class MinionDLLNode():
    cdef:
        Minion mi
        MinionDLLNode prev
        MinionDLLNode next
        bint is_tail
        MinionDLL dll
cdef class MinionDLL():
    cdef:
        int len
        MinionDLLNode head
        MinionDLLNode tail
        MinionDLLNode current

        MinionDLLNode prev(self,MinionDLLNode node)
        MinionDLLNode next(self,MinionDLLNode node)
        MinionDLLNode push(self,Minion mi)
        MinionDLLNode remove_and_get_next(self,MinionDLLNode node)
        void remove_by_link(self,Minion mi)
        void remove_by_search(self,Minion mi)
        bint contains(self,Minion mi)



cdef class World:
    cdef:
        int[:,:] snapshot
        int xsize
        int ysize
        int moment
        int new_id
        double[:,:] mins
        MinionDLL mis
        MinionDLL[:,:] occupy_map
    cdef:
        bint no_age
        bint no_birth
        bint no_eat
        bint no_energy
        bint no_excrete
        bint no_hunt
        int messiness
        bint halluc
        double hidden_mass
    cdef:
        bint record_pedigree
    cdef: #registering
        void register(self, Minion mi)
        void unregister_deads(self)
    cdef:
        void take_mass(self,Minion mi,double amount)
        void loss_mass(self,Minion mi,double amount)
        void take_energy(self,Minion mi,double amount)
        bint loss_energy(self,Minion mi,double amount)
    cdef:
        void childbirth(self,Minion girl,Minion boy)
        bint huntable(self,Minion pred,Minion pray)
        void mk_corpse(self,Minion mi)
        void kill(self,Minion mi,bint corpse)
        void kill_elderly(self)
        void hunt(self,Minion pred,Minion pray)
    cdef:
        void move(self,Minion mi)
        void stretch(self,Minion mi)
        void excrete(self,Minion mi,double amount)
        void digest(self,Minion mi,double amount)
        void try_hunt(self,Minion mi)
        void eat(self,Minion mi)
        void act(self,Minion mi)
        void basal_metabolism(self,Minion mi)
    cdef:
        void control_all(self)
        void act_all(self)
        void basal_metabolism_all(self)
        void render(self)
    cdef:
        double _total_mass(self)
        double _total_minion_mass(self)
        double _total_min(self)
        bint _exhausted(self)
        double _get_avg_r(self)
        double _get_avg_g(self)
        double _get_avg_b(self)
    """
    cpdef:
        void one_step(self)
    cpdef:
        size_t get_xsize(self)
        size_t get_ysize(self)
        double total_mass(self)
        double total_min(self)
    cpdef:
        bint limit(self,int n)
        bint exhausted(self)
    cpdef:
        vector[(int,int)] get_pedigree(self)
        np.ndarray sample_dna(self,int n)
        int get_population(self)
        int get_avg_r(self)
        int get_avg_g(self)
        int get_avg_b(self)
    cpdef:
        np.ndarray body_rect(self,int k)
        np.ndarray geni_rect(self,int k)
    
    """
