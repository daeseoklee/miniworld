
cdef double uniform()
cdef int randint(int n,int m)
cdef int randint_with_probs(int n,int m,double[:] probs)
cdef bint randbool()
cdef double[:] multiple_uniform(int n)
cdef double gaussian()
cdef double[:] multiple_gaussian(int n)
cdef void seed()
cdef void fix_seed()
cdef void assign_random_gaussian_pair(double[:] out, int assign_ix)