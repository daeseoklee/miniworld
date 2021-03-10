from distutils.core import setup
from Cython.Build import cythonize
setup(\
ext_modules=cythonize([\
'./fast_random.pyx'\
],annotate=True),)
setup(\
ext_modules=cythonize([\
'./fast_random_py.pyx'\
],annotate=True),)
setup(\
ext_modules=cythonize([\
'./world_cy.pyx'\
],annotate=True),)
