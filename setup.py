from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("minion.py")
)
setup(
    ext_modules=cythonize("world.py")
)
setup(
    ext_modules=cythonize("play.py")
)
"""
setup(
    ext_modules=cythonize("ex.py")
)
"""