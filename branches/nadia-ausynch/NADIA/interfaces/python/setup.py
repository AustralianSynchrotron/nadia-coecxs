import sys, os,stat,commands
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from Cython.Build.Dependencies import create_extension_list
# 
libNADIAIncludeDir=os.path.join('..','..','include')
libNADIALibDir=os.path.join('..','..','lib')


DoublePrecision='1'

os.environ["LD_RUN_PATH"]=os.path.join('..','..','lib')

if DoublePrecision == '1':
    nadialib = 'NADIAd'
else:
    nadialib = 'NADIA'


modules = create_extension_list([os.path.join('pyNADIA','*.pyx')])

for module in modules:
    module.language="c++"
    module.include_dirs=[".",libNADIAIncludeDir ]
    module.libraries=['tiff', 'm', 'df', 'z','mfhdf', 'stdc++', 'jpeg',nadialib]
    module.library_dirs=[libNADIALibDir]
    module.runtime_library_dirs=[libNADIALibDir]
    module.package="pyNADIA"
    module.define_macros=[('DOUBLE_PRECISION',DoublePrecision)]
    module.gdb_debug=True
    print module.name

setup(
  name="pyNADIA",
  description='Python Bindings for NADIA',
  author='Lenneke Jong',
  author_email='Lenneke.Jong@synchrotron.org.au',
  #packages=["pyNADIA"],
  packages=["pyNADIA"],
  ext_modules=cythonize(modules, compile_time_env={'DOUBLE_PRECISION':DoublePrecision}),
  #ext_modules=cythonize("pyNADIA.pyx"),
  cmdclass = {'build_ext': build_ext},
)
