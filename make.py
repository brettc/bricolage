"""
Attempt to use fabricate to build stuff. Unfortunately, it seems to be
broken no the Mac.
"""
import os
from fabricate import *
# import string

environment = dict(PYTHON='/Users/brett/anaconda')

def optionize(text):
    # template = string.Template(text)
    # sub_text = template.substitute(environment)
    sub_text = text.format(**environment)
    return [_ for _ in map(lambda x: x.strip(), sub_text.split('\n')) if _]

cflags_text = """
    -fno-strict-aliasing 
    -I{PYTHON}/include 
    -arch 
    x86_64
    -DNDEBUG 
    -g 
    -fwrapv 
    -O3 
    -Wall 
    -Wstrict-prototypes
    -I{PYTHON}/lib/python2.7/site-packages/numpy/core/include
    -I{PYTHON}/include/python2.7
    -Wno-c++11-extensions
    -Wno-unused-function
"""

ldflags_text = """
    -bundle
    -undefined
    dynamic_lookup
    -L{PYTHON}/lib
    -L/Users/brett/Dropbox/Code/organismal/organismal
    -arch 
    x86_64
"""

cflags = optionize(cflags_text)
ldflags = optionize(ldflags_text)

def build():
    cythonize('organismal/pubsub2_ext')
    compile('organismal/pubsub2_c')
    compile('organismal/pubsub2_ext')
    # link(objects=['organismal/pubsub2_ext.o', 'organismal/pubsub2_c.o'], build_dir='organismal',
    #      target='organismal/pubsub2_ext.so')

    # link(objects=['organismal/pubsub2_c.o'], build_dir='organismal', target='organismal/pubsub2_c.so')
    # link(objects=['organismal/pubsub2_ext.o'], build_dir='organismal',
    #      target='organismal/pubsub2_ext.so', flags=['-lpubsub2_c'])
    # link(build_dir='organismal', target='organismal/pubsub2_ext.so')

def oname(build_dir, filename):
    return os.path.join(build_dir, os.path.basename(filename))

def compile(source, flags=None):
    run('gcc', '-c', source+'.cpp', '-o', source+'.o', cflags, flags)

def cythonize(src): 
    run('cython', '--cplus', src+'.pyx', '-o', src+'.cpp')

def link(objects, build_dir, target, flags=None):
    # objects = 'organismal/pubsub2_c.o organismal/pubsub2_ext.o'.split()
    run('gcc', objects, '-o', oname(build_dir, target), ldflags, flags)


def clean():
    autoclean()

main()
