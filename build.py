#!/usr/bin/python

import sys
from pprint import pprint as pp
import pdb
import os
import argparse
import subprocess
import time
from datetime import date

sttime = time.time()

parser = argparse.ArgumentParser(description='CFD/Lite build script')
#parser.add_argument('-r', '--rock', default=False, required=False)

#args = parser.parse_args()

build_map = {
    'main':'src/main.f90'
}

dep_map = {}
name_map = {} # short name --> full path
name_map_inv = {} # full path --> short name

fcompiler = 'gfortran'#'mpif90 -f90=pgf90 '
ccompiler = 'gcc'
#fcompiler = 'pgfortran'
#ccompiler = 'pgcc'
#fcompiler ='/fast-nas/apps/pgi-17.4/linux86-64/17.4/bin/pgfortran'
#ccompiler= '/fast-nas/apps/pgi-17.4/linux86-64/17.4/bin/pgcc'

cflags = ' -c -g  -fdefault-real-8 -ffree-line-length-512 -cpp' # gnu
#cflags = ' -c -g -r8 '# pgi

cwd = os.getcwd()
#preprocessDir = cwd + '/src/util'

hdf5='/opt/hdf5'
cgns='/opt/cgnslib_3.2.1_modified'
fftw='/opt/fftw-3.3.4'
sqlite='/opt/sqlite'
cuda='8.0'
cc='cc35'
server='Rock'

xml2 = '/opt/libxml2-2.9.2'

build_cmd = 'cd ..'
os.system(build_cmd)

links= '-ldl -L '+hdf5+'/lib -lhdf5 -L '+cgns+'/lib -lcgns -L '+fftw+'/lib -lfftw3f'# -L'+sqlite+'/lib -lsqlite3 -L '+
# default line length limit 132. this can be increased using -ffree-line-length-<n>
includes = '-I '+cgns+'/include'# -I '+fftw+'/include -I '+sqlite+'/include -I '+ preprocessDir +' -I '+partitioner+'/include '

#links = ' -Bstatic '+links
build_path = os.path.realpath('build')
bin_path = os.path.realpath('bin')

def mkpath(p):
    return os.path.sep.join(p)

def mkcmd(cmd):
    return ' '.join(cmd)

def get_name(t):
    return os.path.splitext(os.path.basename(t))[0]

class SourceFile(object):
    def __init__(self, src_file, dependencies, src_map):
        self.full_name = src_file
        self.full_path = os.path.realpath(src_file)
        self.short_name = get_name(self.full_name)
        self.deps = dependencies
        self.parents = []
        self.src_map = src_map
        self.update = False
        self.compiled = False

        self.optm = ' -O0'
        #print 'No optomization: '+self.full_name

    def up_to_date(self):
        fileobj = self.short_name + '.o'
        if os.path.exists(fileobj):
            if os.path.getmtime(fileobj) >= os.path.getmtime(self.full_path):
                return True
            else:
                return False
        else:
            return False

    def compile(self):
        # verify dependencies have been compiled
        for dep in self.deps:
            self.src_map[dep].compile()
        if not self.compiled:
            if self.update or any([self.src_map[dep].compiled for dep in self.deps]):
                macroFlag = ''#V15+'-DSRC_FILE=\'\"' + self.full_path + '\"\''
                # compile this SourceFile
                compile_cmd = mkcmd([fcompiler, cflags + self.optm + macroFlag, includes, self.full_path])
                print "\t[" + os.path.split(fcompiler)[1] +self.optm+ cflags +"]: " + self.full_path
                os.system(compile_cmd)

                self.compiled = os.path.exists(self.short_name + '.o')
                if not self.compiled:
                  sys.exit("Compilation of '"+self.short_name+"' failed!")
                self.update = False

    def check(self):
        if not self.up_to_date():
            self.update = True

# read list of source files
with open('.srcfiles', 'r') as f:
   src_files = [x.strip() for x in f.readlines()]

# read dependency data as computed by mkdep
with open('.deps', 'r') as f:
    dep_lines = f.readlines()[4:]

# construct dependency map
for line in dep_lines:
    l = line.strip().split(':', 1)
    src_name = get_name(l[0].strip())
    dep_names = [get_name(x.strip()) for x in l[1].strip().split(' ')]
    dep_map[src_name] = dep_names
    #print src_name, '-->', repr(dep_names)

# construct source name maps
src_map = {}
for src_file in src_files:
    src_name = get_name(src_file)
    name_map[src_name] = src_file
    name_map_inv[src_file] = src_name

# construct dependency tree of SourceFile objects
for src_file in src_files:
    src_name = name_map_inv[src_file]
    if src_name in dep_map.keys():
        deps = [name_map[dep_name] for dep_name in dep_map[src_name]]
    else:
        deps = []
    src_map[src_file] = SourceFile(src_file, deps, src_map)

# compute parents for SourceFile objects
for src_file in src_map.keys():
    for dep in src_map[src_file].deps:
        src_map[dep].parents.append(src_file)

def print_tree(src_file, src_map, depth=0):
    s = src_map[src_file]
    print '  '*depth + s.short_name
    for dep in s.deps:
        print_tree(dep, src_map, depth+1)

os.chdir(build_path)

# check for files which have been changed since last compile
for src_file in src_map.keys():
    src_map[src_file].check()
    if src_map[src_file].update:
        print src_file

# print "\n[Compile]  Flags:\t" + fcompiler + " " + cflags + optm + " " + includes + "\n"
for src_file in src_map.keys():
    src_map[src_file].compile()

# compile sqlite bindings for fortran (linked with c code)
#build_cmd = mkcmd([ccompiler, '-c', '-DSGL_UNDERSCORE -I '+sqlite+'/include', '../src/modules/csqlite.c'])
#print build_cmd
#os.system(build_cmd)

link_cmd = mkcmd([fcompiler, '*.o', links, ' -lstdc++ -o %s/runcfd' % bin_path])
print link_cmd
os.system(link_cmd)

endtime = time.time()
tt = endtime-sttime
print "After %.1f seconds..." %(tt)
