import os
import sys

cur_file_dir = os.path.dirname(os.path.realpath(__file__))

# add current file directory so that mypackage_bindings.so is found by python
sys.path.append(cur_file_dir)

# set current file directory as working dir so that mypackage_bindings.so dependancies
# will be found by the linker (mypackage_bindings.so and its deps RPATH are set to $ORIGIN)
os.chdir(cur_file_dir)

# load every symbols of mypackage_bindings into upper mypackage module
from highspy_bindings import *