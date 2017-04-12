# README #

This is a C and Python libary for calculation of the two-point correlation function on 3D cartesian data.

### Installation instructions ###

```
#!bash

cd ~/
git clone https://apellegrino@bitbucket.org/apellegrino/clustering-tree.git
cd clustering-tree
make python
```
Execute this next line in the same directory, or...
```
#!bash
export PYTHONPATH="$PYTHONPATH:$(pwd)"
```
...add this line to the end of your ~/.bashrc file:
```
#!bash
export PYTHONPATH="$PYTHONPATH:/home/[your username here]/clustering-tree"
```
Then finally:
```
#!bash
source ~/.bashrc
```
You can now use the library in python with
```
#!python
import tpcf
```