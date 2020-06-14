#!/bin/bash
f2py3 --f90flags=-ffree-line-length-512 -c sip.f90 -m sip
 
#gfortran -ffree-line-length-512 *.f90