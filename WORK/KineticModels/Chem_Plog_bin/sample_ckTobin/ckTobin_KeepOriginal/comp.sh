#!/bin/sh
rm CKINTERP TRAN chem.bin tran.bin chem.out tran.out
gfortran -o CKINTERP ckinterp.f
gfortran -o TRAN     cklib.f xerror.f tranfit.f
./CKINTERP
./TRAN > tran.out
