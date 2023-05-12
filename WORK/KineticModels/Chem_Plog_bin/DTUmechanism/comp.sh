#!/bin/sh
./clean.sh
gfortran -o CKINTERP ckinterp.f -fallow-argument-mismatch
gfortran -o TRAN cklib.f xerror.f tranfit.f -fallow-argument-mismatch
chmod 755 CKINTERP
chmod 755 TRAN
./CKINTERP
./TRAN > tran.out

