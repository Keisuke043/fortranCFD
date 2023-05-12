#!/bin/sh
rm chem.cti
ck2cti --input=chem.inp --thermo=therm.dat --transport=tran.dat
