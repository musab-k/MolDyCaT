#!/bin/bash
#NB: must be run from moldycat directory 
cp analysis/elasticConstants/c44script.gnu temp.gnu
sed -e "s:DUMMY:${1}:" temp.gnu >> temp2.gnu
gnuplot temp2.gnu
rm temp*
