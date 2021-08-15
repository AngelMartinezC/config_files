#!/bin/bash

echo " " > variable

while read p;
do
  echo $p | awk '{print $4, $6}' | thermodynam | grep "density" | awk '{print $3}' >> variable0
  echo $p | awk '{print $4, $6}' | thermodynam | grep "Gamma1" | awk '{print $3}' >> variable1
  echo $p | awk '{print $4, $6}' | thermodynam | grep "gamma" | awk '{print $3}' >> variable2
  echo $p | awk '{print $4, $6}' | thermodynam | grep "sound" | awk '{print $3}' >> variable3
done < cptrho.l5bi.d.15c

paste variable0 variable1 variable2 variable3 > variable

