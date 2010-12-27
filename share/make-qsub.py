#!/usr/bin/python
import os
import string
for run in range(20333,20368):
  filename = str(run)+".sh"
  f = open(filename, 'w')
  f.write("cd  $PSIP_BATCH\n")
  f.write("boss.exe boss-options/psip-"+str(run)+".txt\n")
