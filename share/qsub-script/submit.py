#!/usr/bin/python
import os
import string
for run in range(20333,20368):
  com = "qsub qsub-script/"+str(run)+".sh"
  print com
  os.system(com)
