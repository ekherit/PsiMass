#!/usr/bin/python
import os
import string
import sys
import fileinput
import re


import os

def callback(run, directory, files):
    r = "run_0+"+str(run)+"_.+\.dst"
    flist = []
    TEMPLATE_DST_FILES=''
    for f in files:
      if re.match(r,f):
        name = os.path.join(directory, f);
        flist.append(name)
        if TEMPLATE_DST_FILES=='': comma=''
        else: comma=',\n'
        TEMPLATE_DST_FILES=TEMPLATE_DST_FILES+comma+'"'+name+'"'
    print TEMPLATE_DST_FILES

run = 25201
print "Proceeding run ", run
os.path.walk("/home/nikolaev/work/BES/test", callback, run)

#for run in range(0025201,0025202):
#  filename = str(run)+".sh"
#  f = open(filename, 'w')
#  f.write("cd  $PSIP_BATCH\n")
#  f.write("boss.exe boss-options/psip-"+str(run)+".txt\n")
