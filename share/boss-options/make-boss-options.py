#!/usr/bin/python
import os
import string
import sys
import fileinput
import re


import os

def configure(source_file_name,  target_file_name,  TEMPLATE_RUN_NUMBER, TEMPLATE_DST_FILES):
#Open source template files with simulation configuration
	source_file = open(source_file_name, 'r')
	target_file = open(target_file_name, 'w')
	for line in source_file:
		line = re.sub("TEMPLATE_DST_FILES", TEMPLATE_DST_FILES, line)
		line = re.sub("TEMPLATE_RUN_NUMBER",TEMPLATE_RUN_NUMBER, line)
		target_file.write(line)
	source_file.close()
	target_file.close()

def proceed(run, directory, files):
    TEMPLATE_RUN_NUMBER=str(run)
    TARGET_FILE = TEMPLATE_RUN_NUMBER+".cfg"
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
    if TEMPLATE_DST_FILES=='': return
    print TEMPLATE_DST_FILES
    configure('template.txt',TARGET_FILE,TEMPLATE_RUN_NUMBER, TEMPLATE_DST_FILES)

for run in range(25244,25337):
  print "Proceeding run ", run
  os.path.walk("../data", proceed, run)
  #create qsub files
  filename = str(run)+".sh"
  f = open(filename, 'w')
  f.write("cd  $PSIP_BATCH\n")
  f.write("boss.exe boss-options/"+str(run)+".cfg\n")

