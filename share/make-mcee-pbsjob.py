#!/usr/bin/python
import os
import string
import sys
import fileinput
import re

ENERGY=1838
RANDOM_SEED=0
EVENT_NUMBER=1
JOB_NAME="test"
#PBS_QUEUE="publicq@torqsrv"
PBS_QUEUE="besq@torqsrv"

def configure(source_file_name,  target_file_name,  JOB_NAME,  ENERGY,  RANDOM_SEED,  EVENT_NUMBER,  RUN_NUMBER):
#Open source template files with simulation configuration
	source_file = open(source_file_name, 'r')
	target_file = open(target_file_name, 'w')
	for line in source_file:
		line = re.sub("TEMPLATE_NAME", JOB_NAME, line)
		line = re.sub("TEMPLATE_BEEM_ENERGY", str(ENERGY/1e3), line)
		line = re.sub("TEMPLATE_RANDOM_SEED", str(RANDOM_SEED), line)
		line = re.sub("TEMPLATE_EVENT_NUMBER", str(EVENT_NUMBER), line)
		line = re.sub("TEMPLATE_PBS_QUEUE", PBS_QUEUE, line)
		line = re.sub("TEMPLATE_RUN_NUMBER", RUN_NUMBER, line)
		target_file.write(line)
	source_file.close()
	target_file.close()

def do_mc(energy,  event_number,  job_number,  run_number):
	random_seed=int(energy*100)/10+job_number
	name="bbyg_ee_"+str(energy)+"_"+str(job_number)
	template_dir="/bes3fs/groups/tauqcd/tauqcdgroup/nikolaev/mc/mcee-template"
	work_dir = name
	if os.path.exists(work_dir):
		print "Directory "+work_dir+" exist! Do nothing!!"
		return
	else:
		print "Creating directory "+work_dir
		os.mkdir(work_dir)
	JOB_NAME=name
	RANDOM_SEED=random_seed
	EVENT_NUMBER=event_number
	ENERGY=energy
	RUN_NUMBER=run_number

#Do substitution for required value
	configure(template_dir+"/template_sim.cfg", work_dir+"/"+name+"_sim.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/template_rec.cfg", work_dir+"/"+name+"_rec.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/template_ana.cfg", work_dir+"/"+name+"_ana.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/common_sim.cfg", work_dir+"/common_sim.cfg", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)
	configure(template_dir+"/pbsjobs", work_dir+"/pbsjobs", JOB_NAME, ENERGY, RANDOM_SEED, EVENT_NUMBER,  RUN_NUMBER)

	print "Starting the "+PBS_QUEUE+" job " + name
	os.system("qsub "+work_dir+"/pbsjobs")
#end of do_mc function

N=50000
jobs = 4
for job in range(1,jobs+1):
	do_mc(1839.0, N, job,  "-20334, -20335, -20339, -20336")
	do_mc(1841.4, N, job,  "-20340, -20341")
	do_mc(1843.2, N, job,  "-20344, -20342, -20343")
	do_mc(1844.1, N, job,  "-20346, -20347")
	do_mc(1848.4, N, job,  "-20350, -20348, -20349")
	do_mc(1838.2, N, job,  "-20353, -20354")
	do_mc(1840.8, N, job,  "-20357, -20355")
	do_mc(1841.8, N, job,  "-20358, -20359")
	do_mc(1843.0, N, job,  "-20361, -20360")
	do_mc(1843.5, N, job,  "-20362, -20363")
	do_mc(1845.1, N, job,  "-20365, -20364")
	do_mc(1846.5, N, job,  "-20366, -20367")
	#this is extra point
	do_mc(1842.0, N, job,  "-20361, -20360")
	do_mc(1842.5, N, job,  "-20361, -20360")
	do_mc(1842.75, N, job,  "-20361, -20360")
