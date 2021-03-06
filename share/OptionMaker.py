#Module for configuring boss options
import os
import sys
import re
import random

#create all file list (used for path.walk
def proceed_create_file_list(filelist, directory, files):
    for file in files:
        filelist += [os.path.join(directory,file)]
    filelist.sort()

#create all file list in directory
def create_file_list(directory):
    files = []
    os.path.walk(directory, proceed_create_file_list, files)
    return files;

def filter_file_list(files, reg):
    r = re.compile(reg)
    filtered_file_list = []
    for file in files:
        if re.match(r,file):
            filtered_file_list += [file]
    return filtered_file_list

#for every run create file list
def create_run_dict(files,reg):
    RunMap = {}
    r = re.compile(reg)
    for file in files:
        m = re.match(r,file)
        if m:
            #g = int(m.group(1))
            #run = "%07d" % int(m.group(1))
            run = int(m.group(1))
            #print run
            if run in RunMap:
                RunMap[run]+=[file]
            else:
                RunMap[run]=[file]
    return RunMap


class RunInfo:
    run  = 0
    files = []
    def __init__(self, Run, Files):
        self.run = Run
        self.files = Files

def group_files(run_map, n) :
    run_info_list = []
    count=0
    tmp_list=[]
    for run, files in run_map.items():
        ri = RunInfo(run,files)
        if count < n :
            tmp_list.append(ri) 
        else :
            run_info_list.append(tmp_list)
            tmp_list=[]
            count=0
            tmp_list.append(ri)
        count=count+1
    return run_info_list

#create comma separated and qouted list of files
def make_files_string(files):
    files_string=''
    comma=''
    for file in files:
        files_string=files_string+comma+'"'+file+'"'
        comma=',\n  '
    return files_string


def lookup( f, d):
    if not os.path.exists(f):
        f = os.path.join(d,f)
    if os.path.exists(f):
        f = os.path.realpath(os.path.abspath(f))
        #print f
    else :
        print 'File ',  f,  ' doest exists'
        sys.exit(1)
    return f



class OptionMaker:
    JPSIKKROOT_DIR       = os.path.abspath(os.environ['JPSIROOT'])
    TEMPLATE_DIR         = os.path.join(JPSIKKROOT_DIR, 'share/template')
    SHARE_DIR            = os.path.join(JPSIKKROOT_DIR, 'share')
    dataDir = ""   #folder where data files will be searched
    templateFile=""# the path to the template file
    decayFile=""   #the decay file
    targetDir="."   #the folder where options file will be located
    jobPrefix=""   #prefix of the job
    fileFilter=".*run_(\d\d\d\d\d\d\d).*.dst"  #dst regex filter
    runFilter=""   #regexp run grouping
    eventNumber=-1 #Event number per job
    runNumber=1    #Run number per job
    runs=""
    jobNumber=1    #number of jobs
    Seed=1         #the random seed
    fileList=[]
    runMap={}
    SelectionMode = True
    SimulationMode = False
    ReconstructionMode = False
    RunInfoList = []
    Energy = 0.0


    def __init__(self, options, args):
        self.decayFile = args[1]
        self.dataDir   = os.path.realpath(os.path.abspath(args[1]))
        self.eventNumber=options.event_number
        self.jobNumber=int(options.job_number)
        self.jobPrefix=options.job_prefix
        self.runs=options.runs
        self.runNumber = options.run_number
        if options.energy != "":
            self.Energy = float(options.energy)

        if len(args) < 1:
            print 'Specify the action: "sel",  "sim",  "rec"'
            sys.exit(1)

        if(args[0] == "selection" or args[0]=="sel"):
            self.SelectionMode=True
            #self.fileFilter=".*run_(\d\d\d\d\d\d\d).*.dst"
            self.fileFilter = ".*(\d{4,7}).*.dst"
            self.fileList = filter_file_list(create_file_list(self.dataDir), self.fileFilter)
            self.templateFile = "selection.cfg"
            if options.run_group == "":
                self.group(".*[^\d](\d{4,7})[^\d].*.dst", self.runNumber)
            else:
                print "Using run_group option: ", options.run_group
                self.group(options.run_group, self.runNumber)

        if args[0] == "simulation" or args[0] == "sim":
            self.SimulationMode = True
            self.templateFile = "simulation.cfg"
            self.decayFile = lookup(self.decayFile,  self.SHARE_DIR)
            print "Making simulation config files."

        if args[0] == "reconstruction" or args[0] == "rec":
            self.ReconstructionMode = True
            self.templateFile = "reconstruction.cfg"
            self.fileFilter = ".*(\d\d\d\d\d\d\d).*.rtraw"
            self.fileList =filter_file_list(create_file_list(self.dataDir),  self.fileFilter)
            #print self.fileList
            print "Making reconstruction config files."
            self.group(".*(\d{5,7}).rtraw",self.runNumber)

        if(args[0] == "pbs"):
            self.PbsMode=True
            #self.fileFilter=".*run_(\d\d\d\d\d\d\d).*.dst"
            self.fileFilter = ".*(\d{4,7}).*.cfg"
            self.fileList = filter_file_list(create_file_list(self.dataDir), self.fileFilter)
            self.templateFile = "pbs.sh"
            self.group(".*[^\d](\d{4,7})[^\d].*.cfg", self.runNumber)

        self.templateFile = lookup(self.templateFile,  self.TEMPLATE_DIR)
        print "Setup template file: ",  self.templateFile

        #if not os.path.exists(self.templateFile):
        #    self.templateFile = os.path.join(self.TEMPLATE_DIR, self.templateFile)
        #    print self.templateFile
        #else:
        #    self.templateFile = os.path.realpath(os.path.abspath(self.templateFile))
        #    print self.templateFile
        #if not os.path.exists(self.templateFile):
        #    print 'Template file ',  options.template_file,  ' doest exists'
        #    sys.exit(1)
	if len(args) >= 2:
		self.targetDir = os.path.realpath(os.path.abspath(args[2]))
        if not os.path.exists(self.targetDir):
            os.mkdir(self.targetDir)
        else:
            if not options.force:
                print "Target dir ",  self.targetDir,  " already exists, exiting"
                sys.exit(1)



    def group(self, run_filter=".*run.*[^\d](\d{4,7})[^\d]+.*.dst", run_number=1):
        self.runFilter=run_filter
        self.runNumber = run_number
        self.runMap = create_run_dict(self.fileList, run_filter)
        self.RunInfoList = group_files(self.runMap,run_number);
        #self.runMap.sort()

    def setup_template_file(self, template_file, target_dir):
        self.templateFile  = template_file
        if not os.path.exists(self.templateFile):
            print "Template file ", self.templateFile, " does not exists"
        self.targetDir  = target_dir
        if not os.path.exists(self.targetDir):
            os.mkdir(self.targetDir)
        else:
            print "Target dir ",  self.targetDir,  " already exists, exiting"
            sys.exit(1)

    def make(self):
        if self.SimulationMode:
            self.make_sim()
            return

        if self.ReconstructionMode:
            self.make_rec()
            return

        if self.SelectionMode:
            self.make_sel()
            return


    def make_sel(self):
        #print self.runMap
        self.run_string()
        for run, files in self.runMap.items():
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.root" % (self.jobPrefix,run)))
            TemplateInputFile = make_files_string(files)
            #define the name of cfg file
            target_file_name = os.path.join(self.targetDir,"%s-%07d.cfg" % (self.jobPrefix, run))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",str(self.eventNumber), line)
                line = re.sub("TEMPLATE_RUN_NUMBER",str(self.runNumber), line)
                line = re.sub("TEMPLATE_CENTER_MASS_ENERGY",str(self.Energy), line)
                target_file.write(line)
            source_file.close()

    def make_sel2(self):
        self.run_string()
        #runinfo_list = group_files(self.runMap,2)
        for rilst in self.RunInfoList:
            runs=[]
            files = []
            for ri in rilst:
                runs.append(ri.run)
                for f in ri.files:
                    files.append(f)
            runs.sort()
            runs_string=""
            for run in runs:
                runs_string+=("_%07d" % run )
            TemplateInputFile = make_files_string(files)
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s%s.root" % (self.jobPrefix,runs_string)))
            target_file_name = os.path.join(self.targetDir,"%s-sel%s.cfg" % (self.jobPrefix, runs_string))
            print "Creating ", target_file_name, "  for runs: ", runs
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",str(self.eventNumber), line)
                target_file.write(line)
            source_file.close()
            target_file.close()

    def make_pbs(self):
				return


    def make_sim(self):
        for job in range(0, self.jobNumber):
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.rtraw" % (self.jobPrefix,job)))
            TemplateInputFile =  self.decayFile 
            #define the name of cfg file
            target_file_name = os.path.join(self.targetDir,"%s-sim-%07d.cfg" % (self.jobPrefix, job))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE", TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE",TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED",TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",str(self.eventNumber), line)
                line = re.sub("TEMPLATE_RUN_NUMBER",self.runs, line)
                target_file.write(line)
            source_file.close()
            target_file.close()

    def make_rec(self):
        for run, files in self.runMap.items():
            #define input and output files in joboptions
            TemplateOutputFile = os.path.abspath(os.path.join(self.targetDir,"%s-%07d.dst" % (self.jobPrefix,run)))
            TemplateInputFile = make_files_string(files)
            #define the name of cfg file
            #print run
            target_file_name = os.path.join(self.targetDir,"%s-rec-%07d.cfg" % (self.jobPrefix, run))
            source_file = open(self.templateFile, 'r')
            target_file = open(target_file_name, 'w')
            TemplateRandomSeed = str(random.randint(0,2**32))
            for line in source_file:
                line = re.sub("TEMPLATE_INPUT_FILE",  TemplateInputFile, line)
                line = re.sub("TEMPLATE_OUTPUT_FILE", TemplateOutputFile, line)
                line = re.sub("TEMPLATE_RANDOM_SEED", TemplateRandomSeed, line)
                line = re.sub("TEMPLATE_EVENT_NUMBER",str(self.eventNumber), line)
                line = re.sub("TEMPLATE_RUN_NUMBER",  str(self.runNumber), line)
                target_file.write(line)
            source_file.close()
            target_file.close()

    def run_string(self):
        first_run=0
        last_run=0
        range_list = []
        runs = self.runMap.keys();
        runs.sort()
        for run in runs:
            if first_run==0 : #initialization
                first_run=run
                last_run=run
                continue
            if run-last_run==1: #expand continouse range
                last_run = run
                continue

            if first_run < last_run: #fill range
                range_list.append(first_run)
                range_list.append(0)
                range_list.append(last_run)
            else: #fill single run
                range_list.append(first_run)
            #reset range
            first_run=run
            last_run=run
        #for r in range_list:
        #    print r


