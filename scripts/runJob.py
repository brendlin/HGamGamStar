#!/usr/bin/env python

import ROOT
import os,sys
import math
import datetime

import HelperTools

#-------------------------------------------------------------------------
def Error(message) :
    print 'Error: %s. Exiting.'%(message)
    sys.exit()

#-------------------------------------------------------------------------
def getOutputName(config) :

    if not config.isDefined('OutputDir') :
        ctime = datetime.datetime.now()
        thetime = ctime.strftime('%Y.%m.%d_%H.%M.%S')
        config.setValue('OutputDir','%s_%s'%(ROOT.TString(options.Alg),ROOT.TString(thetime)))

    submitDir = config.getStr('OutputDir').Data()

    if os.path.exists(submitDir) :
        print 'Output directory %s already exists.'%(submitDir)
        print 'Rerun after deleting it, or specify a new one, like below.'
        print 'OutputDir: NEWDIR'
        sys.exit()

    return submitDir

#-------------------------------------------------------------------------
def getFilesFromCommandLine(Input,InputList) :
    files = []

    if Input :
        files += Input.split(',')

    if InputList :
        for fileName in open(InputList).readlines() :
            fileName = fileName.replace('\n','')
            files.append(fileName)

    return files

#-------------------------------------------------------------------------
def PrintSampleSummary(sh) :

    for sample in sh :

        if issubclass(type(sample),ROOT.SH.SampleGrid) :
            # For SampleGrid, numFiles is not supported.
            inds = sample.meta().castString(ROOT.SH.MetaFields.gridName)
            outds = sample.meta().castString('nc_outputSampleName')
            print 'Sample \"%s\": \n - INDS: %s \n - OUTDS: %s'%(sample.name(),inds,outds)
            continue

        nentries = sample.meta().castDouble(ROOT.SH.MetaFields.numEvents,-1)
        if nentries > 0 :
            print 'Sample \"%s\" has %d entries'%(sample.name(),nentries)
        else :
            print 'Sample \"%s\":'%(sample.name())

        # print first 5 files
        for i in range(min(5,sample.numFiles())) :
            print '    %s'%(sample.fileName(i))

        # print last 3 files
        if sample.numFiles() > 5 :
            print '    ...'
            for i in range(max(5,sample.numFiles()-3),sample.numFiles()) :
                print '    %s'%(sample.fileName(i))

    return

#-------------------------------------------------------------------------
def SetEventsPerWorker(sh,conf) :
    # for Condor or GridEngine

    print 'scanning number of events...'
    ROOT.SH.scanNEvents(sh)

    maxEvents = conf.getNum('nc_EventLoop_EventsPerWorker')
    sh.setMetaDouble(ROOT.EL.Job.optEventsPerWorker,maxEvents)

    # Print out information about the samples
    for sample in sh :
        nEventsInSample = sample.meta().castDouble(ROOT.SH.MetaFields.numEvents, -1)
        nJobs = int( math.ceil(nEventsInSample / float(maxEvents) ) )
        details = 'Sample \"{}\" has {} entries (nc_EventLoop_EventsPerWorker={}, meaning {} jobs).'
        print details.format(sample.name(),nEventsInSample,maxEvents,nJobs)

    return

#-------------------------------------------------------------------------
def SetGridEngineOpts(driver,conf) :
    # Check your own system for what flags are necessary.
    # -V is needed to transfer the environment variables to the jobs
    # There is a default set below if not specified otherwise.

    sge_flags = conf.getStr('nc_EventLoop_SubmitFlags','')

    if not sge_flags :
        sge_flags = '-q default.q,short.q,long.q -l h_vmem=4G -V'

    driver.options().setString(ROOT.EL.Job.optSubmitFlags,sge_flags.Data())

    details = 'Submitting jobs to GridEngine batch system using options \"{}\"'
    details = '(set using nc_EventLoop_SubmitFlags)'

    print details.format(sge_flags)

    return

#-------------------------------------------------------------------------
def SetCondorOpts(driver,conf,options) :

    # Your system might need a specific 'shell init'
    driver.shellInit = ''

    if options.Condor_UseLD_LIBRARY_PATH :
        driver.shellInit += 'export LD_LIBRARY_PATH=%s;'%(os.getenv('LD_LIBRARY_PATH'))

    driver.shellInit += conf.getStr('Condor_shellInit').Data()

    optCondorConf = 'getEnv = True'

    # Possibility to add items to optCondorConf
    if (conf.getStr('optCondorConf', '').Length()) :
        optCondorConf += '\n'
        optCondorConf += conf.getStr('optCondorConf', '').Data()

    driver.options().setString(ROOT.EL.Job.optCondorConf, optCondorConf)

    details = 'Sumitting jobs to Condor batch system using options\n'
    details += 'shellInit: \"{}\" (set using --Condor_shellInit),\n'
    details +=  'optCondorConf: \"{}\" (set using --optCondorConf).'

    print details.format(driver.shellInit,optCondorConf)

    return

#-------------------------------------------------------------------------
def GetGridOptionsDouble() :
    return ['nc_nFiles','nc_nFilesPerJob','nc_nJobs']

def GetGridOptionsString() :
    return ['nc_site','nc_excludedSite','nc_EventLoop_SubmitFlags',
            'nc_mergeOutput','nc_official','nc_voms','nc_destSE']

#-------------------------------------------------------------------------
def SetGridOpts(driver,conf) :

    for opt in GetGridOptionsDouble() :
        if (conf.isDefined(opt)) :
            driver.options().setDouble(opt, conf.getInt(opt))

    for opt in GetGridOptionsString() :
        if (conf.isDefined(opt)) :
            driver.options().setString(opt, conf.getStr(opt).Data())

    if conf.getBool('UsexAODMerge', False) :

        if conf.isDefined('nc_EventLoop_SubmitFlags') :
            Error('You can NOT specify UsexAODMerge and nc_EventLoop_SubmitFlags at the same time')

        # TODO update this to use PathResolver - will also need an update to CMakeLists
        # to install scripts from scripts/ directory

        flag = '--mergeScript=__panda_rootCoreWorkDir/HGamAnalysisFramework/scripts/xaodmerge %OUT %IN'
        driver.options().setString(ROOT.EL.Job.optSubmitFlags,flag)

    return

#-------------------------------------------------------------------------
def main (options,args) :
    
    ####################################
    # Set up configuration class
    ####################################

    print 'Loading compiled code into python....'
    conf = ROOT.HG.Config()
    print 'Loading compiled code into python done.'

    # Add the main configuration file first
    conf.addFile(options.Config)

    # Configuration items from command-line override ones in the config file.
    for option,value in options.__dict__.items() :
        if value == None :
            continue
        conf.setValue(ROOT.TString(option),ROOT.TString(str(value)))

    # Get the submit directory name
    submitDir = getOutputName(conf)

    # Config uses value from first file it's specified in.
    # If specified, read in additional configuration.
    if (conf.isDefined('Include')) :
        for cfg in conf.getStrV('Include') :
            conf.addFile(cfg)

    # Require that a BaseConfig is present.
    if not conf.isDefined('BaseConfig') :
        msg = 'You must specify a base configuration file, using either '
        msg += '--BaseConfig Blah or within your config file'
        Error(msg)

    # Fill unspecified values from default config.
    conf.addFile(ROOT.PathResolverFindCalibFile(conf.getStr('BaseConfig').Data()))

    # Require nc_EventLoop_EventsPerWorker if using BatchGE BatchCondor
    if options.BatchGE or options.BatchCondor :
        if conf.getNum('nc_EventLoop_EventsPerWorker', -1) < 0 :
            Error('To run in batch mode, you must specify --nc_EventLoop_EventsPerWorker')

    # Grid running: Check for a GridTag
    if options.Grid and not conf.isDefined('GridTag') :
        msg = 'To submit to the grid, you MUST define a GridTag of the form: '
        msg += 'user.<UserName>.<UniqueTag>'
        Error(msg)

    # Add output directory name.
    conf.setValue('OutputDir',submitDir)


    ####################################
    # Set up job and algorithm
    ####################################

    myjob = ROOT.EL.Job()

    # Set to branch access mode: 'class' (default), 'branch', 'athena'
    myjob.options().setString(ROOT.EL.Job.optXaodAccessMode,
                              conf.getStr('xAODAccessMode','class').Data())

    # Execute 'alg = ROOT.MyAlgorithm("MyAlgorithm")'
    cmd = 'alg = ROOT.%s("%s")'%(options.Alg,options.Alg)
    print cmd
    exec (cmd)

    # Add the configuration to the algorithm
    alg.setConfig(conf)

    # Add the alg to the job
    myjob.algsAdd(alg)

    # Specify number of events
    if options.NumEvents > 0 :
        myjob.options().setDouble(ROOT.EL.Job.optMaxEvents,options.NumEvents)

    # Set number of events to be skipped
    if options.SkipEvents > 0 :
        myjob.options().setDouble(ROOT.EL.Job.optSkipEvents,options.SkipEvents)


    ####################################
    # Set up samples
    ####################################

    myhandler = ROOT.SH.SampleHandler()

    # Dataset stored on DESY-HH_LOCALGROUPDISK
    if options.GridDirect :

        if not os.getenv('RUCIO_HOME') :
            Error('Please execute localSetupRucioClients and voms-proxy-init')

        griddirect_from = os.getenv('GRIDDIRECT_FROM')
        griddirect_to = os.getenv('GRIDDIRECT_TO')
        localgroupdisk = os.getenv('LOCALGROUPDISK')

        if None in [localgroupdisk,griddirect_from,griddirect_to] :
            msg = 'To use GridDirect, you need to set GRIDDIRECT_FROM, '
            msg += 'GRIDDIRECT_TO, and LOCALGROUPDISK environment variables'
            Error(msg)
            
        griddsets = getFilesFromCommandLine(options.Input,options.InputList)

        for ds in griddsets :
            ROOT.SH.scanRucio(myhandler,ds)

        # last argument is whether to allow partial datasets:
        print 'SH::makeGridDirect GRIDDIRECT_FROM: \"%s\"'%(griddirect_from)
        print 'SH::makeGridDirect GRIDDIRECT_TO: \"%s\"'%(griddirect_to)
        print 'SH::makeGridDirect LOCALGROUPDISK: \"%s\"'%(localgroupdisk)
        ROOT.SH.makeGridDirect(myhandler,localgroupdisk,griddirect_from,griddirect_to,False)


    # Grid samples
    elif options.Grid :

        griddsets = getFilesFromCommandLine(options.Input,options.InputList)

        # # hack since a central cmtConfig is hardcoded
        # sh.setMetaString('nc_cmtConfig', os.getenv('AnalysisBase_PLATFORM'))

        # Set the file pattern, if specified in the config
        if conf.isDefined('nc_grid_filter') :
            myhandler.setMetaString('nc_grid_filter', conf.getStr('nc_grid_filter').Data())

        # Construct sample to run on
        for ds in griddsets :
            ROOT.SH.addGrid(myhandler, ds)

        HelperTools.SetOutputDatasetNames(myhandler,conf.getStr('GridTag').Data())

    # Local file(s), via Input or InputList (file or list of files)
    else :

        sampleName = conf.getStr('SampleName','sample')
        sample = ROOT.SH.SampleLocal(sampleName.Data())

        localfiles = getFilesFromCommandLine(options.Input,options.InputList)

        for localfile in localfiles :
            sample.add(localfile)

        # Add Sample to SampleHandler
        myhandler.add(sample)

    # Set number of events per worker
    if options.BatchGE or options.BatchCondor :
        SetEventsPerWorker(myhandler,conf)

    PrintSampleSummary(myhandler)

    # myhandler.printContent()
    myhandler.setMetaString('nc_tree','CollectionTree') # must be done AFTER scanDir

    # Give SampleHandler to Job
    myjob.sampleHandler(myhandler)


    ####################################
    # Set up Drivers
    ####################################

    if options.Grid :
        driver = ROOT.EL.PrunDriver()
        SetGridOpts(driver,conf)

    elif options.BatchGE :
        driver = ROOT.EL.GEDriver()
        SetGridEngineOpts(driver,conf)

    elif options.BatchCondor :
        driver = ROOT.EL.CondorDriver()
        SetCondorOpts(driver,conf,options)

    else :
        driver = ROOT.EL.DirectDriver()

    if conf.getBool('SubmitAndDownload', False) :
        # this will submit and automatically download the dataset once the job is finished.
        driver.submit(myjob, submitDir)
    else :
        # this will only submit the job -- the user will have to download using rucio.
        driver.submitOnly(myjob, submitDir)

    print 'done'
    return 0

if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser()

    # Give the name of the algorithm
    p.add_option('--Alg',type='string',default='',dest='Alg',help='algorithm (SkimAndSlim, HGamCutflowAndMxAOD, ...')

    # Options to input data. With local running, list of files. With grid running, list of dsets.
    p.add_option('--Input',type='string',default='',dest='Input',help='Input files or datasets (comma-separated)')
    p.add_option('--InputList',type='string',default='',dest='InputList',help='Input file/dataset list')

    # For users with direct access to their grid files
    p.add_option('--GridDirect',action='store_true',default=False,dest='GridDirect',help='Get data using GridDirect option')

    # Other options
    p.add_option('-n','--NumEvents',type='int',default=-1,dest='NumEvents',help='Number of events to run over')
    p.add_option('--SkipEvents',type='int',default=-1,dest='SkipEvents',help='Number of events to skip')

    p.add_option('--BaseConfig',type='string',default=None,dest='BaseConfig',help='Base config file')
    p.add_option('--Config',type='string',default=None,dest='Config',help='config file')
    p.add_option('--OutputDir',type='string',default=None,dest='OutputDir',help='output directory name')

    p.add_option('--BatchGE',action='store_true',default=False,dest='BatchGE',help='use SGE driver')
    p.add_option('--BatchCondor',action='store_true',default=False,dest='BatchCondor',help='use Condor driver')
    p.add_option('--Grid',action='store_true',default=False,dest='Grid',help='use the Grid')

    # Condor Options
    p.add_option('--Condor_shellInit',type='string',default='',dest='Condor_shellInit',help='Condor: any required shellInit commands')
    p.add_option('--optCondorConf',type='string',default='',dest='optCondorConf',help='Configuration lines to be added to the Condor config file')
    p.add_option('--Condor_UseLD_LIBRARY_PATH',action='store_true',default=False,dest='Condor_UseLD_LIBRARY_PATH',help='Export LD_LIBRARY_PATH to job (required sometimes)')

    # SGE + Condor Options
    p.add_option('--nc_EventLoop_EventsPerWorker',type='int',default=-1,dest='nc_EventLoop_EventsPerWorker',help='Number of events per subjob')

    # Lots of grid options:
    for opt in GetGridOptionsDouble() + GetGridOptionsString() + ['nc_grid_filter'] :
        p.add_option('--%s'%(opt),type='string',default=None,dest=opt,help=opt)

    # Grid Output grid tag:
    p.add_option('--GridTag',type='string',default=None,dest='GridTag',help='GridTag (user.<UserName>.<Tag>)' )

    options,args = p.parse_args()

    if not options.Alg :
        Error('You must specify an algorithm (using --Alg)')

    if not options.Config :
        Error('You must specify a config file (using --Config)')

    if (not options.Input) and (not options.InputList) :
        Error('No input (--Input, --InputList) specified')

    if options.Grid and not os.getenv('ATLAS_LOCAL_PANDACLIENT_PATH') :
        Error('You must execute \"lsetup panda\" (or \"localSetupPandaClient\") to run on the grid')

    if options.GridDirect and not os.getenv('RUCIO_HOME') :
        Error('Please execute localSetupRucioClients and voms-proxy-init')

    if options.Grid and options.GridDirect :
        Error('Specified --Grid and --GridDirect. GridDirect option should not be used with the grid')

    main(options,args)