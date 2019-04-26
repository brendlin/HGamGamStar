
#-------------------------------------------------------------------------
def renameSample(sh,oldSampleName,newSampleName) :
    # Trick to rename samples: "merge" them into themselves
    # (This is because they cannot be renamed if they already belong to a SampleHandler)
    # Written in a way that preserves the Nevents sample metadata

    import ROOT

    new_sh = ROOT.SH.SampleHandler()
    renamedsample = ROOT.SH.SampleLocal(newSampleName) # the merged sample

    # Typedef some c++ functions
    MetaVector_Long64_t = ROOT.SH.MetaVector("Long64_t")
    vector_Long64_t = ROOT.vector("Long64_t")

    vec_file_nentries = vector_Long64_t()

    for sample in sh :

        if sample.name() == oldSampleName :
            print 'Renaming \"%s\" to \"%s\"'%(oldSampleName,newSampleName)
            meta = sample.meta()

            for i in range(sample.numFiles()) :
                renamedsample.add(sample.fileName(i))
                vec_file_nentries.push_back(meta.get(ROOT.SH.MetaFields.numEventsPerFile).value[i])

            tot_entries = meta.get(ROOT.SH.MetaFields.numEvents).value
            renamedsample.meta().setDouble(ROOT.SH.MetaFields.numEvents,tot_entries)

        else :
            new_sh.add(sample)

    if vec_file_nentries.size() != renamedsample.makeFileList().size() :
        Error("mergeSamplesByRunRange: the nentries vector and filename vector sizes are different.")

    if renamedsample.numFiles() == 0 :
        Error("renameSample: There are no files in sample %s."%(renamedsample.name()))

    file_nentries_meta = MetaVector_Long64_t(ROOT.SH.MetaFields.numEventsPerFile,vec_file_nentries)
    renamedsample.meta().addReplace(file_nentries_meta)
    new_sh.add(renamedsample)

    return new_sh

#-------------------------------------------------------------------------
def mergeByRunRange(sh,year_tag,tag,periodName,firstRun,lastRun_inclusive) :
    # Meant to re-implement mergeSamples in ToolsJoin.cxx, but in a way that preserves
    # the Nevents sample metadata
    import re
    import ROOT

    # Typedef some c++ functions
    MetaVector_Long64_t = ROOT.SH.MetaVector("Long64_t")
    vector_Long64_t = ROOT.vector("Long64_t")

    sampleName = '%s.%s%s'%(year_tag,periodName,tag)

    new_sh = ROOT.SH.SampleHandler()
    mergedsample = ROOT.SH.SampleLocal(sampleName) # the merged sample
    vec_file_nentries = vector_Long64_t()
    tot_entries = 0

    for sample in sh :

        do_merge = True

        try :
            runNumber = int(sample.name().split('.')[1])
        except ValueError :
            do_merge = False

        do_merge = do_merge and year_tag in sample.name() # e.g. 'data17_13TeV'
        do_merge = do_merge and runNumber >= firstRun
        do_merge = do_merge and runNumber <= lastRun_inclusive

        if do_merge :
            meta = sample.meta()
            tot_entries += meta.get(ROOT.SH.MetaFields.numEvents).value
            for i in range(sample.numFiles()) :
                mergedsample.add(sample.fileName(i))
                vec_file_nentries.push_back(meta.get(ROOT.SH.MetaFields.numEventsPerFile).value[i])

        else :
            new_sh.add(sample)

    if vec_file_nentries.size() != mergedsample.makeFileList().size() :
        Error("mergeSamplesByRunRange: the nentries vector and filename vector sizes are different.")

    file_nentries_meta = MetaVector_Long64_t(ROOT.SH.MetaFields.numEventsPerFile,vec_file_nentries)
    mergedsample.meta().addReplace(file_nentries_meta)
    mergedsample.meta().setDouble(ROOT.SH.MetaFields.numEvents,tot_entries)
    new_sh.add(mergedsample)

    return new_sh

#-------------------------------------------------------------------------
def SetSampleNames(samplehandler,tag='',gridtag='') :
    #
    # This sets the sample name, including the name of the output files.
    # Works best with GridDirect.
    # This is NOT the same as setting nc_outputSampleName, which has to do
    # with the format of the grid output file names.
    #

    import re
    import ROOT

    if tag :
        tag = '.'+tag

    if gridtag :
        gridtag = gridtag + '.'

    map_newnames = dict()

    # Only remove the R-tag if we figured out the MC campaign (see below)
    def removeRtags(nm) :
        nm = re.sub(r'r[0-9][0-9][0-9][0-9]_','',nm)
        nm = re.sub(r'r[0-9][0-9][0-9][0-9][0-9]_','',nm)
        return nm

    merge_data15,merge_data16,merge_data17,merge_data18 = False,False,False,False

    for sample in samplehandler :

        name = sample.name()
        name = name + tag

        # non-Grid data merging is treated below, not here
        if not issubclass(type(sample),ROOT.SH.SampleGrid) :

            merge_data15 = ('data15_13TeV' in name)
            merge_data16 = ('data16_13TeV' in name)
            merge_data17 = ('data17_13TeV' in name)
            merge_data18 = ('data18_13TeV' in name)

            if merge_data15 or merge_data16 or merge_data17 or merge_data18 :
                continue

        if gridtag and issubclass(type(sample),ROOT.SH.SampleGrid) :
            name = gridtag + name
        
        # General
        name = name.replace('.deriv.'       ,'.')
        name = name.replace('.DAOD_HIGG1D2.','.')
        name = name.replace('.DAOD_HIGG1D1.','.')

        name = name.replace('user.amorley.user.amorley','user.amorley')
        name = name.replace('_EXT0','')

        # mc16a / mc16c / mc16d / mc16e
        if 'mc16_13TeV' in name :
            if 'r9364' in name :
                name = name.replace('mc16_13TeV.','mc16a.')
                name = removeRtags(name)
            elif 'r9781' in name :
                name = name.replace('mc16_13TeV.','mc16c.')
                name = removeRtags(name)
            elif 'r10201' in name :
                name = name.replace('mc16_13TeV.','mc16d.')
                name = removeRtags(name)
            elif 'r10724' in name :
                name = name.replace('mc16_13TeV.','mc16e.')
                name = removeRtags(name)
            else :
                name = name.replace('mc16_13TeV.','mc16.')

        name = name.replace('.physics_Main.','.')

        # Generator - Powheg
        name = name.replace('PowhegPythia8EvtGen_NNLOPS'  ,'PowhegPy8_NNLOPS')
        name = name.replace('PowhegPythia8EvtGen'         ,'PowhegPy8'       )
        name = name.replace('PowhegHerwig7EvtGen'         ,'PowhegH7'        )

        # Generator - aMcAtNlo
        name = name.replace('aMcAtNloPythia8EvtGen','aMCnloPy8' )
        name = name.replace('aMcAtNloHppEG'        ,'aMCnloHwpp')
        
        # PDF
        name = name.replace('_NNPDF30_AZNLO_','_')
        name = name.replace('_NNLOPS_nnlo_30_','_')
        name = name.replace('_A14NNPDF23LO_','_')
        name = name.replace('_UEEE5_CTEQ6L1_CT10ME_','_')
        name = name.replace('_NNPDF30_AZNLOCTEQ6L1_','_')

        # Remove e-tags, s-tags, r-tags
        name = re.sub(r'e[0-9][0-9][0-9][0-9]_','',name)
        name = re.sub(r's[0-9][0-9][0-9][0-9]_','',name)

        # Miscellaneous
        name = name.replace('Zy_Zll' ,'Zllgam')
        name = name.replace('_gamgam',''      )

        if issubclass(type(sample),ROOT.SH.SampleGrid) :
            sample.setMetaString("nc_outputSampleName",name)
        else :
            map_newnames[sample.name()] = name


    for k in map_newnames.keys() :

        # Do not rename if it's a grid sample
        if issubclass(type(samplehandler.get(k)),ROOT.SH.SampleGrid) :
            continue

        samplehandler = renameSample(samplehandler,k,map_newnames[k])

    # Merge data runs into separate periods (for file sizes that fit on eos).
    # Not possible for grid datasets.
    if merge_data15 :
        samplehandler = mergeByRunRange(samplehandler,'data15_13TeV',tag,'periodAllYear',0,9999999)

    if merge_data16 :
        samplehandler = mergeByRunRange(samplehandler,'data16_13TeV',tag,'periodAtoG',297730,306714)
        samplehandler = mergeByRunRange(samplehandler,'data16_13TeV',tag,'periodItoL',307124,311481)

    if merge_data17 :
        samplehandler = mergeByRunRange(samplehandler,'data17_13TeV',tag,'periodAtoE',324320,334779)
        samplehandler = mergeByRunRange(samplehandler,'data17_13TeV',tag,'periodFtoK',334842,340453)

    if merge_data18 :
        samplehandler = mergeByRunRange(samplehandler,'data18_13TeV',tag,'periodAtoK',348197,357000)
        samplehandler = mergeByRunRange(samplehandler,'data18_13TeV',tag,'periodLtoR',357001,364486)

    return samplehandler
