
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

    map_newnames = dict()

    # Only remove the R-tag if we figured out the MC campaign (see below)
    def removeRtags(nm) :
        nm = re.sub(r'r[0-9][0-9][0-9][0-9]_','',nm)
        nm = re.sub(r'r[0-9][0-9][0-9][0-9][0-9]_','',nm)
        return nm

    for sample in samplehandler :
        name = sample.name()
        if tag :
            name = '%s.%s'%(name,tag)

        if gridtag and issubclass(type(sample),ROOT.SH.SampleGrid) :
            name = '%s.%s'%(gridtag,name)
        
        # General
        name = name.replace('.deriv.'       ,'.')
        name = name.replace('.DAOD_HIGG1D2.','.')
        name = name.replace('.DAOD_HIGG1D1.','.')

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
            else :
                name = name.replace('mc16_13TeV.','mc16.')

        # shorten data name:
        name = name.replace('data15_13TeV','data15')
        name = name.replace('data16_13TeV','data16')
        name = name.replace('data17_13TeV','data17')
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

        # Trick to rename samples: "merge" them into themselves
        # (This is because they cannot be renamed if they already belong to a SampleHandler)
        print 'Renaming \"%s\" to \"%s\"'%(k,map_newnames[k])
        ROOT.SH.mergeSamples(samplehandler,map_newnames[k],k)

    return
