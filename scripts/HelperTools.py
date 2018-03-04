
#-------------------------------------------------------------------------
def SetOutputDatasetNames(samplehandler,tag) :
    import re
    for sample in samplehandler :
        name = '%s.%s'%(tag,sample.name())
        
        # General
        name = name.replace('.deriv.'       ,'.')
        name = name.replace('mc16_13TeV.'   ,'' )
        name = name.replace('.DAOD_HIGG1D2.','.')
        name = name.replace('.DAOD_HIGG1D1.','.')

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
        name = re.sub(r'r[0-9][0-9][0-9][0-9]_','',name)

        # Miscellaneous
        name = name.replace('Zy_Zll' ,'Zllgam')
        name = name.replace('_gamgam',''      )

        sample.setMetaString("nc_outputSampleName",name)

    return
