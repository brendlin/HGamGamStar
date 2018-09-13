H&rarr;&gamma;&gamma;* MxAOD Code Description
========================
This is a description of how to run the MxAOD code for the H&rarr;&gamma;&gamma;* analysis.

Checking Out the code
--------
If you are running outside of the lxplus network, it is suggested that you create a kerberos ticket by running "kinit username@CERN.CH" and authenticating with your password. Then, set up your area -- for instance:

    mkdir testarea
    cd testarea
    mkdir build run source

To checkout the necessary packages, do:

    cd source
    git clone --recursive ssh://git@gitlab.cern.ch:7999/brendlin/HGamGamStar.git

Compiling
---------

To compile in Rel 21.2.22, do:

    asetup AnalysisBase,21.2.22,here # only needed once per login session
    cd $TestArea/../build
    cmake ../source
    cmake --build .
    
To be able to run executables, run (only once per login session):

    source $TestArea/../build/$CMTCONFIG/setup.sh    

Running
---------
To see more details on how to run, look at `HGamCore/HGamAnalysisFramework/Root/RunUtils.cxx`. Make a list of files you want to run, and specify `InputFileList`. If you have a single file, specify `InputFile`:

    cd $TestArea/../run
    runHggStarCutflowAndMxAOD  InputFileList: myFileList.txt  HGamGamStar/HggStarMxAOD.config
    
The jobs follow the same convention as [the HGam tutorial twiki.](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HGamAnalysisFrameworkTutorial).

Running - New Tools!
---------
A new run tool, based in Python, allows for more user-friendly job running. To run a test job, do:

    runJob.py --InputList MyList.txt --OutputDir Test --Alg HiggsGamGamStarCutflowAndMxAOD --Config HGamGamStar/HggStarMxAOD.config -n 10000

 - **--Alg**: the options are "ZyCutflowAndMxAOD" or "HiggsGamGamStarCutflowAndMxAOD".
 - **Run options supported**: **--BatchGE**, **--BatchCondor** or **--Grid**.
 - **Root files or Datasets?**: If you want to run over datasets on the grid (**--Grid**) or using GridDirect (**--GridDirect**), specify the input using **dataset (DSID) names**. Otherwise, use root file names.
 - **Inputs** can be specified using the command **--Input** (comma-separated list of root files or DSIDs) or **--InputList** (one file/DSID per line, comment a line using #)
 - **-n** (for local tests) - limit yourself to a maximum number of events
 - **--config** (and **--BaseConfig**): same as HGam, e.g. --config HggStarMxAOD.config. BaseConfig is typically set inside the --config, so you don't need to re-specify it here.
 - **--OutputDir**: the output directory where you want the results to go. Must be provided.

Some more details about using this run scheme with Condor:
 - **--nc_EventLoop_EventsPerWorker** is **required** when running on Condor (or another batch system). A good number, I have found, is 100000 (e.g. 100k). Larger number = fewer, longer jobs.
 - **--Condor_shellInit**: In case you want to add any arguments to "shellInit" of your condor job (expert mode)
 - **--optCondorConf**: Another option to add to the condor option configuration (expert mode)
 - **--Condor_UseLD_LIBRARY_PATH**: Some systems refuse to export `LD_LIBRARY_PATH` to your worker nodes, so some code was written to do this for you if it's necessary. (This is a boolean). If you work at DESY, then this is probably necessary.

Putting it all together:

 - You have a list of datasets that you want to run on condor. You make a file Samples.txt that looks like:

    ```
    mc16_13TeV.343981.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_gamgam.deriv.DAOD_HIGG1D1.e5607_e5984_s3126_r9781_r9778_p3404
    mc16_13TeV.343981.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_gamgam.deriv.DAOD_HIGG1D1.e5607_s3126_r9364_r9315_p3404
    # mc16_13TeV.364500.Sherpa_222_NNPDF30NNLO_eegamma_pty_7_15.deriv.DAOD_HIGG1D2.e5928_e5984_s3126_r9781_r9778_p3415
    # mc16_13TeV.364501.Sherpa_222_NNPDF30NNLO_eegamma_pty_15_35.deriv.DAOD_HIGG1D2.e5928_e5984_s3126_r9781_r9778_p3415
    data15_13TeV.00276262.physics_Main.deriv.DAOD_HIGG1D2.r9264_p3083_p3402
    data16_13TeV.00297730.physics_Main.deriv.DAOD_HIGG1D2.r9264_p3083_p3372
    data17_13TeV.00325713.physics_Main.deriv.DAOD_HIGG1D2.f839_m1824_p3372
    ```

 - Then you can run the job on condor with the following command:

    ```
    runJob.py --InputList Samples.txt --OutputDir MyOutputDir --Alg HiggsGamGamStarCutflowAndMxAOD --Config HGamGamStar/HggStarMxAOD.config --BatchCondor --Condor_UseLD_LIBRARY_PATH --GridDirect --nc_EventLoop_EventsPerWorker 100000
    ```

Updating the HGamCore Tag
---------
To update the HGamCore tag, do:
```
cd HGamGamStar/HGamCore
git checkout v1.3.0-h019
git submodule update --init --recursive
# don't forget to commit to the HGamGamStar repository.
```

Unfortunately, the commit hash is saved instead of the tag name, so you will
have to match the has to the tag name if you want to see which one is currently
used. For convenience, here is a list of the tags and their corresponding
hashes used in this package:

| HGamCore tag | HGamCore commit hash |
| ------------ | ----------- |
| v1.1.2-h017  | cb2c8452    |
| v1.2.0-h018  | 04d6a778    |
| v1.3.0-h019  | 76779324    |