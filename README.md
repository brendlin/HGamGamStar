H&rarr;&gamma;&gamma;* MxAOD Code Description
========================
This is a description of how to run the MxAOD code for the H&rarr;&gamma;&gamma;* analysis.

Initial setup
--------

To set up the ATLAS environment, run or put in your `.bash_profile` (sometimes called `.profile`) startup file:
    
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    alias 'setupATLAS=source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
    setupATLAS

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

To compile in Rel 21.2.109, do:

    asetup AnalysisBase,21.2.109,here # only needed once per login session
    cd $TestArea/../build
    cmake ../source
    cmake --build .
    
To be able to run executables, run (only once per login session):

    source $TestArea/../build/$CMTCONFIG/setup.sh    

Running - "The Old Way"
---------
(Note: if you are new to the code, we suggest using the new tools! The "old way" is mentioned here for those who are already used to it, and don't want to change.)

To see more details on how to run, look at `HGamCore/HGamAnalysisFramework/Root/RunUtils.cxx`. Make a list of files you want to run, and specify `InputFileList`. If you have a single file, specify `InputFile`:

    cd $TestArea/../run
    runHggStarCutflowAndMxAOD  InputFileList: myFileList.txt  HGamGamStar/HggStarMxAOD.config
    
The jobs follow the same convention as [the HGam tutorial twiki.](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HGamAnalysisFrameworkTutorial).

Running - New Tools!
---------
A new run tool, based in Python, allows for more user-friendly job running. To run a test job, do:

    cd $TestArea/../run
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

### How (and why) to use the **--GridDirect** option

SampleHandler has a function for converting a DSID into a list of files that are local to a specified LOCALGROUPDISK.
In this way, all you need to do is specify the DSID and a list of files will be made for you.
**Note that --GridDirect can be time consuming -- so the code will save the local files it finds to a text file, and
put it in a directory called GridDirectFiles. When you rerun on the same DSID and specify --GridDirect, the code
will find the file that it saved in the previous run, and therefore save you a lot of time on the second run.**

You will need to export a few variables to your environment so that the code knows how to access your local group disk.
(You can put these in your `.bash_profile` startup script too):

    export LOCALGROUPDISK="DESY-HH_LOCALGROUPDISK"
    export GRIDDIRECT_FROM="root://dcache-atlas-xrootd.desy.de:1094//"
    export GRIDDIRECT_TO="/"

The last two variables relate to how to interpret the file paths on your institute's local group disk. If you do not know the file protocol, try running:

    rucio list-file-replicas --protocols root --rse DESY-HH_LOCALGROUPDISK SOMEDATASET
    
where `SOMEDATASET` is an existing dataset on your local group disk.
For DESY, the output of this command shows file paths that look like `root://dcache-atlas-xrootd.desy.de:1094//some/file/path`.
The `GRIDDIRECT_FROM` and `GRIDDIRECT_TO` environment variables are telling the code to replace `root://dcache-atlas-xrootd.desy.de:1094//` with `/`
in order to convert to a local file path. (Note that DESY is the only environment where this was tested, so please contact the developers if
you can't figure out your system.)

Once you've set up these environment variables, you can use the built-in functionality to
automatically convert the DSIDs in your localgroupdisk to a list of files.

### Running with Condor

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
    cd $TestArea/../run
    runJob.py --InputList Samples.txt --OutputDir MyOutputDir --Alg HiggsGamGamStarCutflowAndMxAOD --Config HGamGamStar/HggStarMxAOD.config --BatchCondor --Condor_UseLD_LIBRARY_PATH --GridDirect --nc_EventLoop_EventsPerWorker 100000 --optCondorConf 'Requirements = ( OpSysAndVer == "CentOS7" )'
    ```

### Rerunning failed / killed Condor jobs

Did any of your Condor jobs fail? You can check by running:

    find MyOutputDir/. | grep fail

In the event that a few jobs failed for "transient reasons" (e.g. there is no inherent bug in the code), you can restart the individual jobs.
To do this, execute:

    cd MyOutputDir/..
    runJob_condor_resubmit_failed MyOutputDir

This will rerun the specific jobs that failed. If they fail again, then you can rerun this script.

Finally, if your jobs stalled for some reason and need to be killed, then **first kill them** and then
you can retry the jobs using the following:

    cd MyOutputDir/..
    runJob_condor_resubmit_killed MyOutputDir

### Running on the Grid

To run on the grid, you must specify a **GridTag** as well as a **ProdTag** via the command-line (or config). An example is below (you can use
**--Input** or **--InputList** to specify the samples):

    runJob.py --Input mc16_13TeV.345961.PowhegPythia8EvtGen_NNLOPS_nnlo_30_ggH125_gamstargam.deriv.DAOD_HIGG1D2.e6740_e5984_s3126_r10201_r10210_p3415 --Alg HiggsGamGamStarCutflowAndMxAOD --Config HGamGamStar/HggStarMxAOD.config --Grid --GridTag user.brendlin --ProdTag ysy00X

You can also specify any of the other normal grid running options via command-line or in your config file
(e.g `--nc_nFilesPerJob`, `--nc_destSE`, etc.).

### Merging the files of an EventLoop job

Once the jobs above are complete, you can merge the files using the `EL::Driver::wait` command.
The easiest way to do this is to launch this via command line, with:

    cd MyOutputDir/..
    runJob_merge MyOutputDir

Note that this will wait for all jobs to finish, and then merge the root files. The resulting MxAODs will be in the directory `MyOutputDir/data-MxAOD`.
If a job failed at any point, an error will be thrown and the merging will be paused. See below for how to re-run the jobs that failed.
(You can simply re-run the above command after your jobs finish successfully, and the merging should start back up where it left off.)

MxAOD Production
=================

Where is the data? There is a script [getSamplesFromGRL.py](https://gitlab.cern.ch/atlas-hgam-sw/HGamCore/blob/master/HGamTools/scripts/getSamplesFromGRL.py)
in HGamCore that creates a list of data DSIDs from some input GRLs. The resulting files are saved in the 
`HGamTools/data/input` directory. We copy the data DSIDs from those files (HIGG1D2), but save them ourselves in [HGamGamStar/data/input](https://gitlab.cern.ch/brendlin/HGamGamStar/tree/master/data/input).

To run all MxAOD production on Condor, 
first **make sure the code is fully committed, tagged in git, and that the event selection description 
[HggStarEventSelection.md](https://gitlab.cern.ch/brendlin/HGamGamStar/blob/master/HggStarEventSelection.md) is fully up-to-date.**
Then do (specifying an appropriate ProdTag):

    prodtag=ysy00X
    for DS in data15_13TeV data16_13TeV data17_13TeV data18_13TeV mc16a_HIGG1D2 mc16d_HIGG1D2 mc16e_HIGG1D2; do
    runJob.py --InputList HGamGamStar/input/$DS.txt --OutputDir ${DS}_${prodtag} --Alg HiggsGamGamStarCutflowAndMxAOD --Config HGamGamStar/HggStarMxAOD.config --BatchCondor --Condor_UseLD_LIBRARY_PATH --GridDirect --nc_EventLoop_EventsPerWorker 100000 --ProdTag $prodtag --optCondorConf 'Requirements = ( OpSysAndVer == "CentOS7" )';
    done;

Wait for all jobs to complete. Then merge using the following commands (it is recommended to run these one-by-one; in case a job failed, follow the rerun procedure outlined above):

    for DS in data15_13TeV data16_13TeV data17_13TeV data18_13TeV mc16a_HIGG1D2 mc16d_HIGG1D2 mc16e_HIGG1D2; do
        runJob_merge ${DS}_${prodtag}
    done;

When all jobs have merged, put the MxAOD output files (located in the directory `data-MxAOD`) into a directory on EOS and the DESY dust (if applicable),
and add the production details to the section below.

Information on MxAOD Productions
----------------
| **MxAOD tag** | **HGamGamStar tag** | **Base Release** | **Event selection** | **Location** |
| ------------- | ----------------- | -------------- | ----------------- | ---------- |
| ysy001 | [ysy001](https://gitlab.cern.ch/brendlin/HGamGamStar/tags/ysy001) | 21.2.25 | [ysy001 HggStarEventSelection.md](https://gitlab.cern.ch/brendlin/HGamGamStar/blob/0bf0779154ff38eade37b64e611707cfd77989a6/HggStarEventSelection.md) | DESY: /nfs/dust/atlas/user/brendlik/eos/ysy/ysy001<br>EOS: /eos/user/b/brendlin/ysy/ysy001  |
| ysy002 | [ysy002](https://gitlab.cern.ch/brendlin/HGamGamStar/tags/ysy002) | 21.2.25 | Not fully documented yet | DESY: /nfs/dust/atlas/user/brendlik/eos/ysy/ysy002<br>EOS: /eos/user/b/brendlin/ysy/ysy002 |
| ysy003 | [ysy003](https://gitlab.cern.ch/brendlin/HGamGamStar/tags/ysy003) | 21.2.25 | [ysy003 HggStarEventSelection.md](https://gitlab.cern.ch/brendlin/HGamGamStar/blob/fc0b714ed0354e37fd8239e7ca732ed25323ca12/HggStarEventSelection.md) | DESY: /nfs/dust/atlas/user/brendlik/eos/ysy/ysy003<br>EOS: /eos/user/b/brendlin/ysy/ysy003  |
| ysy005 | [ysy005](https://gitlab.cern.ch/brendlin/HGamGamStar/tags/ysy005) | 21.2.56 | [ysy005 HggStarEventSelection.md](https://gitlab.cern.ch/brendlin/HGamGamStar/blob/51e7884246726940caf6effe7465154c2343f470/HggStarEventSelection.md) | DESY: /nfs/dust/atlas/user/brendlik/eos/ysy/ysy005<br>EOS: /eos/user/b/brendlin/ysy/ysy005  |
| ysy011 | [ysy011](https://gitlab.cern.ch/brendlin/HGamGamStar/tags/ysy011) | 21.2.99 | Browse HggStarEventSelection.md in tag | DESY: /nfs/dust/atlas/user/brendlik/eos/ysy/ysy011<br>EOS: /eos/user/b/brendlin/ysy/ysy011  |

Updating the HGamCore Tag
================
If you are a *user* who is trying to update the HGamGamStar package, including an update to submodule tags, first make sure your
submodules have no local edits. Then do:

    git pull
    git submodule update --recursive
    
If you are a developer trying to update to a new HGamCore tag, do e.g.:

    cd HGamGamStar/HGamCore
    git checkout v1.3.0-h019
    git submodule update --init --recursive
    # don't forget to commit to the HGamGamStar repository.

Unfortunately, the commit hash is saved instead of the tag name, so you will
have to match the has to the tag name if you want to see which one is currently
used. For convenience, here is a list of the tags and their corresponding
hashes used in this package:

| HGamCore tag | HGamCore commit hash | AnalysisBase release | Notes |
| ------------ | ----------- | ---------- | ---------- |
| v1.1.2-h017  | cb2c8452    | ?          |          |
| v1.2.0-h018  | 04d6a778    | ?          |          |
| v1.3.0-h019  | 76779324    | ?          |          |
| v1.5.5-h021  | b941a3d8    | 21.2.25    |          |
| v1.8.1-h024  | 9f868a8b    | 21.2.56    |          |
| v1.8.15-h024 | f7f097c8    | 21.2.56    |          |
| v1.8.34-h024 | 8551d8d5    | 21.2.56    |          |
| v1.8.47-h024-ttHcp | 5312e5d9    | 21.2.56    | Update for trig-match mem leak |
| master commit | f2e09974   | 21.2.99    | Updated muon and jet recommendations |
| untagged commit | 9d640507 | 21.2.109   | Update our code for p4061, p4062 |
