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

    git clone ssh://git@gitlab.cern.ch:7999/brendlin/HGamGamStar.git
    git clone --recursive ssh://git@gitlab.cern.ch:7999/atlas-hgam-sw/HGamCore.git 

Compiling
---------

To compile in Rel 21.2.4, do:

    asetup AnalysisBase,21.2.4,here # only needed once per login session
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
