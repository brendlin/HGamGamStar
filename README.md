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