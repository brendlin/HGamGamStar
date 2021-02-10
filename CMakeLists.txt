################################################################################
# Package: HGamGamStar
################################################################################

# Declare the package name:
atlas_subdir( HGamGamStar )

# Declare the package dependencies:
atlas_depends_on_subdirs( PUBLIC
                          HGamCore/HGamAnalysisFramework
                          Event/xAOD/xAODCutFlow
                          PhysicsAnalysis/AnalysisCommon/FsrUtils
                          PhysicsAnalysis/HiggsPhys/Run2/HZZ/Tools/ZMassConstraint
                          PhysicsAnalysis/AnalysisCommon/IsolationSelection/IsolationSelection
                          InnerDetector/InDetRecTools/InDetTrackSelectionTool
                          InnerDetector/InDetRecTools/TrackVertexAssociationTool
                          PhysicsAnalysis/ElectronPhotonID/ElectronPhotonShowerShapeFudgeTool
                          )

# Find the needed external(s):
find_package( ROOT COMPONENTS Core RIO Hist Tree )

# build a CINT dictionary for the library
atlas_add_root_dictionary ( HGamGamStarLib HGamGamStarCintDict
                            ROOT_HEADERS HGamGamStar/ZyCutflowAndMxAOD.h
                            HGamGamStar/ZyUnfoldingInput.h 
			    HGamGamStar/ZyyUnfoldingInput.h
                            HGamGamStar/HiggsGamGamStarCutflowAndMxAOD.h
                            HGamGamStar/RadiativeZCutflowAndMxAOD.h
                            HGamGamStar/MergedElectronMxAOD.h
                            HGamGamStar/TrackHandler.h
                            HGamGamStar/ExtraHggStarObjects.h
                            HGamGamStar/MergedElectronID.h
                            HGamGamStar/AngularPosition.h
                            HGamGamStar/TrackModel.h
                            HGamGamStar/MxAODTool.h Root/LinkDef.h
                            EXTERNAL_PACKAGES ROOT
                            )

# build a shared library
atlas_add_library( HGamGamStarLib
                   HGamGamStar/*.h Root/*.cxx ${HGamGamStarCintDict}
                   PUBLIC_HEADERS HGamGamStar
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
                   HGamAnalysisFrameworkLib
                   xAODCutFlow
                   FsrUtilsLib
                   ZMassConstraintLib
                   IsolationSelectionLib
                   InDetTrackSelectionToolLib
                   TrackVertexAssociationToolLib
                   ElectronPhotonShowerShapeFudgeToolLib
                   )


# Install files from the package:
atlas_install_data( data/* )


# Install scripts from the package:
# atlas_install_scripts( scripts/*.sh scripts/*.py )

atlas_add_executable( runHggStarCutflowAndMxAOD util/runHggStarCutflowAndMxAOD.cxx
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} HGamGamStarLib )

atlas_install_scripts( scripts/runJob.py
                       scripts/runJob_condor_resubmit_failed
                       scripts/runJob_condor_resubmit_killed
                       scripts/runJob_merge
                       )
