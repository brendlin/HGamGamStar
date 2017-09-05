################################################################################
# Package: HGamGamStar
################################################################################

# Declare the package name:
atlas_subdir( HGamGamStar )

# Declare the package dependencies:
atlas_depends_on_subdirs( PUBLIC
                          HGamCore/HGamAnalysisFramework
                          Event/xAOD/xAODCutFlow
                          )

# Find the needed external(s):
find_package( ROOT COMPONENTS Core RIO Hist Tree )

# build a CINT dictionary for the library
atlas_add_root_dictionary ( HGamGamStarLib HGamGamStarCintDict
                            ROOT_HEADERS HGamGamStar/HggStarCutflowAndMxAOD.h Root/LinkDef.h
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
                   )

# Install files from the package:
atlas_install_data( data/* )

# Install scripts from the package:
# atlas_install_scripts( scripts/*.sh scripts/*.py )

atlas_add_executable( runHggStarCutflowAndMxAOD util/runHggStarCutflowAndMxAOD.cxx
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} HGamGamStarLib )
