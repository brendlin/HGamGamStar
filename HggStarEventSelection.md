H&rarr;&gamma;*&gamma; Event selection
=================

The event selection steps, and the relevant files where the selection occurs,
are detailed below.
Note: *please update this file* if you update aspects of the event selection!

Each step corresponds to a CutEnum defining the cutflow (see `HiggsGamGamStarCutflowAndMxAOD.h`).
"Step 0" defines the number of events in the xAOD (with the enum name `NxAOD`).

| CutEnum,  Name           | Description  | Location | Weight |
| ------------------------ | ------------ | -------- | ------ |
| 0  NxAOD                 | NxAOD from CutBookkeeper | MxAODTool.cxx | CBK sumOfEventWeights() |
| 1  NDxAOD                | NDxAOD from CutBookkeeper | MxAODTool.cxx | CBK sumOfEventWeights() |
| 2  All Events            | All events run over in DAOD | HiggsGamGamStarCutflowAndMxAOD.cxx | weightInitial() in HgammaAnalysis<br>(mcWeight, pileupWeight, vertexWeight) |
| 3, HIGGS_LEP_DALITZ	   | pass if (Higgs event && dalitz)<br> pass if !(Higgs event) | " " | " " |
| 4, DUPLICATE             | eventHandler()->isDuplicate()  | " " | " " |
| 5, GRL                   | eventHandler()->passGRL()      | " " | " " |
| 6, TRIGGER               | eventHandler()->passTriggers() | " " | " " |
| 7, DQ                    | pass LAr, Tile, SCT DQ         | " " | " " |
| 8, VERTEX                | eventHandler()->passVertex()   | " " | " " |
| **Get ele container**  | m_preSelElectrons, pass OQ, HV,<br> p<sub>T</sub> > 4.5 GeV, &eta;<2.47, | ElectronHandler.cxx | N/A |
| **Get trk container**  | m_preSelTracks, preselected GSF <br> tracks (p<sub>T</sub>>0.3 GeV, &eta;<2.47) <br> nSi+DeadSens &ge; 7, nPix+DeadSens &ge; 2 <br>assoc. to a presel ele | TrackHandler.cxx | N/A |
| **Get muon container** | m_preSelMuons, Medium PID,<br> p<sub>T</sub>>3 GeV, &eta;<2.7 | HggStarMxAOD.config | N/A |
| 9, TWO_SF_LEPTONS        | NpreSelTracks &ge; 2 or<br> NpreSelMuons &ge; 2 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Get &gamma; container** | m_preSelPhotons, OQ, cleaning, HV,<br> Loose ID, p<sub>T</sub>>10 GeV, &eta;<2.37<br> no crack, Ambi(?) | PhotonHandler.cxx | N/A
|10, ONE_LOOSE_GAM         | Nloose &ge; 1 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Photon choice**      | Selected photon is highest-p<sub>T</sub><br> presel photon | " " | N/A |
| **Z boson assignment** | Pick highest vector-sum p<sub>T</sub> SFOS pair<br>(Muons or tracks). Channels: <br>DIMUON=1,<br> RESOLVED_DIELECTRON=2,<br>MERGED_DIELECTRON=3,<br>AMBIGUOUS_DIELECTRON=4 | HggStarVariables.cxx | N/A |
|12, ZBOSON_ASSIGNMENT     | nSFOS &ge; 1 (muon or track pairs) | " " | " " |
| **Get jet container**  | Selected jets | JetHandler.cxx | N/A |
| **Overlap removal**    | Need to remember details | " " | " " |
|13, TWO_SF_LEPTONS_POSTOR | &ge; 1 of the SFOS pairs survives OR | " " | " " |
|14, BAD_MUON              | Reject events with presel BadMuons | " " | " " |
|15, ONE_PHOTON_POSTOR     | Selected photon survives OR | " " | " " |
|16, TRIG_MATCH            | Need to remember details | " " | " " |
|17, GAM_TIGHTID           | Photon passes Tight | " " | " " |
|18, GAM_ISOLATION         | Photon passes Isolation | " " | " " |
|19, ZMASSCUT              | m<sub>ll</sub> < 45 GeV<br> m<sub>trktrk</sub> < 45 GeV for MERGED | " " | " " |
|20, LLGMASSCUT            | 105 < m<sub>ll&gamma;</sub> && m<sub>ll&gamma;</sub> < 160 GeV | " " | " " |
|21, PASSALL               | Everything above passes | " " | " " |
