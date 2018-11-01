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
| 3, HIGGS_LEP_DALITZ	   | pass if (Higgs dalitz event)<br> pass if !(Higgs event)<br>(Dalitz = exactly 2&ell;+1&gamma; children)  | eventIsNonHyyStarHiggs() in<br>HGamGamStar/HggStarVariables.cxx | " " |
| 4, DUPLICATE             | eventHandler()->isDuplicate()  | HGamAnalysisFramework/EventHandler.cxx | " " |
| 5, GRL                   | eventHandler()->passGRL()      | HGamAnalysisFramework/EventHandler.cxx | " " |
| 6, TRIGGER               | eventHandler()->passTriggers() | HGamAnalysisFramework/EventHandler.cxx | " " |
| 7, DQ                    | pass LAr, Tile, SCT DQ         | HGamAnalysisFramework/EventHandler.cxx | " " |
| 8, VERTEX                | eventHandler()->passVertex()   | HGamAnalysisFramework/EventHandler.cxx | " " |
| **Get ele container**  | m_preSelElectrons, pass OQ, HV,<br> p<sub>T</sub> > 4.5 GeV, &#124;&eta;&#124;<2.47, | ElectronHandler.cxx<br>PtPreCutGeV in HggStarMxAOD.config | N/A |
| **Get trk container**  | m_preSelTracks, preselected GSF <br> tracks (p<sub>T</sub>>0.3 GeV, &#124;&eta;&#124;<2.47) <br> nSi+DeadSens &ge; 7, nPix+DeadSens &ge; 2 <br>assoc. to a presel ele | TrackHandler.cxx | N/A |
| **Get muon container** | m_preSelMuons, Medium PID,<br> p<sub>T</sub>>3 GeV, &#124;&eta;&#124;<2.7 | HggStarMxAOD.config<br>&eta; cut is in MuonHandler.cxx | N/A |
| 9, TWO_SF_LEPTONS        | NpreSelTracks &ge; 2 or<br> NpreSelMuons &ge; 2 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Get &gamma; container** | m_preSelPhotons, OQ, cleaning, HV,<br> Loose ID, p<sub>T</sub>>10 GeV, &#124;&eta;&#124;<2.37<br> no crack, AuthorAmbiguous allowed | PhotonHandler.cxx | N/A
|10, ONE_LOOSE_GAM         | Nloose &ge; 1 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Photon choice**      | Selected photon is highest-p<sub>T</sub><br> presel photon | HiggsGamGamStarCutflowAndMxAOD.cxx | N/A |
| **Z boson assignment** | Pick highest vector-sum p<sub>T</sub> SFOS pair<br>(Muons or tracks). Channels: <br>DIMUON=1, RESOLVED_DIELECTRON=2,<br>MERGED_DIELECTRON=3,<br>AMBIGUOUS_DIELECTRON=4. | HggStarVariables.cxx,<br>HiggsGamGamStarCutflowAndMxAOD.cxx | N/A |
|12, ZBOSON_ASSIGNMENT     | nSFOS &ge; 1 (muon or track pairs). | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Get jet container**  | Selected jets | JetHandler.cxx | N/A |
| **Overlap removal**    | Remove e if &Delta;R(e,&gamma;)<0.4<br> Remove jet if &Delta;R(j,&gamma;)<0.4<br>Remove jet if &Delta;R(j,e)<0.2<br>Remove e if &Delta;R(j,e)<0.4<br>Remove &mu; if &Delta;R(&mu;,&gamma;)<0.4<br>Remove &mu; if &Delta;R(&mu;,j)<0.4 | HGamAnalysisFramework/<br>OverlapRemovalHandler.cxx | " " |
|13, TWO_SF_LEPTONS_POSTOR | &ge; 1 of the SFOS pairs survives OR.<br>At this point, &mu;&mu; is preferred;<br>the remaining channels are exclusive. | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|14, BAD_MUON              | Reject events with presel BadMuons | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|15, ONE_PHOTON_POSTOR     | Selected photon survives OR | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|16, TRIG_MATCH            | Currently disabled in<br>HggStarMxAOD.config | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|17, LEP_MEDID             | Resolved e: Medium (resolved e),<br>Merged e: Merged ID,<br>&mu;: Medium | HGamAnalysisFramework/HGamRel21.config<br>for merged objects | " " |
|18, LEP_IP                | Resolved e: d<sub>0</sub>/&sigma;<sub>d0</sub> < 5, &#124;z<sub>0</sub>sin&theta;&#124; < 0.5<br>Merged e: same but only for first track! Fix!<br>&mu;: d<sub>0</sub>/&sigma;<sub>d0</sub> < 3, &#124;z<sub>0</sub>sin&theta;&#124; < 0.5 | Defaults in ElectronHandler.cxx, MuonHandler.cxx | " " |
|19, LEP_ISO               | Resolved e: CloseByCorrected Loose<br>Merged e: Loose<br>&mu;: GradientLoose (no CloseBy correction!) | Merged e: specially done in<br>HiggsGamGamStarCutflowAndMxAOD.cxx.<br>Resolved e/&mu;: HGamAnalysisFramework/HGamRel21.config | " " |
|20, GAM_TIGHTID           | Photon passes Tight | HGamAnalysisFramework/HGamRel21.config | " " |
|21, GAM_ISOLATION         | Photon passes FixedCutLoose | HGamAnalysisFramework/HGamRel21.config | " " |
|22, ZMASSCUT              | m<sub>ll</sub> < 45 GeV<br> m<sub>trktrk</sub> < 45 GeV for MERGED | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|23, LLGMASSCUT            | 105 < m<sub>ll&gamma;</sub> && m<sub>ll&gamma;</sub> < 160 GeV | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|24, PASSALL               | Everything above passes | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
