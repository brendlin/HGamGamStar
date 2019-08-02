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
| **Get ele container**  | m_preSelElectrons, pass OQ, HV,<br> p<sub>T</sub> > 4.5 GeV, &#124;&eta;&#124;<2.47 | ElectronHandler.cxx<br>PtPreCutGeV in HggStarMxAOD.config | N/A |
| **Get trk container**  | m_preSelTracks, Assoc. to a presel elec<br>**See note below for other cuts** | TrackHandler.cxx | N/A |
| **Get muon container** | m_preSelMuons, Medium PID,<br> p<sub>T</sub>>3 GeV, &#124;&eta;&#124;<2.7 | HggStarMxAOD.config<br>&eta; cut is in MuonHandler.cxx | N/A |
| 9, TWO_SF_LEPTONS        | NpreSelTracks &ge; 2 or<br> NpreSelMuons &ge; 2 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Get &gamma; container** | m_preSelPhotons, OQ, cleaning, HV,<br> Loose ID, p<sub>T</sub>>20 GeV, &#124;&eta;&#124;<2.37<br> no crack, AuthorAmbiguous allowed | PhotonHandler.cxx | N/A
|10, ONE_LOOSE_GAM         | Nloose &ge; 1 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Photon choice**      | Selected photon is highest-p<sub>T</sub><br> presel photon | HiggsGamGamStarCutflowAndMxAOD.cxx | N/A |
| **Z boson assignment** | Pick highest vector-sum p<sub>T</sub> SFOS pair<br>(&mu;s or tracks). Channels: <br>DIMUON=1, RESOLVED_DIELECTRON=2,<br>MERGED_DIELECTRON=3.<br>&mu; preferred. Use preselection below. | HggStarVariables.cxx,<br>HiggsGamGamStarCutflowAndMxAOD.cxx | N/A |
| **Presel for Z-assignment above** | Res e: VeryLooseLH,<br>lead p<sub>T</sub>>13 GeV<br>Mrgd e: Rhad<0.1, NtrkPassBL&geq;1,<br>p<sub>T</sub>>20 GeV<br>&mu;: lead p<sub>T</sub>>11 GeV | " " | N/A |
|12, ZBOSON_ASSIGNMENT     | nSFOS &ge; 1 (muon or track pairs). | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
| **Get jet container**  | Selected jets | JetHandler.cxx | N/A |
| **Overlap removal**    | Remove e if &Delta;R(e,&gamma;)<0.4<br> Remove jet if &Delta;R(j,&gamma;)<0.4<br>Remove jet if &Delta;R(j,e)<0.2<br>Remove e if &Delta;R(j,e)<0.4<br>Remove &mu; if &Delta;R(&mu;,&gamma;)<0.4<br>Remove &mu; if &Delta;R(&mu;,j)<0.4 | HGamAnalysisFramework/<br>OverlapRemovalHandler.cxx | " " |
|13, TWO_SF_LEPTONS_POSTOR | &ge; 1 of the SFOS pairs survives OR.<br>At this point, &mu;&mu; is preferred;<br>the remaining channels are exclusive. | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|14, BAD_MUON              | Reject events with presel BadMuons | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|15, ONE_PHOTON_POSTOR     | Selected photon survives OR | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|16, TRIG_MATCH            | Objects must match for a trigger<br>that fired (see config) | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|17, LEP_MEDID             | Resolved e: Medium,<br>Merged e: Merged **TMVA** ID,<br>&mu;: Medium | HGamAnalysisFramework/HGamRel21.config<br>for merged objects | " " |
|18, LEP_IP                | Res/Mrgd e: d<sub>0</sub>/&sigma;<sub>d0</sub> < 5, &#124;z<sub>0</sub>sin&theta;&#124; < 0.5<br>&mu;: d<sub>0</sub>/&sigma;<sub>d0</sub> < 3, &#124;z<sub>0</sub>sin&theta;&#124; < 0.5 | Defaults in ElectronHandler.cxx, MuonHandler.cxx | " " |
|19, LEP_ISO               | Resolved e: CloseByCorrected FCLoose<br>Merged e: FCLoose<br>&mu;: CloseByCorrected FCLoose_FixedRad | Merged e: specially done in<br>HiggsGamGamStarCutflowAndMxAOD.cxx.<br>Resolved e/&mu;: HGamAnalysisFramework/HGamRel21.config | " " |
|20, GAM_TIGHTID           | Photon passes Tight | HGamAnalysisFramework/HGamRel21.config | " " |
|21, GAM_ISOLATION         | Photon passes FixedCutLoose | HGamAnalysisFramework/HGamRel21.config | " " |
|22, ZMASSCUT              | m<sub>ll</sub> < 45 GeV<br> (see code for merged mass definiton) | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|23, LLGMASSCUT            | 105 < m<sub>ll&gamma;</sub> && m<sub>ll&gamma;</sub> < 160 GeV | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|24, LLMASSCUT             | Res/Mrgd e: m<sub>ll</sub> < 2.5 GeV or<br>m<sub>ll</sub> > 3.5 GeV (J/&Psi; peak veto)<br>&mu;: m<sub>ll</sub> < 2.9 GeV or<br>m<sub>ll</sub> > 3.3 GeV (J/&Psi; peak veto); <br> m<sub>ll</sub> < 9.1 GeV or<br>m<sub>ll</sub> > 10.6 GeV (&Upsilon; peak veto)    | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|25, DILEP_PT_FRAC         | p<sub>T</sub><sup>ll</sup>/m<sub>ll&gamma;</sub> > 0.3 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|26, GAM_PT_FRAC           | p<sub>T</sub><sup>&gamma;</sup>/m<sub>ll&gamma;</sub> > 0.3 | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |
|27, PASSALL               | Everything above passes | HiggsGamGamStarCutflowAndMxAOD.cxx | " " |

Notes on track selection above
-----------------------

The track selection has been harmonized to match the selection in the DAOD. This introduces an **"index requirement"** as
described below:

If the track is the **best-matched track** to an electron, it **automatically passes the index requirement** (this is for resolved electrons.)

Tracks considered for merged electrons must be **associated to a >5 GeV electron, have nSi+DeadSens &ge; 7 and passBLayerRequirement**.
Implicitly they have a p<sub>T</sub>>0.5 GeV (as they are associated to an electron) and &#124;&eta;&#124;<2.50.
If the track is the first (based on electron-track-matching ordering) to pass these requirements, then it
passes the index requirement, and it is assigned to `vtxTrkIndex1`. The next track (`electron->trackParticle(i)`) that (a) passes the selection, and
(b) is opposite-charge to the first track, is considered to have passed the index requirement, and assigned `vtxTrkIndex2`. No other tracks are accepted
(unless they are accepted due to their association with another track).

**Tracks that pass this index requirement are further required** to have p<sub>T</sub>>0.5 GeV (should be redundant), &#124;&eta;&#124;<2.50 (should be redundant), nSi+DeadSens &ge; 7 (redundant).

How to select events passing up to a certain cutflow point
-----------------------

Running on MxAODs, in order to provide a simple way to select events passing a certain stage of the event selection, one can require
e.g. (if one wants all events passing the trigger matching step):

    HGamEventInfoAuxDyn.cutFlow > TRIG_MATCH;

**Warning**: Be careful that if you use `16` instead of `TRIG_MATCH` you are in danger of your code becoming outdated with updates to the event selection!

Note also that the cutflow is **sequential**, so there is no way (using only the cutflow variable) to specify
e.g. something that passses photon ID but fails lepton ID (because they are not assessed in that order).

Convenience variables
-----------------------

For convenience, a few predefined booleans are also saved:

| Variable (HGamEventInfoAuxDyn.*) | Definition | Reason |
| -------- | ---------- | ------ |
| isPassedObjPreselection | m_cutFlow > TRIG_MATCH    | For bkg CRs - basic event selection passes, but objects are still "Loose" |
| isPassedObjSelection    | m_cutFlow > GAM_ISOLATION | All object selection passes at this point |
| isPassedEventSelection  | m_cutFlow >= PASSALL      | All event selection passes (including m<sub>ll</sub> and m<sub>ll&gamma;</sub> cuts) |
