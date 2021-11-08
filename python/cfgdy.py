import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi
from FWCore.ParameterSet.VarParsing import VarParsing
import os
#INPUT PARSING SECTION
options = VarParsing ('analysis')

options.register( 'isSignal',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "True when running over signal MC samples"
               )

options.register( 'nEvents',
                  -1,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "Number of events. Default is all events (-1)"
               )

options.register( 'genTrainData',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "True when generating a training set for neural network"
                   )       
           
options.register( 'trainFile',
                  'nn.txt',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name of file for training set for neural network"
                   )
options.register ('gitTag',
                  'applePie',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Git Tag")

options.register( 'isMC',
                  True,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "True if is MC dataset")

options.register( 'ISmcatnlo',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "True if is MC @ NLO dataset")

options.register( 'doFast',
                  False,
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "True to run a condensed analysis for HiggsCombine"
               )

options.register( 'reco',
                  '80X',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "MC or data CMSSW reco")

options.register( 'era',
                  '2016',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Year of Run II")

options.register( 'checkZ',
          True,
          VarParsing.multiplicity.singleton,
          VarParsing.varType.bool,
          "True when running over Drell-Yan MC samples"
           )
options.parseArguments()

#LOCAL VARIABLE DEFINITIONS
muonID =' userInt("highPtID") == 1'

process = cms.Process("Analysis")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#setup global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
from Configuration.AlCa.autoCond import autoCond

process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3') #
if not options.isMC: process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvents) )
process.source = cms.Source ("PoolSource",
      fileNames = cms.untracked.vstring (
#         'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/003F809B-022C-174A-ABFC-A5F8D27C101A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/0073D522-807B-3744-8A19-6FA4077AA7CE.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/01376E26-2F41-1145-834F-DCBC2358E3EC.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/02D18F78-E048-7646-8F8E-E714BA38DC3E.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/04834C44-478D-CC4F-BC64-56FC081E3BCB.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/05EEC364-4D12-DB47-BAF6-18CC75A0A00B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/08ED2FF0-4345-024A-B1C8-A7F412D2E4FC.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/08FCC054-81BC-2C4F-AF2D-16D2B9D53164.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/0A75830D-A212-294D-9DE1-25359228C4E7.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/128129FD-30DB-FE4A-9B43-C383A61A492F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/136E3643-3C2F-524A-B2BC-E764624D65AE.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/1420AC0A-AB2F-8F4B-A5FC-8E040B670492.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/172F250E-8FDC-1C45-93C1-79ABE8E343B4.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/1F1FD277-CF2A-2844-8D34-C6505323E8FA.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/21D6AB94-83FF-6A4C-96DC-6502880C52CB.root'
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/001109DA-8045-D042-B06E-3D532A77474B.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0426496B-695F-3644-8F91-D11095F1308B.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0CCABC08-81D1-9B4E-8DF3-8D0ACFBBC1D0.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0FAD41E4-9E0D-2546-A139-394DB38288D6.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1098D0E1-44F1-8948-A251-E29014065F9C.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/11AE3311-392A-F34E-95AB-28044688B17F.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1440DF45-6915-6443-8CF8-B9830DD9BD30.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/161BBFBA-6308-9C42-99BB-048261E227D7.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/16E07E8E-B3C4-9544-AACB-575E2DDBC25F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/18DFDC70-4BAE-6242-824F-AE7991634BB4.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1B69C17B-2AB4-2449-9866-D13A60C12717.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1CEA7BD3-9104-264E-99EB-F2529EA4382E.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1E05A44D-A30C-2B4C-B3D9-0429BBEFCA29.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1E67C683-E5A1-9D46-8201-EF5BB8D45BE9.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2111E83D-E509-2643-9511-A36C506EA331.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/238A560F-E6CF-2D4C-82F4-86B3FE7B9926.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2721CD77-2441-BF43-9F5A-4957317104E1.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/29425CF4-43B6-2040-95BC-A023A4B03F39.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2ADA69BB-A7D1-A44B-B7A3-FF8A6DD6E4CA.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2B5CEBDC-01BD-D24A-BF1E-B669C51B56E9.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/30A447EC-AF5E-ED43-864A-9286C5DAA270.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3133B703-B302-6A41-92F2-D546632164A8.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/36257AE3-6B47-0842-8CC6-D875E3D9CA04.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3F586030-F596-A34F-8138-D1E0E0A16ABE.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/444D63DC-4A96-494F-A5D3-CD6E137F027C.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/4672FEED-C9D8-E045-99B1-114C67466DA9.root', 
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/493D0BA5-13E2-2142-9137-D2BDBB89331E.root' 
) 
)

#process.source = cms.Source ("PoolSource",
#     fileNames = cms.untracked.vstring ('file:/hdfs/cms/user/evans908/wrSkims/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/skimmedOut_1.root')
#)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
process.options.allowUnscheduled = cms.untracked.bool(False)

process.TFileService = cms.Service("TFileService", 
                        fileName = cms.string(options.outputFile)
)  

process.badGlobalMuonTagger = cms.EDFilter("BadGlobalMuonTagger",
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonPtCut = cms.double(20),
    selectClones = cms.bool(False),
    taggingMode = cms.bool(True),
    verbose     = cms.untracked.bool(False)
)
process.cloneGlobalMuonTagger = process.badGlobalMuonTagger.clone(
    selectClones = cms.bool(True)
)

process.removeBadAndCloneGlobalMuons = cms.EDProducer("MuonRefPruner",
    input = cms.InputTag("slimmedMuons"),
    toremove = cms.InputTag("badGlobalMuonTagger", "bad"),
    toremove2 = cms.InputTag("cloneGlobalMuonTagger", "bad")
)

process.tunePMuons = cms.EDProducer("TunePMuonProducer",
        src = cms.InputTag("removeBadAndCloneGlobalMuons")
        #src = cms.InputTag("slimmedMuons")
)

### muon ID and isolation
# make a collection of TuneP muons which pass isHighPt ID
process.tuneIDMuons = cms.EDFilter("PATMuonSelector",
                               src = cms.InputTag("tunePMuons"),
                               cut = cms.string(muonID),
)
#HERE WE RUN A MODULE FROM SAM HARPER WHICH INSERTS HEEP CUT INFO INTO THE PAT ELECTRON USER DATA
#we setup the HEEP ID V7.0 and enable VID via the following function
#and then add it to a new collection of pat::Electrons
#there is the option to call the new collection "slimmedElectrons" (useStdName=True)
#otherwise it calls them "heepElectrons"
#it creates a sequence "process.heepSequence" which we add to our path

#from HEEP.VID.tools import addHEEPV70ElesMiniAOD

from python.tools import addHEEPV70ElesMiniAOD
addHEEPV70ElesMiniAOD(process,useStdName=False)



process.options.allowUnscheduled = cms.untracked.bool(True)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),    
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(True),
                                   #status = cms.untracked.vint32( 3 )
                                   )

#process.p = cms.Path(
#    process.heepSequence*
#    process.heepIdExample) #our analysing example module, replace with your module
#

process.muonSelectionSeq = cms.Sequence(cms.ignore(process.badGlobalMuonTagger) * cms.ignore(process.cloneGlobalMuonTagger) * process.removeBadAndCloneGlobalMuons * process.tunePMuons * process.tuneIDMuons)

muonPaths = cms.vstring("HLT_Mu50_v", "HLT_TkMu50_v")
electronPaths = cms.vstring("HLT_Ele27_WPTight_Gsf_v", "HLT_Photon175_v")
if options.era == '2016':
    muonPaths = cms.vstring("HLT_Mu50_v", "HLT_TkMu50_v")
    electronPaths = cms.vstring("HLT_Ele27_WPTight_Gsf_v", "HLT_Ele115_CaloIdVT_GsfTrkIdT_v", "HLT_Photon175_v")
elif options.era == '2017':
    muonPaths = cms.vstring("HLT_Mu50_v", "HLT_OldMu100", "HLT_TkMu100")
    electronPaths = cms.vstring("HLT_Ele35_WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v")
elif options.era == '2018':
    muonPaths = cms.vstring("HLT_Mu50", "HLT_OldMu100", "HLT_TkMu100")
    electronPaths = cms.vstring("HLT_Ele32WPTight_Gsf_v","HLT_Photon200_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v")

#process.analysis = cms.EDAnalyzer('WR_MASS_PLOT',
#                       tracks = cms.untracked.InputTag('ctfWithMaterialTracks'),
#                       genParticles = cms.InputTag("prunedGenParticles"),
#                       AK4recoCHSJets = cms.InputTag("slimmedJets"),
#                       highMuons = cms.InputTag("tuneIDMuons"),
#                       highElectrons = cms.InputTag("heepElectrons"),
#                       genInfo = cms.InputTag("generator"),
#                       vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                       trainFile = cms.untracked.string(options.trainFile),
#                       isSignal = cms.untracked.bool(options.isSignal),
#                       genTrainData = cms.untracked.bool(options.genTrainData)
#)

process.analysis = cms.EDAnalyzer('cut_flow2',
                        tracks = cms.untracked.InputTag('ctfWithMaterialTracks'),
                        genParticles = cms.InputTag("prunedGenParticles"),
                        AK4recoCHSJets = cms.InputTag("slimmedJets"),
                        highMuons = cms.InputTag("tuneIDMuons"),
                        highElectrons = cms.InputTag("heepElectrons"),
                        trigResults = cms.InputTag("TriggerResults","","HLT"),
                        genInfo = cms.InputTag("generator"),
                        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                        trainFile = cms.untracked.string(options.trainFile),
                        isSignal = cms.untracked.bool(options.isSignal),
                        genTrainData = cms.untracked.bool(options.genTrainData),
                        electronPathsToPass = electronPaths,
                        recoMET = cms.InputTag("slimmedMETs")

)


process.load('FWCore.Modules.printContent_cfi')

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    
process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)")
)
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
process.heepIDVarValueMaps.elesMiniAOD  = 'selectedElectrons'
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.heepElectrons.src = cms.InputTag('selectedElectrons')

####EE L1 Prefiring Correction ####
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
if options.era == '2016':
  process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
      DataEra = cms.string("2016BtoH"),
      UseJetEMPt = cms.bool(True),
      PrefiringRateSystematicUncty = cms.double(0.2),
      SkipWarnings = False)
elif options.era == '2017':
  process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
      DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
      UseJetEMPt = cms.bool(True),
      PrefiringRateSystematicUncty = cms.double(0.2),
      SkipWarnings = False)
elif options.era == '2018':
  process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
      DataEra = cms.string("2018AtoD"), #Use 2016BtoH for 2016
      UseJetEMPt = cms.bool(True),
      PrefiringRateSystematicUncty = cms.double(0.2),
      SkipWarnings = False)


process.totalPath = cms.Path(process.selectedElectrons * process.heepSequence
                           * process.muonSelectionSeq * process.analysis)# * process.printTree)
#process.totalPath = cms.Path(process.selectedElectrons * process.heepSequence
#                           * process.muonSelectionSeq * process.prefiringweight * process.analysis * process.printTree)

