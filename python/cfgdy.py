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

        #WR1600_N800
'file:/home/meier307/MINIAOD_21.root',
'file:/home/meier307/MINIAOD_22.root',
'file:/home/meier307/MINIAOD_23.root',
'file:/home/meier307/MINIAOD_24.root',
'file:/home/meier307/MINIAOD_25.root',
'file:/home/meier307/MINIAOD_26.root',
'file:/home/meier307/MINIAOD_27.root',
        #ttbar
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/079C3FC4-8835-394B-8E54-1C67DFAE7E8D.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0B4D0775-CC78-904D-A4B0-6B755608ABB5.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0D86B53B-2397-EE40-9A96-8115D6A754C2.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0F1A6FA9-1D64-3E4D-B6CE-69EC5C1462C9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1C576984-A8C6-B348-97FB-EEDC216ABDBD.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2498CC8E-233D-EE4B-91E1-467862BB453A.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/29C48882-452F-C14F-AFC9-C107D0623F83.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3AA76E9B-B61B-0647-957D-B698FA7C972A.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/52713434-77E4-EA4F-822F-3B3AE54C4E03.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/80E49E44-DFA3-C945-A509-521815C9808B.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/9E03FC31-C3F0-7641-B2CB-09E05C120375.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/A0218B26-3D0D-AD4B-8492-2D4192006B42.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/A75E8694-953F-FF43-B6C1-3DF07860CB4F.root',
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/C768B41A-863B-8C4F-85C8-D190C6B5A014.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/D49C6A22-12EA-3143-A46F-EECE13A1426D.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/D774DA06-04F6-5E45-B02E-70ECD0DD697F.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/DA0E374F-795A-BC44-8B7F-7E7BA3A69484.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/EF17AB4C-742D-A347-BCDB-D5C5E0EC2936.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/F77D5822-1E86-2D4F-BF44-1E7A483BC77B.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/006D6DFA-0A89-FB40-8E1B-B291EF3BB4A0.root',
#   'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0088149D-CE23-504B-8CD4-BFB55DABFDA5.root', 
 
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/00AB5775-0B55-4642-87B5-0588532E4A59.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0154A45F-B0B5-1747-BA1C-1F0534A317CA.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0159FCCB-9B6E-D241-A6EC-05A346036443.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/020DEED5-29A9-134B-B814-DFDD659379B5.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/02BC770B-EB0B-364F-AC14-3A5A8BB35A8B.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0323360A-9EB5-3943-9EF9-56FDD43D4962.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0377D569-1E63-8548-91DA-DB9DAC15BA85.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0406D9B7-8760-6D44-80CA-F08E39553F75.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/046A7C3C-3FB2-F340-931F-C49AB1342061.root',
 #  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/04B6756A-D517-3E49-831C-A378D2E13A95.root',

        #drell-yan
#  'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/00000/003F809B-022C-174A-ABFC-A5F8D27C101A.root',
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
        #WW
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
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/29425CF4-43B6-2040-95BC-A023A4B03F39.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2ADA69BB-A7D1-A44B-B7A3-FF8A6DD6E4CA.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2B5CEBDC-01BD-D24A-BF1E-B669C51B56E9.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/30A447EC-AF5E-ED43-864A-9286C5DAA270.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3133B703-B302-6A41-92F2-D546632164A8.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/36257AE3-6B47-0842-8CC6-D875E3D9CA04.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3F586030-F596-A34F-8138-D1E0E0A16ABE.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/444D63DC-4A96-494F-A5D3-CD6E137F027C.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/4672FEED-C9D8-E045-99B1-114C67466DA9.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/493D0BA5-13E2-2142-9137-D2BDBB89331E.root',
        #WZ
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/2FB6E6CF-3693-804E-933A-46B8DDEEA0D9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/31D1E0D2-53B2-1344-A43F-57F56701ED96.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/BAA059EC-D068-7D49-BD8D-A97A64FB7D34.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/E1E0486A-2629-A647-9980-7A19528FA7B0.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/100000/F74A2D91-9BB3-A24C-87F8-6B07C19E5E03.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/014F41AB-EA1B-A44C-A70D-8215AEE91342.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/016C85A5-CF13-514D-B024-D7911FEBF192.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/03C3EE48-411B-2E45-AD88-25DFCC3005A2.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/042FB9CE-9D7C-6F42-83E8-4BD3FA52F6A5.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/0444A3C3-EBC9-9E4B-ADDA-ED94CDE9BA24.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/04E45868-9F78-1F4A-A2C4-7B9CEC311189.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/0A74A225-7F01-CA44-86C6-ED95296216DC.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/0CD69F01-E403-8B44-833C-EFAA4F3EBEC9.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/183FA7F1-215E-C84D-8D87-B9EE52482B53.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/18BB6C29-461A-7F44-97ED-ED7CDA8DBFC8.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/1A44549F-9531-E140-AFC3-0682C926272B.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/1AB27B6F-8BA3-0A4A-BDFD-4FF386817F7D.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/1D395002-44DC-B946-A601-FCF565B7C815.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/28CF2D92-FF73-8447-9332-511568AF78C8.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/2A419830-0036-9B40-9C33-DAFEEF359166.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/2CC7B2AA-0813-DB4D-841B-F7ED61B52C8A.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/30233419-BBB0-404A-B52C-39D74A5749DA.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/3113091D-9D8D-1C44-B4BD-2CC41B022258.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/35C08F9D-3259-054D-AD10-6ACC4553A026.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/380F2BD1-FE27-1B4F-8701-99B915450966.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/431BB646-59D9-0445-B57E-057B1F0B23FB.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/458C040D-8259-5D40-A677-C833B821B331.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/478434CF-02F3-804B-97A4-99B9BD0303FB.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4A6838E6-176F-F548-9D8F-D66DBECEDC67.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4B4B8123-DE2D-EC43-AB8E-9042333C0479.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4B4E86FD-EDFB-B348-A4B7-9DEBF059A406.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4B78E0A8-D3A4-5941-A6F5-E7ED121D865D.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4C3F097F-711C-8743-A54D-86DCFCCE8623.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/4C5B9636-11DF-8A46-82C5-74BE15E47E44.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/52114A4F-1C58-894C-9E02-F37724E3CA25.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/539283EE-F265-DD47-863A-48FEC367401B.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/54D7F475-5E63-8746-ACAE-6E72B38FA52F.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/555B437C-36F5-214D-8596-2FD19BD7C155.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/WZ_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/60000/5681D495-A460-3E47-B5B7-931D016337BE.root',
        # ZZ
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/29C2DE0A-F874-474F-8DBF-9A5B99143AB9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/34ACA25E-0750-7C41-B1F3-9C2B6BB9AD7F.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/3A9094F2-88DF-1F4E-AFE6-5833BB8E1E92.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/3ECA98DC-3D54-8D41-A423-4A4BE9F0F720.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/42A63E59-CFF0-244E-894B-BE34BED02600.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/47876B1C-0E6D-5B41-8830-137588D309C1.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/4CD613E0-7A9E-8C4F-B0FD-9E2F4B5538E3.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/4F4A4ADF-106B-2849-917F-FD98A3F00AE0.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/520EC7AD-9B15-6A4A-B64D-01C87C40A4A0.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/5DE8DDCB-CD91-8F43-A117-ACA21C066965.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/5FDEA9D3-9B12-FE47-8E38-3FC2DDCF42D7.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/62510930-4BAD-0740-AC28-C25AD08524C4.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/632AEB28-7969-0948-9334-13A41B55F96B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/7C037D72-D14D-D543-8A88-9736FA86A233.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/80772B2D-4E8B-024E-BB43-56EBE1F656CD.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/8480B388-A24E-8D4C-982E-854A72559E81.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/854EED46-4E59-A442-87B5-8972B44CC1D0.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/89261A80-005C-A349-B11E-940C6922D441.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/8FED1DC0-32EA-CF44-B71A-DCEA5633306A.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/98E05B95-FE95-464E-89B5-01B8E09F884D.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/A45BDF38-BA60-7E4B-9CBC-E82EE139E2DA.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/A502A600-2DA9-1F44-A9A5-6348F9162ED5.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/A684104F-8334-9E4C-B14A-EA2BEFF72A5A.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/ABBB8509-E223-9247-80F5-0DC4ECC74FE1.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/B65A632D-415F-5542-9117-6980F52DFB89.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/CF5E487A-D91B-DA46-8DCB-558D573795A5.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/D23A7846-59A6-4543-8C25-12C4914C6B53.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/D2D82EBA-51B2-9E46-B85A-FEDE757BCEC1.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/ZZ_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/D5F494A9-BF67-754E-BC29-5DEBF7E14F81.root'



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

