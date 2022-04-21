import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi
from FWCore.ParameterSet.VarParsing import VarParsing
from PhysicsTools.PatAlgos.tools.jetTools import *
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
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_1.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_2.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_3.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_4.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_5.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_6.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_7.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_8.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_9.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_10.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_12.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_13.root', 
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_14.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_15.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_16.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_17.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_18.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_19.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_20.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_21.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_22.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_23.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_24.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_25.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_26.root',
# 'file:/local/cms/user/meier307/CMS/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/MINIAOD_27.root',

        #ttbar
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/079C3FC4-8835-394B-8E54-1C67DFAE7E8D.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0B4D0775-CC78-904D-A4B0-6B755608ABB5.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0D86B53B-2397-EE40-9A96-8115D6A754C2.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/0F1A6FA9-1D64-3E4D-B6CE-69EC5C1462C9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/1C576984-A8C6-B348-97FB-EEDC216ABDBD.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/2498CC8E-233D-EE4B-91E1-467862BB453A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/29C48882-452F-C14F-AFC9-C107D0623F83.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3AA76E9B-B61B-0647-957D-B698FA7C972A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/52713434-77E4-EA4F-822F-3B3AE54C4E03.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/80E49E44-DFA3-C945-A509-521815C9808B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/9E03FC31-C3F0-7641-B2CB-09E05C120375.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/A0218B26-3D0D-AD4B-8492-2D4192006B42.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/A75E8694-953F-FF43-B6C1-3DF07860CB4F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/C768B41A-863B-8C4F-85C8-D190C6B5A014.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/D49C6A22-12EA-3143-A46F-EECE13A1426D.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/D774DA06-04F6-5E45-B02E-70ECD0DD697F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/DA0E374F-795A-BC44-8B7F-7E7BA3A69484.root',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/EF17AB4C-742D-A347-BCDB-D5C5E0EC2936.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/F77D5822-1E86-2D4F-BF44-1E7A483BC77B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/006D6DFA-0A89-FB40-8E1B-B291EF3BB4A0.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0088149D-CE23-504B-8CD4-BFB55DABFDA5.root', 
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/00AB5775-0B55-4642-87B5-0588532E4A59.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0154A45F-B0B5-1747-BA1C-1F0534A317CA.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0159FCCB-9B6E-D241-A6EC-05A346036443.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/020DEED5-29A9-134B-B814-DFDD659379B5.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/02BC770B-EB0B-364F-AC14-3A5A8BB35A8B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0323360A-9EB5-3943-9EF9-56FDD43D4962.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0377D569-1E63-8548-91DA-DB9DAC15BA85.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0406D9B7-8760-6D44-80CA-F08E39553F75.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/046A7C3C-3FB2-F340-931F-C49AB1342061.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/04B6756A-D517-3E49-831C-A378D2E13A95.root',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/04B76B45-C745-0C44-B1EC-7A4E7BF3A5EF.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/06BA0964-D212-D440-9579-757CA881439A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/07AE5D91-A45E-8C4E-B071-262FE652EFF3.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/07B1BACA-E06A-764C-BC48-FCF7E70EC395.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0808FF87-B371-8A42-AA30-F052DCBF7969.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/08528C35-25BA-A349-B313-DFE941D172C7.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0869D9D7-F525-534E-AF02-812D67A52BC8.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/087C82DE-6E25-574F-8D49-EAAC16DB02DC.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/08DABABD-8369-B448-B328-E26144FEE196.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/09078C9E-41F2-8D46-B146-F760D41E7B8E.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/09595344-C123-174C-BCD1-893E17886C99.root',
  
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/09A39EEA-3CFE-A646-A00B-9BA4799BCEE9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/09B738E5-526E-3B4F-8B93-E834F4221B58.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0A7FC591-3BCF-7541-AA21-67C819F5D77A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0AA27EEC-66E4-6A40-A03E-D0E8498F8B90.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0AA67CD6-3AF3-3F4D-8CB5-AC7928A00326.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0AF2DF62-69C0-8848-A31C-14EAAABE8746.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0B3459D0-1242-5240-89EB-2D09B63F0B5B.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0B362470-DBEB-EC42-B3B7-D1CC149A979C.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0B9D5269-469E-7048-BC3D-2B6A5CBAE3E5.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0C01E2DC-E8A1-9E49-BEAE-8D0266D7B2EB.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0C103BD4-3F8F-AB41-B05E-73C17C114A65.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0C389E2A-4138-1C42-9455-32818E3F13C4.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0C5B50CC-3BC0-544A-8AFD-B47118D7E556.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0D461118-889C-9443-B239-3CB4A6B83154.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0DD8E98A-8AE3-F242-A63F-344BD5216527.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0DF4858D-DF5C-2644-9E53-CC6F94DE151C.root',

'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0DF8AA1A-32AF-054B-BB5D-630A22121EAB.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/0EAC6E13-EF61-4744-9B05-F4CF2687E553.root',
'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/104CE126-7C8C-4E4F-BF11-A1E704E4D5CB.root',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/105BDBB5-B780-BD46-82F2-B66CEA0FBB55.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/11381273-2BE2-EC41-9125-FC7EFD93A3F2.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/117C9220-9DA4-2C47-A7B0-00C288460587.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/124552A8-AC8E-484C-804A-0A3F060D0114.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1313FAA7-09C4-284F-8172-F894722F7788.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/13766D4D-5B9B-6449-A41A-CDCE839080CD.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1385437B-E07C-1E49-9305-BC9D54E085FB.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/13C4A89C-CE62-EF49-9757-505A6F58D79F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/14F8E041-7C82-C54D-A009-BD319674FFDB.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1543208D-7586-314F-A2D7-7F16F5889F46.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/15536481-5AAB-8B49-AFF7-294CC63F8DC9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1626CA12-EEFA-FC4C-82C7-9EE290417FAB.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/173167DD-345E-F547-8AFC-B58D26365B0B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/174B5A4D-E88A-914D-8132-34E04FDA4A21.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/17DBAC2F-F831-544C-B078-3C2B45BCA646.root',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/18FAF6FD-7AD5-FB4B-B1CA-ECD1D1D6CAF9.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/191980BA-EACF-3143-A5C2-750496852C53.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/193D2A76-3AEA-3B43-92BE-E7C8D802641A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/19C6F86F-EDDB-4E44-A4B3-F30347089B7F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1A2C3F23-192D-F749-B0EE-52E13735DC59.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1A8D988B-EB21-324E-A69B-9D065E7EF98A.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1AC73AD6-1756-E848-91A7-630026EA2C9C.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1B128757-6D4B-2F40-89D6-E67DED6E975D.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1BC32114-03AC-BF4D-8419-33E240576E18.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1BDBE4B3-B51E-574D-9914-FE339665FDB7.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1E19C2B8-2FEF-7640-B1D9-FC4622DB3E62.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1E734ECA-E36C-FB48-BD5D-FF7D19543605.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/1EF8FDBB-0C9D-5A4F-935B-AD96B7271790.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/20252306-6370-BB47-9C4E-FA64A1661DDB.root',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2151F69E-CF98-AB4A-A399-74D2FC9874DE.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2184C195-D721-9844-B772-2221FF4CA73B.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/21C1647B-91B4-C647-B33F-51E04673D16F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/222589E1-C385-2743-A899-9BB51BD796B0.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/22F1F190-B11B-3C4D-945F-C7D5327BD19F.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/230788AD-BBDE-6D45-8BFD-E24BB4D8E3DB.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/23604C2F-29F8-9541-A218-59956EE3E955.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/245F3B24-E798-EC4E-94D5-B9F04E6CBDB8.root',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/24621C12-76F9-CC4A-9C4B-244795A99623.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/24EB4F79-A1C0-9444-816F-C084B115C978.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2516EC22-167D-BE45-85FD-9669B587211B.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/25CD3FD9-BC9F-3746-84D2-2E2C0EFC19E9.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/25DE819B-E954-7742-BFCC-3DD5964A5D6E.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/26316D40-6BE3-A148-9667-FEC224DA23EE.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/26D83EA1-FC5A-AD4E-8E1A-F53DF8B6232A.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/26F46F80-1095-2546-9452-3943CFCDBAF4.root ',

# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/26F5BB5E-6DDE-3141-B3C5-BC466A1F8663.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2A5BB1B7-0919-6F45-AE35-B5DC52F00AB9.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2A73A57B-F2E4-DD41-8F83-D64F7B3FEE07.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2AAC3050-0853-E242-ABB6-1C30D92E2BCB.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2BF3039C-DE6C-7046-8AC0-94837E81B6FF.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2C2428A9-7121-AF4A-A9EC-B1F4EF04DDFC.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2C45CC5E-B615-204B-BFAA-8BCA73093196.root ',
# 'root://cmsxrootd.fnal.gov//store/mc/RunIIAutumn18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/120000/2C8B102B-87DE-954B-8F1B-AEEDE1058E05.root ',


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
                                   printPtEtaPhi = cms.untracked.bool(True),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(False),
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

# process.muonpPtFilter = cms.EDFilter("MCSingleParticleFilter",
#     MaxEta = cms.untracked.vdouble(2.4, 2.4),
#     Status = cms.untracked.vint32(1,  1),
#     MinEta = cms.untracked.vdouble(-2.4, -2.4),
#     MinPt = cms.untracked.vdouble(200, 200),
#     ParticleID = cms.untracked.vint32(13, -13)
# )

process.muonpPtFilter = cms.EDFilter("muonFilter",
      ptMin = cms.double(400),
      genParticles = cms.InputTag("prunedGenParticles"),
)

process.analysis = cms.EDAnalyzer('NNstudies',
                        tracks = cms.untracked.InputTag('ctfWithMaterialTracks'),
                        genParticles = cms.InputTag("prunedGenParticles"),
                        AK4genCHSJets = cms.InputTag("slimmedGenJets"),
                        highMuons = cms.InputTag("tuneIDMuons"),
                        highElectrons = cms.InputTag("heepElectrons"),
                        trigResults = cms.InputTag("TriggerResults","","HLT"),
                        genInfo = cms.InputTag("generator"),
                        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                        trainFile = cms.untracked.string(options.trainFile),
                        isSignal = cms.untracked.bool(options.isSignal),
                        genTrainData = cms.untracked.bool(options.genTrainData),
                        electronPathsToPass = electronPaths,
                        recoMET = cms.InputTag("slimmedMETs"),
                        packedGenParticles = cms.InputTag("packedGenParticles"),
                        packedPFCandidates = cms.InputTag("packedPFCandidates")
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
                           * process.muonSelectionSeq * process.analysis ) #* process.printTree)
# process.totalPath = cms.Path(process.selectedElectrons * process.heepSequence
#                            * process.analysis)# * process.printTree)
#process.totalPath = cms.Path(process.selectedElectrons * process.heepSequence
#                           * process.muonSelectionSeq * process.prefiringweight * process.analysis * process.printTree)

