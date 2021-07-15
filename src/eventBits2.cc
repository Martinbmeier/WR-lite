#include "ExoAnalysis/WR_lite/interface/eventBits2.h"

#include <list>
#include <string>


eventBits2::eventBits2() {
	//EVENT WEIGHT
	eventWeight = 0.0;
	count = 0.0;

	twoElectrons=false;
	twoMuons=false;
	muonElectron=false;
	muonTau=false;
	electronTau=false;
	hasPVertex=false;

	failedMatch=false;


	leadLeptonPt=-1000;
	leadLeptonType=-1000;
	subleadLeptonPt=-1000;
	subleadLeptonType=-1000;

	failedGenPtEta=false;

	quark1Eta=-1000;
	quark1Phi=-1000;
	quark1Mass=-1000;
	quark1Pt=-1000;

	quark2Eta=-1000;
	quark2Phi=-1000;
	quark2Mass=-1000;
	quark2Pt=-1000;

	//EVENT GEN VALUES
	electronTrigger=false;

	lepton1Eta=-1000;
	lepton1Phi=-1000;
	lepton1Mass=-1000;
	lepton1Pt=-1000;
	lepton1Id=0;

	lepton2Eta=-1000;
	lepton2Phi=-1000;
	lepton2Mass=-1000;
	lepton2Pt=-1000;
	lepton2Id=0;

	muonGenCount=0;
	electronGenCount=0;

	unmatchedPhi.clear();
	unmatchedEta.clear();
	unmatchedPt.clear();
	unmatchedDR.clear();

	//Electrons	
	Electron1Eta=-1000;
	Electron1Phi=-1000;
	Electron1Pt=-1000;
	Electron1dR=-1000;
	Electron1PtRatio=-1000;
	Electron1TW=false;
	Electron1chargeMatch=false;

	
	Electron2Eta=-1000;
	Electron2Phi=-1000;
	Electron2Pt=-1000;
	Electron2dR=-1000;
	Electron2PtRatio=-1000;
	Electron2TW=false;
	Electron2chargeMatch=false;

	electronRecoCount=0;

	//Muons
	Muon1Eta=-1000;
	Muon1Phi=-1000;
	Muon1Pt=-1000;
	Muon1dR=-1000;
	Muon1PtRatio=-1000;
	Muon1TW=false;
	Muon1chargeMatch=false;

	
	Muon2Eta=-1000;
	Muon2Phi=-1000;
	Muon2Pt=-1000;
	Muon2dR=-1000;
	Muon2PtRatio=-1000;
	Muon2TW=false;
	Muon2chargeMatch=false;


	muonRecoCount=0;

	subElectronleadElectronRecodr2=-1000;
	subMuonleadMuonRecodr2=-1000;

	
}

void eventBits2::clear() {
	//EVENT WEIGHT
	eventWeight = 0.0;
	count = 0.0;

	twoElectrons=false;
	twoMuons=false;
	muonElectron=false;
	muonTau=false;
	electronTau=false;
	hasPVertex=false;

	failedMatch=false;

	leadLeptonPt=-1000;
	leadLeptonType=-1000;
	subleadLeptonPt=-1000;
	subleadLeptonType=-1000;

	failedGenPtEta=false;

	quark1Eta=-1000;
	quark1Phi=-1000;
	quark1Mass=-1000;
	quark1Pt=-1000;

	quark2Eta=-1000;
	quark2Phi=-1000;
	quark2Mass=-1000;
	quark2Pt=-1000;


	//EVENT GEN VALUES
	electronTrigger=false;

	lepton1Eta=-1000;
	lepton1Phi=-1000;
	lepton1Mass=-1000;
	lepton1Pt=-1000;
	lepton1Id=0;

	lepton2Eta=-1000;
	lepton2Phi=-1000;
	lepton2Mass=-1000;
	lepton2Pt=-1000;
	lepton2Id=0;

	muonGenCount=0;
	electronGenCount=0;

	unmatchedPhi.clear();
	unmatchedEta.clear();
	unmatchedPt.clear();
	unmatchedDR.clear();

	//Electrons	
	Electron1Eta=-1000;
	Electron1Phi=-1000;
	Electron1Pt=-1000;
	Electron1dR=-1000;
	Electron1PtRatio=-1000;
	Electron1TW=false;
	Electron1chargeMatch=false;

	
	Electron2Eta=-1000;
	Electron2Phi=-1000;
	Electron2Pt=-1000;
	Electron2dR=-1000;
	Electron2PtRatio=-1000;
	Electron2TW=false;
	Electron2chargeMatch=false;

	electronRecoCount=0;

	//Muons
	Muon1Eta=-1000;
	Muon1Phi=-1000;
	Muon1Pt=-1000;
	Muon1dR=-1000;
	Muon1PtRatio=-1000;
	Muon1TW=false;
	Muon1chargeMatch=false;

	
	Muon2Eta=-1000;
	Muon2Phi=-1000;
	Muon2Pt=-1000;
	Muon2dR=-1000;
	Muon2PtRatio=-1000;
	Muon2TW=false;
	Muon2chargeMatch=false;

	muonRecoCount=0;

	subElectronleadElectronRecodr2=-1000;
	subMuonleadMuonRecodr2=-1000;
	
}


	
