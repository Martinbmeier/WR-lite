#ifndef eventBits2_h
#define eventBits2_h

#include <stdint.h>
#include <list>
#include <string>
#include <iostream>
#include <vector>

class eventBits2 {

public:

	eventBits2();
	void clear();


	bool twoElectrons;
	bool twoMuons;
	bool muonElectron;
	bool muonTau;
	bool electronTau;
	bool failedGenPtEta;
	bool hasPVertex;

	bool failedMatch;


	double leadLeptonPt;
	double leadLeptonType;
	double subleadLeptonPt;
	double subleadLeptonType;


	double quark1Eta;
	double quark1Phi;
	double quark1Mass;
	double quark1Pt;

	double quark2Eta;
	double quark2Phi;
	double quark2Mass;
	double quark2Pt;

	//EVENT WEIGHT
	double eventWeight;
	double count;

	//EVENT GEN VALUES
	bool electronTrigger;

	double lepton1Eta;
	double lepton1Phi;
	double lepton1Mass;
	double lepton1Pt;
	int lepton1Id;
	bool lepton1Matched;

	double lepton2Eta;
	double lepton2Phi;
	double lepton2Mass;
	double lepton2Pt;
	int lepton2Id;
	bool lepton2Matched;

	int muonGenCount;
	int electronGenCount;

	std::vector<double> unmatchedPhi = std::vector<double>(0);
	std::vector<double> unmatchedEta= std::vector<double>(0);
	std::vector<double> unmatchedPt= std::vector<double>(0);
	std::vector<double> unmatchedDR= std::vector<double>(0);

	//Electrons	
	double Electron1Eta;
	double Electron1Phi;
	double Electron1Pt;
	double Electron1dR;
	double Electron1PtRatio;
	bool Electron1TW;
	bool Electron1chargeMatch;

	
	double Electron2Eta;
	double Electron2Phi;
	double Electron2Pt;
	double Electron2dR;
	double Electron2PtRatio;
	bool Electron2TW;
	bool Electron2chargeMatch;

	int electronRecoCount;

	//Muons
	double Muon1Eta;
	double Muon1Phi;
	double Muon1Pt;
	double Muon1dR;
	double Muon1PtRatio;
	bool Muon1TW;
	bool Muon1chargeMatch;

	
	double Muon2Eta;
	double Muon2Phi;
	double Muon2Pt;
	double Muon2dR;
	double Muon2PtRatio;
	bool Muon2TW;
	bool Muon2chargeMatch;

	int muonRecoCount;

	double subElectronleadElectronRecodr2;
	double subMuonleadMuonRecodr2;


private:



};
#endif