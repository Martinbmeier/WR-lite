
// -*- C++ -*-
//
// Package:    ExoAnalysis/gen_match
// Class:      gen_match
//
/**\class gen_match gen_match.cc ExoAnalysis/gen_match/plugins/gen_match.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Evans
//         Created:  Mon, 14 Oct 2019 19:43:16 GMT
//
//

// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "HEEP/VID/interface/CutNrs.h"
#include "HEEP/VID/interface/VIDCutCodes.h"

#include "TLorentzVector.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <cmath>

#include "TH1.h"
#include "TDirectory.h"


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/eventInfo.h"
#include "ExoAnalysis/WR_lite/interface/eventHistos2.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class gen_match : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit gen_match(const edm::ParameterSet&);
		~gen_match();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		double dR2(double eta1, double eta2, double phi1, double phi2);
		double dPhi(double phi1, double phi2);
		bool tWfinder(const edm::Event&, const reco::GenParticle* );
		bool passElectronTrig(const edm::Event& , eventBits2& );
		//double transverseSphericity(math::XYZTLorentzVector p1, math::XYZTLorentzVector p2, math::XYZTLorentzVector p3);
		//void saveElectronData(eventBits2 * iCut[ii], double matched1Mass, double matched2Mass);
		//void saveMuonData(eventBits2 * iCut[ii], double matched1Mass, double matched2Mass);
		
		
		eventHistos2 m_histoMaker;


		//neuralNet networkResolved = neuralNet("/home/kronh006/Version3/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/data/Resolved");
		//neuralNet networkSuperResolved = neuralNet("/home/kronh006/Version3/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/data/SuperResolved");

		// ----------member data ---------------------------

		edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
		edm::EDGetToken m_genParticleToken;
		edm::EDGetToken m_highMuonToken;
		edm::EDGetToken m_highElectronToken;
		edm::EDGetToken m_AK4recoCHSJetsToken;
		edm::EDGetToken m_genEventInfoToken;
		edm::EDGetToken m_offlineVerticesToken;
		std::vector<std::string>  m_electronPathsToPass;
		edm::EDGetToken m_trigResultsToken;
		
		std::string m_dataSaveFile;
		bool m_isSignal;
		bool m_genTrainData;


		std::string  cSV_bTag1      = "pfDeepCSVJetTags:probb";
		std::string  cSV_bTag2      = "pfDeepCSVJetTags:probbb";

		
		edm::Service<TFileService> fs; 
		
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
gen_match::gen_match(const edm::ParameterSet& iConfig)
	:
	tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	m_genParticleToken(consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"))),
	m_highMuonToken (consumes<std::vector<pat::Muon>> (iConfig.getParameter<edm::InputTag>("highMuons"))),
	m_highElectronToken (consumes<std::vector<pat::Electron>> (iConfig.getParameter<edm::InputTag>("highElectrons"))),
	m_AK4recoCHSJetsToken (consumes<std::vector<pat::Jet>> (iConfig.getParameter<edm::InputTag>("AK4recoCHSJets"))),
	m_genEventInfoToken (consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"))),
	m_offlineVerticesToken (consumes<std::vector<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("vertices"))),
	m_dataSaveFile (iConfig.getUntrackedParameter<std::string>("trainFile")),
	m_isSignal (iConfig.getUntrackedParameter<bool>("isSignal"))
	//m_genTrainData (iConfig.getUntrackedParameter<bool>("genTrainData"))

{
   //now do what ever initialization is needed

   m_electronPathsToPass  = iConfig.getParameter<std::vector<std::string> >("electronPathsToPass");
   m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
}


gen_match::~gen_match()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

// ------------ method called for each event  ------------
void
gen_match::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	  eventBits2 cuts1; //pt > 10 , eta < 2.4
	  eventBits2 cuts2; //pt > 10 , eta < 2.4 , isHighPtMuon
	  eventBits2 cuts3; //pt > 10 , eta < 2.4 , isolationR03
	  eventBits2 cuts4; //pt > 10 , eta < 2.4 , isHighPtMuon, isolationR03

	  std::vector<eventBits2> iCut; 
	  iCut.push_back(cuts1);
	  iCut.push_back(cuts2);
	  iCut.push_back(cuts3);
	  iCut.push_back(cuts4);

   
   eventInfo myEvent; 
	
	edm::Handle<GenEventInfoProduct> eventInfo;
	iEvent.getByToken(m_genEventInfoToken, eventInfo);

	edm::Handle<std::vector<pat::Muon>> highMuons;
	iEvent.getByToken(m_highMuonToken, highMuons);

	edm::Handle<std::vector<pat::Electron>> highElectrons;
	iEvent.getByToken(m_highElectronToken, highElectrons);

	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByToken(m_genParticleToken, genParticles);

	std::vector<pat::Muon> prunedMuons;

for(int ii=0; ii<4; ii++){

	prunedMuons.clear();

	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){
		if(ii==0){
			if(fabs(iMuon->eta()) < 2.4 && iMuon->tunePMuonBestTrack()->pt() > 10 ){ prunedMuons.push_back(*iMuon);}
		}
		if(ii==1){
			if(fabs(iMuon->eta()) < 2.4 && iMuon->tunePMuonBestTrack()->pt() > 10 && (iMuon->isHighPtMuon(*myEvent.PVertex))){ prunedMuons.push_back(*iMuon);}
		}
		if(ii==2){
			if(fabs(iMuon->eta()) < 2.4 && iMuon->tunePMuonBestTrack()->pt() > 10 && (iMuon->isolationR03().sumPt/iMuon->pt() < .1)){ prunedMuons.push_back(*iMuon);}
		}
		if(ii==3){
			if(fabs(iMuon->eta()) < 2.4 && iMuon->tunePMuonBestTrack()->pt() > 10 && (iMuon->isHighPtMuon(*myEvent.PVertex)) && (iMuon->isolationR03().sumPt/iMuon->pt() < .1)){ prunedMuons.push_back(*iMuon);}
		}
		//const pat::Muon goodMuon = *iMuon;
		//prunedMuons.push_back(*iMuon);
	}
  
	iCut[ii].count = eventInfo->weight()/fabs(eventInfo->weight());
	iCut[ii].eventWeight = eventInfo->weight();
	
	
	edm::Handle<std::vector<reco::Vertex>> vertices;
	iEvent.getByToken(m_offlineVerticesToken, vertices);
	if(!vertices.isValid()) {
		throw cms::Exception("Vertex collection not valid!");
	}

	iCut[ii].hasPVertex = myEvent.PVselection(vertices);


// get gen information 
		int lepton1Cuts = 0;
		int lepton2Cuts = 0;
		int lepton1ID = 0;
		int lepton2ID = 0;

		const reco::GenParticle* lepton1   = 0;
		const reco::GenParticle* lepton2   = 0;

		for (std::vector<reco::GenParticle>::const_iterator iParticle = genParticles->begin(); iParticle != genParticles->end(); iParticle++) {
			if( ! iParticle->isHardProcess() ){ continue; }  //ONLY HARD PROCESS AND NOT INCOMING
			if( abs( iParticle->pdgId() ) == 13 || abs( iParticle->pdgId() ) == 11 || abs( iParticle->pdgId() ) == 15) {//HERE'S A LEPtON
				if(fabs(iParticle->eta()) > 2.4 || iParticle->pt() < 10){
					if(lepton1Cuts==0){
						lepton1Cuts=1;
					} else if(lepton2Cuts==0){
						lepton2Cuts=1;
					}
					iCut[ii].failedGenPtEta = true;
					continue;
				}

				if(abs(iParticle->pdgId())==13){iCut[ii].muonGenCount++;}
				if(abs(iParticle->pdgId())==11){iCut[ii].electronGenCount++;}
				if(!tWfinder(iEvent, &(*iParticle))){ continue; }		//Only final state muons and electrons from t->W
				if(lepton1==0){
					lepton1 = &(*iParticle);
					lepton1ID = abs(iParticle->pdgId());
				} else if(lepton2==0){
					lepton2 = &(*iParticle);
					lepton2ID = abs(iParticle->pdgId());
				}
			}
		}

		if ((lepton1 == 0) || (lepton2 == 0)) {
			iCut[ii].failedGenPtEta = true; 
		} else if((lepton1ID == 11 && lepton2ID == 13) || (lepton1ID == 13 && lepton2ID == 11 )) {
			iCut[ii].muonElectron = true;
		} else if(lepton1ID == 11 && lepton2ID == 11){
			iCut[ii].twoElectrons = true;
		} else if(lepton1ID == 13 && lepton2ID == 13){
			iCut[ii].twoMuons = true;
		} else if((lepton1ID == 11 && lepton2ID == 15) || (lepton1ID == 15 && lepton2ID == 11)) {
			iCut[ii].electronTau = true;
		} else if((lepton1ID == 13 && lepton2ID == 15 ) || ( lepton1ID == 15 && lepton2ID == 13)) {
			iCut[ii].muonTau = true;
		}

		//gen values
		if(!iCut[ii].failedGenPtEta){
		iCut[ii].lepton1Eta = lepton1->eta();
		iCut[ii].lepton1Phi = lepton1->phi();
		iCut[ii].lepton1Mass = lepton1->p4().mass();
		iCut[ii].lepton1Pt = lepton1->pt();
		iCut[ii].lepton1Id = lepton1->pdgId();
		
		iCut[ii].lepton2Eta = lepton2->eta();
		iCut[ii].lepton2Phi = lepton2->phi();
		iCut[ii].lepton2Mass = lepton2->p4().mass();
		iCut[ii].lepton2Pt = lepton2->pt();
		iCut[ii].lepton2Id = lepton2->pdgId();
		}



//start reco

	edm::Handle<std::vector<pat::Jet>> recoJetsAK4;  
	iEvent.getByToken(m_AK4recoCHSJetsToken, recoJetsAK4);  


	//const pat::Jet* leadJet = 0;
	//const pat::Jet* subleadJet = 0;
	int jetCount = 0;

	//bool twoBjets = false;


	//Get jets with maximum pt
		for(std::vector<pat::Jet>::const_iterator iJet = recoJetsAK4->begin(); iJet != recoJetsAK4->end(); iJet++) {
			//Make sure jets are not around leptons
			//if ( fabs(iJet->eta()) > 2.4) continue;
			double NHF  =           iJet->neutralHadronEnergyFraction();
			double NEMF =           iJet->neutralEmEnergyFraction();
			double CHF  =           iJet->chargedHadronEnergyFraction();
			double CEMF =           iJet->chargedEmEnergyFraction();
			double NumConst =       iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
			double MUF      =       iJet->muonEnergyFraction();
			double EUF      =       iJet->electronEnergyFraction();
			double CHM      =       iJet->chargedMultiplicity(); 
			double BJP		 =       iJet->bDiscriminator(cSV_bTag1) + iJet->bDiscriminator(cSV_bTag2); 
			//APPLYING TIGHT QUALITY CUTS
			if (NHF > .9) continue;
			if (NEMF > .9) continue;
			if (NumConst <= 1) continue;
			if (MUF >= .8) continue; //MAKE SURE THE AREN'T MUONS
			if (EUF >= .8) continue; //MAKE SURE THE AREN'T ELECTRONS
			//ADDITIONAL CUTS BECAUSE OF TIGHT ETA CUT
			if (CHF == 0) continue;
			if (CHM == 0) continue;
			//if (CEMF > .99) continue;
			if (CEMF > .90)  continue; 
			if (BJP < 0.4184) continue;
			
			if (jetCount == 0) {
				//leadJet = &(*(iJet));
			} else if (jetCount == 1) {
				//subleadJet = &(*(iJet));
			}
			jetCount++;
		}
		if(jetCount==2){//twoBjets=true;
			//m_eventType->Fill("2 tagged b jets",1);
		}

	int mu1Match = 0;
	int mu2Match = 0;
	int el1Match = 0;
	//int el2Match = 0;

	const pat::Muon* matchedMuon1;
	const pat::Muon* matchedMuon2;

	const pat::Electron* matchedElectron1;
	//const pat::Electron* matchedElectron2;



//muon reconstruction

if(iCut[ii].twoMuons){

	double match1DR = 1000;
	double match2DR = 1000;

   for(std::vector<pat::Muon>::const_iterator iMuon = prunedMuons.begin(); iMuon != prunedMuons.end(); iMuon++){


   	//if(fabs(iMuon->eta()) > 2.4 || iMuon->tunePMuonBestTrack()->pt() < 10 || !(iMuon->isHighPtMuon(*myEvent.PVertex)) || (iMuon->isolationR03().sumPt/iMuon->pt() > .1)){ continue;}// || !(iMuon->isHighPtMuon(*myEvent.PVertex)) || (iMuon->isolationR03().sumPt/iMuon->pt() > .1)){ continue;} //preliminary cu

   	iCut[ii].muonRecoCount++;

		match1DR = sqrt(dR2(iMuon->eta(), iCut[ii].lepton1Eta, iMuon->phi(), iCut[ii].lepton1Phi));
		match2DR = sqrt(dR2(iMuon->eta(), iCut[ii].lepton2Eta, iMuon->phi(), iCut[ii].lepton2Phi));
		//iCut[ii].recoMuonCount++;

		if(mu1Match == 0 && match1DR < 0.3 && match1DR < match2DR) {
			mu1Match++;
			iCut[ii].lepton1Matched=true;
			matchedMuon1 = &(*(iMuon));
			iCut[ii].Muon1dR=match1DR;
			if(matchedMuon1->pdgId()==iCut[ii].lepton1Id){ iCut[ii].Muon1chargeMatch = true ;}
    		else{ iCut[ii].Muon1chargeMatch = false ;}
			iCut[ii].Muon1PtRatio = iCut[ii].lepton1Pt/matchedMuon1->pt();
			//iCut[ii].Muon1TW=tWfinder(iEvent,lepton1);
			iCut[ii].Muon1Pt=matchedMuon1->pt();
			iCut[ii].Muon1Eta=matchedMuon1->eta();
			iCut[ii].Muon1Phi=matchedMuon1->phi();
		}
		if(mu2Match == 0 && match2DR < 0.3 && match2DR < match1DR)  {
			matchedMuon2 = &(*(iMuon));
			mu2Match++; 
			iCut[ii].lepton2Matched=true;
			iCut[ii].Muon2dR=match2DR;
			if(matchedMuon2->pdgId()==iCut[ii].lepton2Id){ iCut[ii].Muon2chargeMatch = true ;}
    		else{ iCut[ii].Muon2chargeMatch = false ;}
			iCut[ii].Muon2PtRatio = iCut[ii].lepton2Pt/matchedMuon2->pt();
			//iCut[ii].Muon2TW=tWfinder(iEvent,lepton2);
			iCut[ii].Muon2Pt=matchedMuon2->pt();
			iCut[ii].Muon2Eta=matchedMuon2->eta();
			iCut[ii].Muon2Phi=matchedMuon2->phi();
		}
		if(mu2Match==0 && mu1Match==0){		
			iCut[ii].unmatchedPhi.push_back(iMuon->phi());
			iCut[ii].unmatchedEta.push_back(iMuon->eta());
			iCut[ii].unmatchedPt.push_back(iMuon->pt());
			if(match1DR<match2DR){
				iCut[ii].unmatchedDR.push_back(match1DR);
			}
			else{iCut[ii].unmatchedDR.push_back(match2DR);}
	}
	}
	if(mu1Match==0 || mu2Match==0 || iCut[ii].failedGenPtEta ){
		iCut[ii].failedMatch=true;
	}

   
}

/*
//electron reconstruction
if(iCut[ii].twoElectrons){

	iCut[ii].electronTrigger=false;
   if (passElectronTrig(iEvent, iCut[ii])){ iCut[ii].electronTrigger=true;}

   double match1DR = 1000;
	double match2DR = 1000;

	//for all reco electrons, loop through gen electrons to find spatial matches
	for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++) {
		if(fabs(iElectron->eta()) > 2.4) {continue;}
		if(iElectron->pt() < 10 ) {continue;}

		iCut[ii].electronRecoCount++;

    	match1DR = sqrt(dR2(iElectron->eta(), iCut[ii].lepton1Eta, iElectron->phi(), iCut[ii].lepton1Phi));
		match2DR = sqrt(dR2(iElectron->eta(), iCut[ii].lepton2Eta, iElectron->phi(), iCut[ii].lepton2Phi));

		if(el1Match == 0 && match1DR < 0.3 && match1DR<match2DR ) {
			el1Match++;
			iCut[ii].lepton1Matched=true;
			matchedElectron1 = &(*(iElectron));
			iCut[ii].Electron1dR=match1DR;
			if(matchedElectron1->pdgId()==iCut[ii].lepton1Id){ iCut[ii].Electron1chargeMatch = true ;}
    		else{ iCut[ii].Electron1chargeMatch = false ;}
			iCut[ii].Electron1PtRatio = iCut[ii].lepton1Pt/matchedElectron1->pt();
			iCut[ii].Electron1TW=tWfinder(iEvent,lepton1);
			iCut[ii].Electron1Pt=matchedElectron1->pt();
			iCut[ii].Electron1Eta=matchedElectron1->eta();
			iCut[ii].Electron1Phi=matchedElectron1->phi();
		}
		if(el2Match == 0 && match2DR < 0.3 && match2DR<match1DR)  {
			matchedElectron2 = &(*(iElectron));
			el2Match++; 
			iCut[ii].lepton2Matched=true;
			iCut[ii].Electron2dR=match2DR;
			if(matchedElectron2->pdgId()==iCut[ii].lepton2Id){ iCut[ii].Electron2chargeMatch = true ;}
    		else{ iCut[ii].Electron2chargeMatch = false ;}
			iCut[ii].Electron2PtRatio = iCut[ii].lepton2Pt/matchedElectron2->pt();
			iCut[ii].Electron2TW=tWfinder(iEvent,lepton2);
			iCut[ii].Electron2Pt=matchedElectron2->pt();
			iCut[ii].Electron2Eta=matchedElectron2->eta();
			iCut[ii].Electron2Phi=matchedElectron2->phi();
		}
		if(el1Match==0 && el2Match==0){		
			iCut[ii].unmatchedPhi.push_back(iElectron->phi());
			iCut[ii].unmatchedEta.push_back(iElectron->eta());
			iCut[ii].unmatchedPt.push_back(iElectron->pt());
			if(match1DR<match2DR){
				iCut[ii].unmatchedDR.push_back(match1DR);
			}
			else{iCut[ii].unmatchedDR.push_back(match2DR);}
		}

	}

	if(el1Match==0 || el2Match==0 || iCut[ii].failedGenPtEta ){
	iCut[ii].failedMatch=true;
	}

}
*/

//muon/electron reconstruction

//muon
if(iCut[ii].muonElectron  || iCut[ii].muonTau) {  //|| iCut[ii].electronTau ) {

	double match1DR = 1000;
	double match2DR = 1000;

	if(!iCut[ii].electronTau){
   for(std::vector<pat::Muon>::const_iterator iMuon = prunedMuons.begin(); iMuon != prunedMuons.end(); iMuon++){

   	//if(fabs(iMuon->eta()) > 2.4 || iMuon->tunePMuonBestTrack()->pt() < 10 || !(iMuon->isHighPtMuon(*myEvent.PVertex)) || (iMuon->isolationR03().sumPt/iMuon->pt() > .1)){ continue;} //preliminary cut
   	iCut[ii].muonRecoCount++;

   	if( abs(lepton1ID)==13){
			match1DR = sqrt(dR2(iMuon->eta(), iCut[ii].lepton1Eta, iMuon->phi(), iCut[ii].lepton1Phi));
			if(mu1Match == 0 && match1DR < 0.3) {
				mu1Match++;
				iCut[ii].lepton1Matched=true;
				matchedMuon1 = &(*(iMuon));
				iCut[ii].Muon1dR=match1DR;
				if(matchedMuon1->pdgId()==iCut[ii].lepton1Id){ iCut[ii].Muon1chargeMatch = true ;}
    			else{ iCut[ii].Muon1chargeMatch = false ;}
				iCut[ii].Muon1PtRatio = iCut[ii].lepton1Pt/matchedMuon1->pt();
				//iCut[ii].Muon1TW=tWfinder(iEvent,lepton1);
				iCut[ii].Muon1Pt=matchedMuon1->pt();
				iCut[ii].Muon1Eta=matchedMuon1->eta();
				iCut[ii].Muon1Phi=matchedMuon1->phi();
			}

		}
		else if(abs(lepton2ID)==13){
			match2DR = sqrt(dR2(iMuon->eta(), iCut[ii].lepton2Eta, iMuon->phi(), iCut[ii].lepton2Phi));
			if(mu1Match == 0 && match2DR < 0.3)  {
				matchedMuon1 = &(*(iMuon));
				mu1Match++; 
				iCut[ii].lepton2Matched=true;
				iCut[ii].Muon1dR=match2DR;
				if(matchedMuon1->pdgId()==iCut[ii].lepton2Id){ iCut[ii].Muon1chargeMatch = true ;}
    			else{ iCut[ii].Muon2chargeMatch = false ;}
				iCut[ii].Muon1PtRatio = iCut[ii].lepton2Pt/matchedMuon1->pt();
				//iCut[ii].Muon1TW=tWfinder(iEvent,lepton2);
				iCut[ii].Muon1Pt=matchedMuon1->pt();
				iCut[ii].Muon1Eta=matchedMuon1->eta();
				iCut[ii].Muon1Phi=matchedMuon1->phi();
			}
		}
		
		if(mu1Match==0){
		iCut[ii].unmatchedPhi.push_back(iMuon->phi());
		iCut[ii].unmatchedEta.push_back(iMuon->eta());
		iCut[ii].unmatchedPt.push_back(iMuon->pt());
		if( abs(lepton1ID)==13){
		iCut[ii].unmatchedDR.push_back(match1DR);
		}
		else{iCut[ii].unmatchedDR.push_back(match2DR);}
		}

		}

		if(mu1Match==0 || iCut[ii].failedGenPtEta ){
			iCut[ii].failedMatch=true;
		}
	}

	//electron
	if(!iCut[ii].muonTau){
	iCut[ii].electronTrigger=false;
   if (passElectronTrig(iEvent, iCut[ii])){ iCut[ii].electronTrigger=true; }


	//for all reco electrons, loop through gen electrons to find spatial matches
	for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++){	
		if(fabs(iElectron->eta()) > 2.4) {continue;}
		if(iElectron->pt() < 10 ) {continue;}
		iCut[ii].electronRecoCount++;

		if(abs(lepton1ID)==11){
			double match1DR = sqrt(dR2(iElectron->eta(), iCut[ii].lepton1Eta, iElectron->phi(), iCut[ii].lepton1Phi));

				if(el1Match == 0 && match1DR < 0.3) {
					el1Match++;
					iCut[ii].lepton1Matched=true;
					matchedElectron1 = &(*(iElectron));
					iCut[ii].Electron1dR=match1DR;
					if(matchedElectron1->pdgId()==iCut[ii].lepton1Id){ iCut[ii].Electron1chargeMatch = true ;}
    				else{ iCut[ii].Electron1chargeMatch = false ;}
					iCut[ii].Electron1PtRatio = iCut[ii].lepton1Pt/matchedElectron1->pt();
					//iCut[ii].Electron1TW=tWfinder(iEvent,lepton1);
					iCut[ii].Electron1Pt=matchedElectron1->pt();
					iCut[ii].Electron1Eta=matchedElectron1->eta();
					iCut[ii].Electron1Phi=matchedElectron1->phi();
				}
		}
		else if(abs(lepton2ID)==11){
			double match2DR = sqrt(dR2(iElectron->eta(), iCut[ii].lepton2Eta, iElectron->phi(), iCut[ii].lepton2Phi));
				if(el1Match == 0 && match2DR < 0.3)  {
					matchedElectron1 = &(*(iElectron));
					el1Match++; 
					iCut[ii].lepton2Matched=true;
					iCut[ii].Electron1dR=match2DR;
					if(matchedElectron1->pdgId()==iCut[ii].lepton2Id){ iCut[ii].Electron1chargeMatch = true ;}
    				else{ iCut[ii].Electron1chargeMatch = false ;}
					iCut[ii].Electron1PtRatio = iCut[ii].lepton2Pt/matchedElectron1->pt();
					//iCut[ii].Electron1TW=tWfinder(iEvent,lepton2);
					iCut[ii].Electron1Pt=matchedElectron1->pt();
					iCut[ii].Electron1Eta=matchedElectron1->eta();
					iCut[ii].Electron1Phi=matchedElectron1->phi();
				}
		}
		if(el1Match==0){
			iCut[ii].unmatchedPhi.push_back(iElectron->phi());
			iCut[ii].unmatchedEta.push_back(iElectron->eta());
			iCut[ii].unmatchedPt.push_back(iElectron->pt());
			if( abs(lepton1ID)==11){
			iCut[ii].unmatchedDR.push_back(match1DR);
			}
			else{iCut[ii].unmatchedDR.push_back(match2DR);}
			}

	}
	if(el1Match==0 || iCut[ii].failedGenPtEta ){
	iCut[ii].failedMatch=true;
	}

	}	
	}

m_histoMaker.fill(iCut[ii],ii);
}
//m_histoMaker.fill(cuts1,0);
//m_histoMaker.fill(cuts2,1);
//m_histoMaker.fill(cuts3,2);
//m_histoMaker.fill(cuts4,3);

//std::cout << iCut[ii]Muon1Pt<<std::endl;

}



//HELPERS
double gen_match::dR2(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    double dphi = dPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
}
double gen_match::dPhi(double phi1, double phi2) {
    double raw_dphi = phi1 - phi2;
    if (fabs(raw_dphi) < ROOT::Math::Pi()) return raw_dphi;
    double region = std::round(raw_dphi / (2.*ROOT::Math::Pi()));
    return raw_dphi - 2.*ROOT::Math::Pi()*region;
}

//check if a lepton can be traced back to a W boson and then to a top
bool gen_match::tWfinder(const edm::Event& iEvent, const reco::GenParticle* lepton) {

		//edm::Handle<std::vector<reco::GenParticle>> genParticles;
		//iEvent.getByToken(m_genParticleToken, genParticles);

    		bool ttbar=false;
    		int iStatus = 0;

    		if((abs(lepton->pdgId())==11 || abs(lepton->pdgId())==13 ) && lepton->status()!=1){
    			int n = lepton->numberOfDaughters();
    			//std::cout<<lepton->pdgId();
    			//std::cout<<" -> ";
    			//std::cout<<"number of daughters: ";
    			//std::cout<<n<<std::endl;
    			//std::cout<<"pdgIDs -> "<<std::endl;
    			   for(int j = 0; j < n; ++ j) {
       				int dauId = abs(lepton->daughter(j)->pdgId());
       				//std::cout<<dauId<<std::endl;
			
    					if(dauId==11 || dauId==13){return false; }//std::cout<<"FAILED"<<std::endl; std::cout<<"------------------------------"<<std::endl; } //the lepton is not the final state lepton, it decays into another muon or electron
					}
			}

    		const reco::Candidate* iParticle = lepton->mother();

    		//while(iParticle->pdgId()!=2212){

    		while(iStatus!=4){  //status=4 is the initial proton
    			iStatus = iParticle->status();
    			
    			if(abs(iParticle->pdgId())==24){ //found W
    				//while(iParticle->pdgId()!=2212){
    				while(iStatus!=4){
    				iParticle = iParticle->mother();
    				iStatus = iParticle->status();

    			   	if(abs(iParticle->pdgId())==6){ ttbar=true; 
    			   		break;
    			   	}
					}

    			}

    			iParticle = iParticle->mother();
    		}

		if(ttbar==true){//std::cout <<"found T->W lepton!"<<std::endl; 
		return true;}
		else{//std::cout <<"no T->W lepton :("<<std::endl; 
		return false;}
}

bool gen_match::passElectronTrig(const edm::Event& iEvent, eventBits2& iCut) {
  bool passTriggers = false;

  //std::cout <<"checking electron trigger paths "<<std::endl;
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(m_trigResultsToken, triggerResults);

  //std::cout <<"grabbing trigger names"<<std::endl;
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults); //SEGFAULT
  //std::cout <<"grabbing trigger results "<<std::endl;
  //std::cout <<"looping over paths to pass"<<std::endl;
  for(size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    for(auto& pathName : m_electronPathsToPass){
      if((name.find(pathName) != std::string::npos )){
        if(triggerResults->accept(i)){
          passTriggers = true;
        }
      }
    }
  }
  //if(passTriggers) iCut[ii].electronTrigger = 1.;
  //else iCut[ii].electronTrigger = 0.;

  return passTriggers;
}


// ------------ method called once each job just before starting event loop  ------------
void
gen_match::beginJob() {
	edm::Service<TFileService> fs; 

	fs->mkdir("cuts1");
	fs->mkdir("cuts2");
	fs->mkdir("cuts3");
	fs->mkdir("cuts4");

	m_histoMaker.book(fs->mkdir("cuts1/twoMuons"),0,0);
	//m_iCut[ii].book(fs->mkdir("twoElectrons"),1);
	m_histoMaker.book(fs->mkdir("cuts1/muonElectron"),0,1);
	m_histoMaker.book(fs->mkdir("cuts1/muonTau"),0,2);
	m_histoMaker.book(fs->mkdir("cuts1/allEvents"),0,3);
	//m_iCut[ii].book(fs->mkdir("electronTau"),4);

	m_histoMaker.book(fs->mkdir("cuts2/twoMuons"),1,0);
	//m_cuts2.book(fs->mkdir("twoElectrons"),6);
	m_histoMaker.book(fs->mkdir("cuts2/muonElectron"),1,1);
	m_histoMaker.book(fs->mkdir("cuts2/muonTau"),1,2);
	m_histoMaker.book(fs->mkdir("cuts2/allEvents"),1,3);
	//m_cuts2.book(fs->mkdir("electronTau")9);

	m_histoMaker.book(fs->mkdir("cuts3/twoMuons"),2,0);
	//m_cuts3.book(fs->mkdir("twoElectrons"),11);
	m_histoMaker.book(fs->mkdir("cuts3/muonElectron"),2,1);
	m_histoMaker.book(fs->mkdir("cuts3/muonTau"),2,2);
	m_histoMaker.book(fs->mkdir("cuts3/allEvents"),2,3);
	//m_cuts3.book(fs->mkdir("electronTau"),14);

	m_histoMaker.book(fs->mkdir("cuts4/twoMuons"),3,0);
	//m_cuts4.book(fs->mkdir("twoElectrons"),16);
	m_histoMaker.book(fs->mkdir("cuts4/muonElectron"),3,1);
	m_histoMaker.book(fs->mkdir("cuts4/muonTau"),3,2);
	m_histoMaker.book(fs->mkdir("cuts4/allEvents"),3,3);
	//m_cuts4.book(fs->mkdir("electronTau"),19);

}


// ------------ method called once each job just after ending the event loop  ------------
void
gen_match::endJob() {
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
gen_match::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(gen_match);
