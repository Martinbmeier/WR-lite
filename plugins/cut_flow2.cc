
// -*- C++ -*-
//
// Package:    ExoAnalysis/cut_flow2
// Class:      cut_flow2
//
/**\class cut_flow2 cut_flow2.cc ExoAnalysis/cut_flow2/plugins/cut_flow2.cc
 Description: cut flow, event selection, and extraction of variables for NN
 Implementation:
     [Notes on implementation]
*/
//
//  Author:  Andrew Evans, adapted by Martin Meier
//         Created:  Aug 2021
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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "HEEP/VID/interface/CutNrs.h"
#include "HEEP/VID/interface/VIDCutCodes.h"

#include "TLorentzVector.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <cmath>

#include "TH1.h"
#include "TMath.h"
#include "TDirectory.h"


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/eventInfo.h"
//#include "ExoAnalysis/WR_lite/interface/eventHistos2.h"
#include "ExoAnalysis/WR_lite/interface/cutFlowHistos.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class cut_flow2 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit cut_flow2(const edm::ParameterSet&);
		~cut_flow2();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		double dR2(double eta1, double eta2, double phi1, double phi2);
		double dPhi(double phi1, double phi2);
		bool tWfinder(const edm::Event&, const reco::GenParticle* );
		bool passElectronTrig(const edm::Event&);
		void csvTable(double genMuonPt, double genElectronPt, const pat::Muon*, const pat::Electron*, const pat::Jet*, const pat::Jet* , const pat::MET);
		//double transverseSphericity(math::XYZTLorentzVector p1, math::XYZTLorentzVector p2, math::XYZTLorentzVector p3);
		//void saveElectronData(eventBits2 * iBit, double matched1Mass, double matched2Mass);
		//void saveMuonData(eventBits2 * iBit, double matched1Mass, double matched2Mass);
		
		
		cutFlowHistos m_histoMaker;
		TH1D* m_eventsWeight;
	  TH2D* m_cosJets;
	  TH2D* m_cosJet1electron;
	  TH2D* m_cosJet2electron; 
	  TH2D* m_cosJet1muon; 
	  TH2D* m_cosJet2muon; 
	  TH2D* m_cosLeptons; 
	  TH2D* m_dRjets;
	  TH2D* m_dRLeptons;
	  TH2D* m_dRJet1muon;
	  TH2D* m_dRJet2muon;
	  TH2D* m_dRJet1electron;
	  TH2D* m_dRJet2electron;

		//neuralNet networkResolved = neuralNet("/home/kronh006/Version3/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/data/Resolved");
		//neuralNet networkSuperResolved = neuralNet("/home/kronh006/Version3/CMSSW_10_4_0_patch1/src/ExoAnalysis/WR_lite/data/SuperResolved");

		// ----------member data ---------------------------

		edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
		edm::EDGetToken m_genParticleToken;
		edm::EDGetToken m_recoMETToken;
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

		//std::ofstream myfile;
		
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
cut_flow2::cut_flow2(const edm::ParameterSet& iConfig)
	:
	tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	m_genParticleToken(consumes<std::vector<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("genParticles"))),
	m_recoMETToken(consumes<std::vector<pat::MET>> (iConfig.getParameter<edm::InputTag>("recoMET"))),
	m_highMuonToken (consumes<std::vector<pat::Muon>> (iConfig.getParameter<edm::InputTag>("highMuons"))),
	m_highElectronToken (consumes<std::vector<pat::Electron>> (iConfig.getParameter<edm::InputTag>("highElectrons"))),
	m_AK4recoCHSJetsToken (consumes<std::vector<pat::Jet>> (iConfig.getParameter<edm::InputTag>("AK4recoCHSJets"))),
	m_genEventInfoToken (consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genInfo"))),
	m_offlineVerticesToken (consumes<std::vector<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("vertices"))),
	m_dataSaveFile (iConfig.getUntrackedParameter<std::string>("trainFile")),
	m_isSignal (iConfig.getUntrackedParameter<bool>("isSignal"))

{
   //now do what ever initialization is needed

   m_electronPathsToPass  = iConfig.getParameter<std::vector<std::string> >("electronPathsToPass");
   m_trigResultsToken = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("trigResults"));
}


cut_flow2::~cut_flow2()
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
cut_flow2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


	bool oneElectronMuon=false; 

	bool electronTrigger=false;  //cut1
	bool oneHeepElectron=false;  //cut2
	bool oneMuonHighpT=false;    //cut3
	bool twoJets=false;				   //cut4
	bool angularSeparation=false;//cut5
	bool electronHighPt=false;	 //cut6
	bool oneBTag=false;				   //cut7
	bool twoBTag=false;				   //cut8
	bool muonIsolation1=false;	 //cut9
	bool muonIsolation2=false;	 //cut10

	eventBits2 iBit; 
   
  eventInfo myEvent; 
	
	edm::Handle<GenEventInfoProduct> eventInfo;
	iEvent.getByToken(m_genEventInfoToken, eventInfo);

	edm::Handle<std::vector<pat::Muon>> highMuons;
	iEvent.getByToken(m_highMuonToken, highMuons);

	edm::Handle<std::vector<pat::Electron>> highElectrons;
	iEvent.getByToken(m_highElectronToken, highElectrons);

	edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByToken(m_genParticleToken, genParticles);

	edm::Handle<std::vector<pat::MET>> recoMET;
	iEvent.getByToken(m_recoMETToken, recoMET);
	
  
	float eventCount = eventInfo->weight()/fabs(eventInfo->weight());
	//double eventWeight = eventInfo->weight();
	double eventWeight = 1;
	
	
	edm::Handle<std::vector<reco::Vertex>> vertices;
	iEvent.getByToken(m_offlineVerticesToken, vertices);
	if(!vertices.isValid()) {
		throw cms::Exception("Vertex collection not valid!");
	}

	iBit.hasPVertex = myEvent.PVselection(vertices);






//start reco

	const pat::MET Met = recoMET->front();

	edm::Handle<std::vector<pat::Jet>> recoJetsAK4;  
	iEvent.getByToken(m_AK4recoCHSJetsToken, recoJetsAK4);  


	int jetCount = 0;
	int btagcount = 0;

	const pat::Jet* Jet1=0;
	const pat::Jet* Jet2=0;




	//Get jets with maximum pt
		for(std::vector<pat::Jet>::const_iterator iJet = recoJetsAK4->begin(); iJet != recoJetsAK4->end(); iJet++) {

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
			if(BJP > 0.4184){ btagcount++; }		
			if(jetCount==0){Jet1=&(*(iJet));}
			if(jetCount==1){Jet2=&(*(iJet));}
			jetCount++;
		}

		if(jetCount>1){twoJets=true;}
		if(btagcount>0){oneBTag=true;}
		if(btagcount>1){twoBTag=true;}




//gen lepton info

	double genMuonpT=-1000;
	//double newGenMuonPt;

	double genElectronpT=-1000;
	//double newGenElectronPt;

	bool genMuon = false;
	bool genElectron = false;


for (std::vector<reco::GenParticle>::const_iterator iParticle = genParticles->begin(); iParticle != genParticles->end(); iParticle++) {
	if( ! iParticle->isHardProcess() ){ continue; }
	if(abs(iParticle->pdgId())==13 && !genMuon){
		//if(!tWfinder(iEvent, &(*iParticle))){ continue; } //could check if the gen particle comes from a top->W->lepton
		genMuon = true;
		genMuonpT=iParticle->pt();
		//if(newGenMuonPt>genMuonpT){genMuonpT=newGenMuonPt;} //get the highest pt gen muon if multiple
	}
	if(abs(iParticle->pdgId())==11 && !genElectron){
		genElectron=true;
		genElectronpT = iParticle->pt();
		//if(newGenElectronPt>genElectronpT){genElectronpT=newGenElectronPt;} //get the highest pt gen electron if multiple
	}
}

if(genElectron && genMuon){oneElectronMuon=true;}


//check electron trigger
if (passElectronTrig(iEvent)){ electronTrigger=true; }

//muon/electron reconstruction

	//muon reco

	  const pat::Muon* recoMuon=0;

		//double recoMuonpT = -1000;
		//double newrecoMuonpT=-1000;
		bool foundMuon=false;


   	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){

   		//if(fabs(iMuon->eta()) > 2.4 || iMuon->tunePMuonBestTrack()->pt() < 10 || !(iMuon->isHighPtMuon(*myEvent.PVertex)) || (iMuon->isolationR03().sumPt/iMuon->pt() > .1)){ continue;} //preliminary cut
   		
   		if(iMuon->isHighPtMuon(*myEvent.PVertex) && !oneMuonHighpT){oneMuonHighpT=true; recoMuon=&(*(iMuon));}

   	

   		//if(newrecoMuonpT>recoMuonpT){recoMuonpT=newrecoMuonpT; recoMuon=&(*(iMuon)); foundMuon=true;} //get the highest pt muon if multiple

		}

		if(oneMuonHighpT){
			if(recoMuon->passed(reco::Muon::TkIsoLoose)){muonIsolation1=true;}
			if(recoMuon->passed(reco::Muon::TkIsoTight)){muonIsolation1=true; muonIsolation2=true; }
		}

	//electron reco

		const pat::Electron* recoElectron=0;
		bool foundRecoElectron=false;
		//double recoElectronpT = -1000;
		//double newrecoElectronpT = -1000;

			for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++){	
				//if(fabs(iElectron->eta()) > 2.4) {continue;}
				//if(iElectron->pt() < 10 ) {continue;}
				
				//recoElectronpT=iElectron->pt();
				if(!foundRecoElectron){
					oneHeepElectron=true;
					recoElectron=&(*(iElectron));
					foundRecoElectron=true;

					if(recoElectron->pt()>75){electronHighPt=true;}
				}
   			//if(newrecoElectronpT>recoElectronpT){recoElectronpT=newrecoElectronpT; recoElectron=&(*(iElectron)); oneHeepElectron=true;} //get the highest pt electron if multiple

			}


			if(oneHeepElectron && foundMuon && twoJets){
				double dileptonSeparation=sqrt(dR2(recoMuon->eta(), recoElectron->eta(), recoMuon->phi(), recoElectron->phi()));
			   double muonJet1Sep=sqrt(dR2(Jet1->eta(), recoMuon->eta(), Jet1->phi(), recoMuon->phi()));
			   double muonJet2Sep=sqrt(dR2(Jet2->eta(), recoMuon->eta(), Jet2->phi(), recoMuon->phi()));
			   double electronJet1Sep=sqrt(dR2(Jet1->eta(), recoElectron->eta(), Jet1->phi(), recoElectron->phi()));
				double electronJet2Sep=sqrt(dR2(Jet2->eta(), recoElectron->eta(), Jet2->phi(), recoElectron->phi()));
				double jetSeparation=sqrt(dR2(Jet2->eta(), Jet1->eta(), Jet2->phi(), Jet1->phi()));
				if(dileptonSeparation>0.4 && muonJet1Sep>0.4 && muonJet2Sep>0.4 && electronJet1Sep > 0.4 && electronJet2Sep>0.4 && jetSeparation>0.4){angularSeparation=true;} // jet/lepton separation cut
			} 

		
	m_eventsWeight->Fill(0.5, eventCount);

	double em_ratio=genMuonpT/genElectronpT;

if(oneElectronMuon){// || !oneElectronMuon){
	m_histoMaker.fill(genMuonpT,0,eventWeight);
	if(electronTrigger){
		m_histoMaker.fill(genMuonpT,1,eventWeight);
		if(oneHeepElectron){
			m_histoMaker.fill(genMuonpT,2,eventWeight);
			if(oneMuonHighpT){
				m_histoMaker.fill(genMuonpT,3,eventWeight);
				if(twoJets){
					m_histoMaker.fill(genMuonpT,4,eventWeight);
					if(angularSeparation){
						m_histoMaker.fill(genMuonpT,5,eventWeight);
						if(electronHighPt){
							m_histoMaker.fill(genMuonpT,6,eventWeight);

								//csvTable(genMuonpT,genElectronpT,recoMuon,recoElectron,Jet1,Jet2,Met);  //fill a csv table with variables for the NN 

								m_cosJets->Fill(TMath::Cos(deltaPhi(Jet2->phi(),Jet1->phi())),em_ratio,1);
								m_cosLeptons->Fill(TMath::Cos(deltaPhi(recoMuon->phi(),recoElectron->phi())),em_ratio,1);


							if(Jet1->pt()>Jet2->pt()){
								m_cosJet1electron->Fill(TMath::Cos(deltaPhi(Jet1->phi(),recoElectron->phi())),em_ratio,1);
								m_cosJet2electron->Fill(TMath::Cos(deltaPhi(Jet2->phi(),recoElectron->phi())),em_ratio,1);
								m_cosJet1muon->Fill(TMath::Cos(deltaPhi(Jet1->phi(),recoMuon->phi())),em_ratio,1);
								m_cosJet2muon->Fill(TMath::Cos(deltaPhi(Jet2->phi(),recoMuon->phi())),em_ratio,1);

							}
							else{
								m_cosJet1electron->Fill(TMath::Cos(deltaPhi(Jet2->phi(),recoElectron->phi())),em_ratio,1);
								m_cosJet2electron->Fill(TMath::Cos(deltaPhi(Jet1->phi(),recoElectron->phi())),em_ratio,1);
								m_cosJet1muon->Fill(TMath::Cos(deltaPhi(Jet2->phi(),recoMuon->phi())),em_ratio,1);
								m_cosJet2muon->Fill(TMath::Cos(deltaPhi(Jet1->phi(),recoMuon->phi())),em_ratio,1);
							}

							m_dRjets->Fill(deltaR(Jet1->eta(),Jet1->phi(),Jet2->eta(),Jet2->phi()),em_ratio,1);
							m_dRLeptons->Fill(deltaR(recoMuon->eta(),recoMuon->phi(),recoElectron->eta(),recoElectron->phi()),em_ratio,1);

							if(Jet1->pt()>Jet2->pt()){
								m_dRJet1muon->Fill(deltaR(Jet1->eta(),Jet1->phi(),recoMuon->eta(),recoMuon->phi()),em_ratio,1);
								m_dRJet2muon->Fill(deltaR(Jet2->eta(),Jet2->phi(),recoMuon->eta(),recoMuon->phi()),em_ratio,1);

								m_dRJet1electron->Fill(deltaR(Jet1->eta(),Jet1->phi(),recoElectron->eta(),recoElectron->phi()),em_ratio,1);
								m_dRJet2electron->Fill(deltaR(Jet2->eta(),Jet2->phi(),recoElectron->eta(),recoElectron->phi()),em_ratio,1);
							}
							else{
								m_dRJet1muon->Fill(deltaR(Jet2->eta(),Jet2->phi(),recoMuon->eta(),recoMuon->phi()),em_ratio,1);
								m_dRJet2muon->Fill(deltaR(Jet1->eta(),Jet1->phi(),recoMuon->eta(),recoMuon->phi()),em_ratio,1);

								m_dRJet1electron->Fill(deltaR(Jet2->eta(),Jet2->phi(),recoElectron->eta(),recoElectron->phi()),em_ratio,1);
								m_dRJet2electron->Fill(deltaR(Jet1->eta(),Jet1->phi(),recoElectron->eta(),recoElectron->phi()),em_ratio,1);
							}

						
							if(oneBTag){
								m_histoMaker.fill(genMuonpT,7,eventWeight);
								if(twoBTag){
									m_histoMaker.fill(genMuonpT,8,eventWeight);
									if(muonIsolation1){
										m_histoMaker.fill(genMuonpT,9,eventWeight);
										if(muonIsolation2)
											m_histoMaker.fill(genMuonpT,10,eventWeight);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


}





//HELPERS
double cut_flow2::dR2(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    double dphi = dPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
}
double cut_flow2::dPhi(double phi1, double phi2) {
    double raw_dphi = phi1 - phi2;
    if (fabs(raw_dphi) < ROOT::Math::Pi()) return raw_dphi;
    double region = std::round(raw_dphi / (2.*ROOT::Math::Pi()));
    return raw_dphi - 2.*ROOT::Math::Pi()*region;
}

//check if a lepton can be traced back to a W boson and then to a top
bool cut_flow2::tWfinder(const edm::Event& iEvent, const reco::GenParticle* lepton) {

    		bool ttbar=false;
    		int iStatus;

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

bool cut_flow2::passElectronTrig(const edm::Event& iEvent) {
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

  return passTriggers;
}

void cut_flow2::csvTable(double genMuonPt, double genElectronPt, const pat::Muon* muon, const pat::Electron* electron, const pat::Jet* jet1, const pat::Jet* jet2, const pat::MET Met) {

std::ofstream myfile;
myfile.open("neuralNetData.csv",std::ios_base::app);
myfile << muon->phi() << ", "
       << muon->eta() << ", "
       << electron->pt() << ", "
       << electron->phi() << ", "
       << electron->eta() << ", "
       << jet1->pt() << ", "
       << jet1->phi() << ", "
       << jet1->eta() << ", "
       << jet2->pt() << ", "
       << jet2->phi() << ", "
       << jet2->eta() << ", "
       << Met.pt() << ", "
       << Met.phi() <<", "
       << genMuonPt/genElectronPt << "\n ";

myfile.close();

}


// ------------ method called once each job just before starting event loop  ------------
void
cut_flow2::beginJob() {

	std::ofstream myfile;

	myfile.open("neuralNetData.csv",std::ios_base::app);
	myfile<<"muon phi, muon eta, electron pt, electron phi, electron eta, jet 1 pt, jet 1 phi, jet 1 eta, jet 2 pt, jet 2 phi, jet 2 eta, MET pt, MET phi, gen muon/electron pt ratio\n";
	myfile.close();

	edm::Service<TFileService> fs; 

	TFileDirectory countFolder = fs->mkdir("event_count");

	TFileDirectory variableCorrelations = fs->mkdir("variable_correlations");
	
	m_histoMaker.book(fs->mkdir("cuts1"),0);

	m_histoMaker.book(fs->mkdir("cuts2"),1);

	m_histoMaker.book(fs->mkdir("cuts3"),2);
	
	m_histoMaker.book(fs->mkdir("cuts4"),3);

	m_histoMaker.book(fs->mkdir("cuts5"),4);

	m_histoMaker.book(fs->mkdir("cuts6"),5);

	m_histoMaker.book(fs->mkdir("cuts7"),6);

	m_histoMaker.book(fs->mkdir("cuts8"),7);

	m_histoMaker.book(fs->mkdir("cuts9"),8);

	m_histoMaker.book(fs->mkdir("cuts10"),9);

	m_histoMaker.book(fs->mkdir("cuts11"),10);


	m_eventsWeight = {countFolder.make<TH1D>("eventsWeight","number of events weighted", 1, 0.0, 1)};

	m_cosJets = {variableCorrelations.make<TH2D>("cosJets","muon electron pt ratio vs. cos(delta Phi) for two jets",50,-1,1,50,0,2.5)};
	m_cosJet1electron = {variableCorrelations.make<TH2D>("cosJet1Electron","muon electron pt ratio vs. cos(delta Phi) for jet1 and electron",50,-1,1,50,0,2.5)};
	m_cosJet2electron = {variableCorrelations.make<TH2D>("cosJet2Electron","muon electron pt ratio vs. cos(delta Phi) for jet2 and electron",50,-1,1,50,0,2.5)};
	m_cosJet1muon = {variableCorrelations.make<TH2D>("cosJet1Muon","muon electron pt ratio vs. cos(delta Phi) for jet1 and muon",50,-1,1,50,0,2.5)};
	m_cosJet2muon = {variableCorrelations.make<TH2D>("cosJet2Muon","muon electron pt ratio vs. cos(delta Phi) for jet2 and muon",50,-1,1,50,0,2.5)};
	m_cosLeptons = {variableCorrelations.make<TH2D>("cosleptons","muon electron pt ratio vs. cos(delta Phi) for electron and muon",50,-1,1,50,0,2.5)};

	m_dRjets = {variableCorrelations.make<TH2D>("dRjets","muon electron pt ratio vs. deltaR for jets",50,0,7.5,50,0,2.5)};
	m_dRLeptons = {variableCorrelations.make<TH2D>("dRLeptons","muon electron pt ratio vs. deltaR for muon and electron",50,0,7.5,50,0,2.5)};
	m_dRJet1muon = {variableCorrelations.make<TH2D>("dRJet1muon","muon electron pt ratio vs. deltaR for jet1 and muon",50,0,7.5,50,0,2.5)};
	m_dRJet2muon = {variableCorrelations.make<TH2D>("dRJet2muon","muon electron pt ratio vs. deltaR for jet2 and muon",50,0,7.5,50,0,2.5)};
	m_dRJet1electron = {variableCorrelations.make<TH2D>("dRJet1electron","muon electron pt ratio vs. deltaR for jet1 and electron",50,0,7.5,50,0,2.5)};
	m_dRJet2electron = {variableCorrelations.make<TH2D>("dRJet2electron","muon electron pt ratio vs. deltaR for jet2 and electron",50,0,7.5,50,0,2.5)};






	
}


// ------------ method called once each job just after ending the event loop  ------------
void
cut_flow2::endJob() {
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
cut_flow2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(cut_flow2);
