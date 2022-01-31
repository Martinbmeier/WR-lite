
// -*- C++ -*-
//
// Package:    ExoAnalysis/ttbarstudy
// Class:      cut_flow2
//
/**\class ttbarstudy ttbarstudy.cc ExoAnalysis/ttbarstudy/plugins/ttbarstudy.cc
 Description: ttbarstudy
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

#include "ExoAnalysis/WR_lite/interface/eventInfo.h"
#include "ExoAnalysis/WR_lite/interface/ttbarHistos.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ttbarstudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit ttbarstudy(const edm::ParameterSet&);
		~ttbarstudy();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		double dR2(double eta1, double eta2, double phi1, double phi2);
		double dPhi(double phi1, double phi2);
		bool tWfinder(const edm::Event&, const reco::GenParticle* );
		bool passElectronTrig(const edm::Event&);
		void csvTable(double genMuonPt, double genElectronPt, const pat::Muon*, const pat::Electron*, const pat::Jet*, const pat::Jet*, const pat::Jet*, const pat::Jet*, math::XYZTLorentzVector combinedJets, const pat::MET, double weight);
		void csvTable(const pat::Muon*, const pat::Electron*, const pat::Jet*, const pat::Jet*, const pat::Jet*, const pat::Jet*, math::XYZTLorentzVector combinedJets, const pat::MET, double weight);
		double misIDrate(const pat::Muon*);
		
		
		ttbarHistos m_histoMaker;

		TH1D* m_eventsWeight;

	  // TH2D* m_cosJets;
	  // TH2D* m_cosJet1electron;
	  // TH2D* m_cosJet2electron; 
	  // TH2D* m_cosJet1muon; 
	  // TH2D* m_cosJet2muon; 
	  // TH2D* m_cosLeptons; 
	  // TH2D* m_dRjets;
	  // TH2D* m_dRLeptons;
	  // TH2D* m_dRJet1muon;
	  // TH2D* m_dRJet2muon;
	  // TH2D* m_dRJet1electron;
	  // TH2D* m_dRJet2electron;

	  // TH1D* m_deltaPhiJet2Electron;
	  // TH1D* m_deltaPhiLeptons;
	  // TH1D* m_deltaPhiJet1Muon;
	  // TH1D* m_deltaPhiJets;
	  // TH1D* m_deltaPhiJet1Electron;
	  // TH1D* m_deltaPhiJet2Muon;

	  // TH2D* m_cosMetJet1;
	  // TH2D* m_cosMetJet2;
	  // TH2D* m_cosMetElectron;
	  // TH2D* m_cosMetMuon;

	  // TH1D* m_deltaPhiMetJet1;
	  // TH1D* m_deltaPhiMetJet2;
	  // TH1D* m_deltaPhiMetMuon;
	  // TH1D* m_deltaPhiMetElectron; 

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

		double binEdges[17] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 250.0, 300.0, 350.0, 400.0, 1000.0};
		double rates[16] = {0.001, 0.002, 0.003, 0.005, 0.01, 0.012, 0.015, 0.018, 0.02, 0.021, 0.024, 0.026, 0.029, 0.04, 0.06, 0.07};

		
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
ttbarstudy::ttbarstudy(const edm::ParameterSet& iConfig)
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


ttbarstudy::~ttbarstudy()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ttbarstudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	//bool background = true;
	//background = !m_isSignal;

	bool electronTrigger=false;  
	// bool twoJets=false;			  
	// bool angularSeparation=false;
	// bool angularSeparation2B=false;
	// bool angularSeparation1B=false;
	// bool oneBTag=false;			  
	// bool twoBTag=false;		
	// bool dileptonMass=false;	  

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
	double eventWeight = eventInfo->weight();
	//double eventCountmiss = 0;
	
	edm::Handle<std::vector<reco::Vertex>> vertices;
	iEvent.getByToken(m_offlineVerticesToken, vertices);
	if(!vertices.isValid()) {
		throw cms::Exception("Vertex collection not valid!");
	}

	iBit.hasPVertex = myEvent.PVselection(vertices);


//start reco

	const pat::MET Met = recoMET->front();

	math::XYZTLorentzVector combinedJetsP4 = {0., 0., 0., 0.};

	edm::Handle<std::vector<pat::Jet>> recoJetsAK4;  
	iEvent.getByToken(m_AK4recoCHSJetsToken, recoJetsAK4);  


	// int jetCount = 0;
	// // int btagcount = 0;

	// // const pat::Jet* bJet1=0;
	// // const pat::Jet* bJet2=0;

	// const pat::Jet* Jet1=0;
	// const pat::Jet* Jet2=0;


	// //Get jets with maximum pt
	// 	for(std::vector<pat::Jet>::const_iterator iJet = recoJetsAK4->begin(); iJet != recoJetsAK4->end(); iJet++) {

	// 		combinedJetsP4=combinedJetsP4+iJet->p4();

	// 		double NHF  =           iJet->neutralHadronEnergyFraction();
	// 		double NEMF =           iJet->neutralEmEnergyFraction();
	// 		double CHF  =           iJet->chargedHadronEnergyFraction();
	// 		double CEMF =           iJet->chargedEmEnergyFraction();
	// 		double NumConst =       iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
	// 		double MUF      =       iJet->muonEnergyFraction();
	// 		double EUF      =       iJet->electronEnergyFraction();
	// 		double CHM      =       iJet->chargedMultiplicity();
	// 		// double BJP		 =       iJet->bDiscriminator(cSV_bTag1) + iJet->bDiscriminator(cSV_bTag2);
	// 		//APPLYING TIGHT QUALITY CUTS
	// 		if (NHF > .9) continue;
	// 		if (NEMF > .9) continue;
	// 		if (NumConst <= 1) continue;
	// 		if (MUF >= .8) continue; //MAKE SURE THE AREN'T MUONS
	// 		if (EUF >= .8) continue; //MAKE SURE THE AREN'T ELECTRONS
	// 		//ADDITIONAL CUTS BECAUSE OF TIGHT ETA CUT
	// 		if (CHF == 0) continue;
	// 		if (CHM == 0) continue;
	// 		//if (CEMF > .99) continue;
	// 		if (CEMF > .90)  continue;

	// 		if(jetCount==0){
	// 			Jet1=&(*(iJet));
	// 		}
	// 		if(jetCount==1){
	// 			Jet2=&(*(iJet));
	// 		}
	// 		jetCount++;

			// if(BJP > 0.4184){ 
	

			// 	if(btagcount==0){
			// 		bJet1=&(*(iJet));
			// 	}
			// 	if(btagcount==1){
			// 		bJet2=&(*(iJet));
			// 	}
			// 	btagcount++; 
			// }
		//}

		// if(jetCount>1){twoJets=true;}
		// if(btagcount>0){oneBTag=true; oneBTag=false;}
		// if(btagcount>1){twoBTag=true; twoBTag=false;}



//gen lepton info

	double genMuonpT=-1000;
	double genElectronpT=-1000;

	int genLeptonCount = 0;

	for (std::vector<reco::GenParticle>::const_iterator iParticle = genParticles->begin(); iParticle != genParticles->end(); iParticle++) {
		if( ! iParticle->isHardProcess() ){ continue; }
		if( ! tWfinder(iEvent, &(*iParticle))){ continue; }  //could check if the gen particle comes from a top->W->lepton
		if(abs(iParticle->pdgId())==13 || abs(iParticle->pdgId())==11){
			genLeptonCount+=1;
			if(abs(iParticle->pdgId())==13 && genMuonpT<0){genMuonpT=iParticle->pt();}
			else if(abs(iParticle->pdgId())==11 && genElectronpT<0){genElectronpT=iParticle->pt();}
		}	
	}


	//	//check electron trigger
if (passElectronTrig(iEvent)){ electronTrigger=true; }

//muon/electron reconstruction

	// int leptonCount = 0;


	// //muon reco

	//   const pat::Muon* recoMuon1=0;
	//   const pat::Muon* recoMuon2=0;
	//   int muonCount = 0;
	//   double recoMuonpT = -1000;

 //   	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){

 //   		if(!(iMuon->isHighPtMuon(*myEvent.PVertex))) continue; // || !iMuon->passed(reco::Muon::TkIsoTight)) continue; //preliminary cut
   		
 //   		if(muonCount==0){ 
 //   				recoMuon1=&(*(iMuon)); muonCount += 1;
 //   		}
 //   		else if(muonCount==1){ 
 //   				recoMuon2=&(*(iMuon)); muonCount += 1; 
 //   		}

 //   		leptonCount += 1;
   		
	// 	}

	//electron reco

		// const pat::Electron* recoElectron1=0;
		// const pat::Electron* recoElectron2=0;
		// int electronCount = 0;

		// //double recoElectronpT = -1000;
		// //double newrecoElectronpT = -1000;

		// 	for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++){	

		// 		const vid::CutFlowResult* vidResult =  iElectron->userData<vid::CutFlowResult>("heepElectronID_HEEPV70");
		// 		const bool heepIDVID = vidResult->cutFlowPassed();
		// 		if (heepIDVID == false) continue;
				
		// 		//recoElectronpT=iElectron->pt();
		// 		if(electronCount==0){ 
		// 				recoElectron1=&(*(iElectron)); electronCount+=1; 
		// 		}
		// 		else if(electronCount==1){ 
		// 			recoElectron2=&(*(iElectron)); electronCount+=1;
		// 		}

		// 		leptonCount+=1;
				
		// 	}

			m_eventsWeight->Fill(0.5, eventCount);

			if(genLeptonCount == 2 && electronTrigger){
			
				m_histoMaker.fill(genMuonPT, genElectronPt, eventCount);

		}


	// if(oneBTag == true && angularSeparation1B == true){
	// 	m_histoMaker.fill(recoMuonpT,invMassAS_mu_mu,invMassSS_mu_mu,invMassOS_mu_mu,invMassSS_e_e,invMassOS_e_e,invMassSS_mu_e,invMassOS_mu_e,1,eventCount);
	// }
	// if(twoBTag == true && angularSeparation2B == true){
	// 	m_histoMaker.fill(recoMuonpT,invMassAS_mu_mu,invMassSS_mu_mu,invMassOS_mu_mu,invMassSS_e_e,invMassOS_e_e,invMassSS_mu_e,invMassOS_mu_e,2,eventCount);
	// }

	
	
}





//HELPERS
double ttbarstudy::dR2(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    double dphi = dPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
}
double ttbarstudy::dPhi(double phi1, double phi2) {
    double raw_dphi = phi1 - phi2;
    if (fabs(raw_dphi) < ROOT::Math::Pi()) return raw_dphi;
    double region = std::round(raw_dphi / (2.*ROOT::Math::Pi()));
    return raw_dphi - 2.*ROOT::Math::Pi()*region;
}

//check if a lepton can be traced back to a W boson and then to a top
bool ttbarstudy::tWfinder(const edm::Event& iEvent, const reco::GenParticle* lepton) {

    		bool ttbar=false;
    		int iStatus;

    		const reco::Candidate* iParticle = lepton->mother();

    		while(iStatus!=4){  //status=4 is the initial proton
    			iStatus = iParticle->status();
    			
    			if(abs(iParticle->pdgId())==24){ //found W
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

		if(ttbar==true){return true;}
		else{return false;}
}

bool ttbarstudy::passElectronTrig(const edm::Event& iEvent) {
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

double ttbarstudy::misIDrate(const pat::Muon* muon){

	double muonPT = muon->pt();

	for(int i=0; i<16; i++){
		if(muonPT >= binEdges[i] && muonPT < binEdges[i+1]){
			return rates[i];
		}
	}
	if(muonPT > binEdges[16]){return 0.1;}
	return 0.0;
}

void ttbarstudy::csvTable(double genMuonPt, double genElectronPt, const pat::Muon* muon, const pat::Electron* electron, const pat::Jet* bjet1, const pat::Jet* bjet2, const pat::Jet* jet1, const pat::Jet* jet2, math::XYZTLorentzVector combinedJets, const pat::MET Met, double weight) {

std::ofstream myfile;
myfile.open("neuralNetDataWW1.csv",std::ios_base::app);
myfile << muon->pt() << ", "
		 << muon->phi() << ", "
       << muon->eta() << ", "
       << electron->pt() << ", "
       << electron->phi() << ", "
       << electron->eta() << ", "
       << bjet1->pt() << ", "
       << bjet1->phi() << ", "
       << bjet1->eta() << ", "
       << bjet2->pt() << ", "
       << bjet2->phi() << ", "
       << bjet2->eta() << ", "
       << jet1->pt() << ", "
       << jet1->phi() << ", "
       << jet1->eta() << ", "
       << jet2->pt() << ", "
       << jet2->phi() << ", "
       << jet2->eta() << ", "
       << combinedJets.pt() <<", "
       << combinedJets.phi() <<", "
       << combinedJets.eta() <<", "
       << combinedJets.mass() <<", "
       << Met.pt() << ", "
       << Met.phi() <<", "
       << genElectronPt << ", "
       << genMuonPt << ", "
       << weight << ", "
       << genMuonPt/genElectronPt << "\n ";

myfile.close();

}

//no gen info here
void ttbarstudy::csvTable(const pat::Muon* muon, const pat::Electron* electron, const pat::Jet* bjet1, const pat::Jet* bjet2, const pat::Jet* jet1, const pat::Jet* jet2, math::XYZTLorentzVector combinedJets, const pat::MET Met, double weight) {

std::ofstream myfile;
myfile.open("neuralNetDataTT1.csv",std::ios_base::app);
myfile << muon->pt() << ", "
		 << muon->phi() << ", "
       << muon->eta() << ", "
       << electron->pt() << ", "
       << electron->phi() << ", "
       << electron->eta() << ", "
       << bjet1->pt() << ", "
       << bjet1->phi() << ", "
       << bjet1->eta() << ", "
       << bjet2->pt() << ", "
       << bjet2->phi() << ", "
       << bjet2->eta() << ", "
       << jet1->pt() << ", "
       << jet1->phi() << ", "
       << jet1->eta() << ", "
       << jet2->pt() << ", "
       << jet2->phi() << ", "
       << jet2->eta() << ", "
       << combinedJets.pt() <<", "
       << combinedJets.phi() <<", "
       << combinedJets.eta() <<", "
       << combinedJets.mass() <<", "
       << Met.pt() << ", "
       << Met.phi() <<", "
       << weight << "\n ";

myfile.close();

}


// ------------ method called once each job just before starting event loop  ------------
void
ttbarstudy::beginJob() {

	std::ofstream myfile;

	// myfile.open("neuralNetDataZZ1.csv",std::ios_base::app);
	// myfile<<"muon pt, muon phi, muon eta, electron pt, electron phi, electron eta, bjet 1 pt, bjet 1 phi, bjet 1 eta, bjet 2 pt, bjet 2 phi, bjet 2 eta, jet 1 pt, jet 1 phi, jet 1 eta, jet 2 pt, jet 2 phi, jet 2 eta, combined jets pt, combined jets phi, combined jets eta, combined jets mass, MET pt, MET phi, gen electron pt, gen muon pt, gen muon/electron pt ratio, event weight\n";
	// myfile.close();

	edm::Service<TFileService> fs;

	TFileDirectory countFolder = fs->mkdir("event_count");
	
	m_histoMaker.book(fs->mkdir("Analysis"));


	m_eventsWeight = {countFolder.make<TH1D>("eventsWeight","number of events weighted", 1, 0.0, 1)};
}


// ------------ method called once each job just after ending the event loop  ------------
void
ttbarstudy::endJob() {
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttbarstudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(ttbarstudy);
