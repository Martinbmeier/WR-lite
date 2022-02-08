
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
		void csvTable(double genMuonPt, int binNumber, const pat::Muon*, const pat::Electron*, const pat::Jet*, const pat::Jet*, math::XYZTLorentzVector combinedJets, const pat::MET, double weight);
		int binNumber(const pat::Muon*);
		
		
		//cutFlowHistos m_histoMaker;

		TH1D* m_eventsWeight;

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
		//double center[17] = {10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 235.0, 275.0, 325.0, 375.0, 700.0,};
		
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
void
cut_flow2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	bool electronTrigger=false;  		  
	//bool angularSeparation=false; 

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
	// double eventWeight = eventInfo->weight();
	
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

	const pat::Jet* bJet1 = 0;
	// const pat::Jet* bJet2=0;

	const pat::Jet* Jet1 = 0;
	//const pat::Jet* Jet2=0;


	//Get jets with maximum pt
		for(std::vector<pat::Jet>::const_iterator iJet = recoJetsAK4->begin(); iJet != recoJetsAK4->end(); iJet++) {

			combinedJetsP4=combinedJetsP4+iJet->p4();

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

			if(BJP > 0.4184){ 
				if(bJet1==0){
					bJet1=&(*(iJet));
				}
				else if(Jet1==0){
					Jet1=&(*(iJet));
				}
			}

			if(BJP < 0.4184){
				if(Jet1==0){
					Jet1=&(*(iJet));
				}
			}
		}



// gen lepton info

double genMuonpT=-1000;
	
for (std::vector<reco::GenParticle>::const_iterator iParticle = genParticles->begin(); iParticle != genParticles->end(); iParticle++) {
	if( ! iParticle->isHardProcess() ){ continue; }
	//if( ! tWfinder(iEvent, &(*iParticle))){ continue; }  //could check if the gen particle comes from a top->W->lepton
	if(abs(iParticle->pdgId())==13){
		if(genMuonpT<0){genMuonpT=iParticle->pt();}	
	}	
}


// check electron trigger
if (passElectronTrig(iEvent)){ electronTrigger=true; }

//muon/electron reconstruction

	int leptonCount = 0;


	//muon reco

	  const pat::Muon* recoMuon1=0;

   	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){

   		if(!(iMuon->isHighPtMuon(*myEvent.PVertex))) continue; // || !iMuon->passed(reco::Muon::TkIsoTight)) continue; //preliminary cut

   		if(recoMuon1==0){
   				recoMuon1=&(*(iMuon));
   		}
   		
   		leptonCount += 1;	
		}

	//electron reco

		const pat::Electron* recoElectron1=0;

			for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++){	

				const vid::CutFlowResult* vidResult =  iElectron->userData<vid::CutFlowResult>("heepElectronID_HEEPV70");
				const bool heepIDVID = vidResult->cutFlowPassed();
				if(heepIDVID == false){continue;}
				if(iElectron->pt()<35){continue;}
				
				if(recoElectron1==0){ 
						recoElectron1=&(*(iElectron)); 
				}

				leptonCount+=1;
				
			}

		m_eventsWeight->Fill(0.5, eventCount);

		if(leptonCount == 2 && electronTrigger){

				if(recoMuon1!=0 && recoElectron1!=0 && bJet1!=0 && Jet1!=0){
					double dileptonSeparation=sqrt(dR2(recoMuon1->eta(), recoElectron1->eta(), recoMuon1->phi(), recoElectron1->phi()));
			   	double muonJet1Sep=sqrt(dR2(bJet1->eta(), recoMuon1->eta(), bJet1->phi(), recoMuon1->phi()));
			   	double muonJet2Sep=sqrt(dR2(Jet1->eta(), recoMuon1->eta(), Jet1->phi(), recoMuon1->phi()));
			   	double electronJet1Sep=sqrt(dR2(bJet1->eta(), recoElectron1->eta(), bJet1->phi(), recoElectron1->phi()));
					double electronJet2Sep=sqrt(dR2(Jet1->eta(), recoElectron1->eta(), Jet1->phi(), recoElectron1->phi()));
					double jetSeparation=sqrt(dR2(bJet1->eta(), Jet1->eta(), bJet1->phi(), Jet1->phi()));		

	
				if(dileptonSeparation>0.4 && muonJet1Sep>0.4 && muonJet2Sep>0.4 && electronJet1Sep > 0.4 && electronJet2Sep>0.4 && jetSeparation>0.4){	
					csvTable(genMuonpT,binNumber(recoMuon1),recoMuon1,recoElectron1,bJet1,Jet1,combinedJetsP4,Met,eventCount);	
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

int cut_flow2::binNumber(const pat::Muon* muon){

	double muonPT = muon->pt();

	for(int i=0; i<16; i++){
		if(muonPT >= binEdges[i] && muonPT < binEdges[i+1]){
			return i+1;
		}
	}
	if(muonPT > binEdges[16]){return 18;}
	return 0;
}

void cut_flow2::csvTable(double genMuonPt, int binNumber, const pat::Muon* muon, const pat::Electron* electron, const pat::Jet* bjet1, const pat::Jet* jet1, math::XYZTLorentzVector combinedJets, const pat::MET Met, double weight) {

std::ofstream myfile;
myfile.open("neuralNetDataTT2.csv",std::ios_base::app);
myfile << muon->pt() << ", "
		 << muon->phi() << ", "
       << muon->eta() << ", "
       << electron->pt() << ", "
       << electron->phi() << ", "
       << electron->eta() << ", "
       << bjet1->pt() << ", "
       << bjet1->phi() << ", "
       << bjet1->eta() << ", "
       << jet1->pt() << ", "
       << jet1->phi() << ", "
       << jet1->eta() << ", "
       << combinedJets.pt() <<", "
       << combinedJets.phi() <<", "
       << combinedJets.eta() <<", "
       << combinedJets.mass() <<", "
       << Met.pt() << ", "
       << Met.phi() <<", "
       << weight << ", "
       << binNumber << ", "
       << genMuonPt << "\n ";


myfile.close();

}


// ------------ method called once each job just before starting event loop  ------------
void
cut_flow2::beginJob() {

	std::ofstream myfile;

	// myfile.open("neuralNetDataTT1.csv",std::ios_base::app);
	// myfile<<"muon pt, muon phi, muon eta, electron pt, electron phi, electron eta, bjet pt, bjet phi, bjet eta, jet pt, jet phi, jet eta, combined jets pt, combined jets phi, combined jets eta, combined jets mass, MET pt, MET phi, event weight, bin number, gen muon pt\n";
	// myfile.close();

	edm::Service<TFileService> fs;

	TFileDirectory countFolder = fs->mkdir("event_count");
	
	//m_histoMaker.book(fs->mkdir("Analysis"));  //2jets

	m_eventsWeight = {countFolder.make<TH1D>("eventsWeight","number of events weighted", 1, 0.0, 1)};

	
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
