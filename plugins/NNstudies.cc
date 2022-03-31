
// -*- C++ -*-
//
// Package:    ExoAnalysis/NNstudies
// Class:      NNstudies
//
/**\class NNstudies NNstudies.cc ExoAnalysis/NNstudies/plugins/NNstudies.cc
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
#include "Math/GenVector/Boost.h"

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

class NNstudies : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit NNstudies(const edm::ParameterSet&);
		~NNstudies();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		double dR2(double eta1, double eta2, double phi1, double phi2);
		double dPhi(double phi1, double phi2);
		bool tWfinder(const edm::Event&, const reco::GenParticle* );
		bool tfinder(const edm::Event&, const reco::GenParticle* );
		bool passElectronTrig(const edm::Event&);
		void csvTable(const reco::GenParticle*, const reco::GenParticle*, const reco::GenParticle*, const reco::GenParticle*, const reco::GenParticle*, const reco::GenParticle*, int binNumber, const pat::Muon*, const pat::Electron*, const pat::Jet*, const pat::Jet*, math::XYZTLorentzVector combinedJets, const pat::MET, double weight);
		int binNumber(const pat::Muon*);
		
		
		cutFlowHistos m_histoMaker;

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
NNstudies::NNstudies(const edm::ParameterSet& iConfig)
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


NNstudies::~NNstudies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
NNstudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	//bool electronTrigger=false;  		  
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


// gen lepton info

const reco::GenParticle* genMuon = 0;
const reco::GenParticle* genElectron = 0;
const reco::GenParticle* muNu = 0;
const reco::GenParticle* eNu = 0;

const reco::GenParticle* bquark = 0;
const reco::GenParticle* antibquark = 0;
const reco::GenParticle* tquark = 0;
const reco::GenParticle* antitquark = 0;

int leptonCount = 0;

for (std::vector<reco::GenParticle>::const_iterator iParticle = genParticles->begin(); iParticle != genParticles->end(); iParticle++) {
	if( ! iParticle->isHardProcess() ){ continue; }
	if(iParticle->pdgId()==6 && tquark==0){ tquark =  &(*(iParticle)); }
	if(iParticle->pdgId()==-6 && antitquark==0){ antitquark =  &(*(iParticle)); }
	if(tfinder(iEvent, &(*iParticle))){ // check if the gen particle comes from a top
		if(iParticle->pdgId()==5){ bquark = &(*(iParticle)); }
		if(iParticle->pdgId()==-5){ antibquark = &(*(iParticle)); }

		if(tWfinder(iEvent, &(*iParticle))){  // check if the gen particle comes from a top->W->lepton
			if(abs(iParticle->pdgId())==13){
				leptonCount += 1;
				if(genMuon==0){genMuon = &(*(iParticle));}
			}	
			else if(abs(iParticle->pdgId())==11){
				leptonCount += 1;
				if(genElectron==0){genElectron = &(*(iParticle));}
			}
			else if(abs(iParticle->pdgId())==12){
				if(eNu==0){eNu = &(*(iParticle));}
			}
			else if(abs(iParticle->pdgId())==14){
				if(muNu==0){muNu = &(*(iParticle));}
			}

		}
	}
}


//start reco

	const pat::MET Met = recoMET->front();

	math::XYZTLorentzVector combinedJetsP4 = {0., 0., 0., 0.};

	edm::Handle<std::vector<pat::Jet>> recoJetsAK4;  
	iEvent.getByToken(m_AK4recoCHSJetsToken, recoJetsAK4);  

	const pat::Jet* bJet1 = 0;
	// // const pat::Jet* bJet2=0;

	const pat::Jet* Jet1 = 0;
	// //const pat::Jet* Jet2=0;

	const pat::Jet* bJet=0;
	const pat::Jet* antibJet=0;

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
			if(BJP > 0.4184) {
				if(bJet == 0 && sqrt(dR2(iJet->eta(), bquark->eta(), iJet->phi(), bquark->phi())) < 0.3 ){
					bJet=&(*(iJet));
				}
				else if(antibJet == 0 && sqrt(dR2(iJet->eta(), antibquark->eta(), iJet->phi(), antibquark->phi())) < 0.3 ){
					antibJet=&(*(iJet));
				}
			}
			// if(BJP > 0.4184){ 
			// 	if(bJet1==0){
			// 		bJet1=&(*(iJet));
			// 	}
			// 	else if(Jet1==0){
			// 		Jet1=&(*(iJet));
			// 	}
			// }

			// if(BJP < 0.4184){
			// 	if(Jet1==0){
			// 		Jet1=&(*(iJet));
			// 	}
			// }
		}


// check electron trigger
//if (passElectronTrig(iEvent)){ electronTrigger=true; }

//muon/electron reconstruction

	//muon reco

	  const pat::Muon* recoMuon1=0;

   	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){

   		if(!(iMuon->isHighPtMuon(*myEvent.PVertex))) continue; // || !iMuon->passed(reco::Muon::TkIsoTight)) continue; //preliminary cut

   		if(recoMuon1==0){
   				recoMuon1=&(*(iMuon));
   		}
   		
   		// leptonCount += 1;	
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

				// leptonCount += 1;
				
			}

		m_eventsWeight->Fill(0.5, eventCount);

		if(leptonCount == 2){

			if(genMuon!=0 && genElectron!=0 && antibJet!=0 && bJet!=0){
				if(genMuon->pdgId()>0){bJet1=bJet; Jet1=antibJet; }
				if(genMuon->pdgId()<0){bJet1=antibJet; Jet1=bJet;}

				if(muNu!=0 && eNu!=0 && recoMuon1!=0 && recoElectron1!=0 && tquark
					!= 0 && antitquark != 0){	
					csvTable(genMuon,genElectron,muNu,eNu,tquark,antitquark,binNumber(recoMuon1),recoMuon1,recoElectron1,bJet1,Jet1,combinedJetsP4,Met,eventCount);
					// ROOT::Math::Boost boostTT;
				  // boostTT.SetComponents(combinedJetsP4.BoostToCM());
				  // auto tquarkP4 = boostTT(tquark->p4());
				  // auto antitquarkP4 = boostTT(antitquark->p4());
					// m_histoMaker.fill( (tquarkP4 + antitquarkP4 ).pt(), genMuon->pt());
					m_histoMaker.fill( (tquark->p4() + antitquark->p4() ).pt(), genMuon->pt());
				}
			}


		}
	}




//HELPERS
double NNstudies::dR2(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    double dphi = dPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
}
double NNstudies::dPhi(double phi1, double phi2) {
    double raw_dphi = phi1 - phi2;
    if (fabs(raw_dphi) < ROOT::Math::Pi()) return raw_dphi;
    double region = std::round(raw_dphi / (2.*ROOT::Math::Pi()));
    return raw_dphi - 2.*ROOT::Math::Pi()*region;
}

//check if a lepton can be traced back to a W boson and then to a top
bool NNstudies::tWfinder(const edm::Event& iEvent, const reco::GenParticle* lepton) {

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

bool NNstudies::tfinder(const edm::Event& iEvent, const reco::GenParticle* quark) {

    		bool ttbar=false;
    		int iStatus;

    		const reco::Candidate* iParticle = quark->mother();
    		iStatus = iParticle->status();

    		while(iStatus!=4){  //status=4 is the initial proton

    			   	if(abs(iParticle->pdgId())==6){ 
    			   		ttbar=true; 
    			   		break;
    			   	}
    			   iParticle = iParticle->mother();
    			   iStatus = iParticle->status();
    		}

		if(ttbar==true){return true;}
		else{return false;}
}

bool NNstudies::passElectronTrig(const edm::Event& iEvent) {
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

int NNstudies::binNumber(const pat::Muon* muon){

	double muonPT = muon->pt();

	for(int i=0; i<16; i++){
		if(muonPT >= binEdges[i] && muonPT < binEdges[i+1]){
			return i+1;
		}
	}
	if(muonPT > binEdges[16]){return 17;}
	return 0;
}

void NNstudies::csvTable(const reco::GenParticle* genMuon, const reco::GenParticle* genElectron, const reco::GenParticle* muNu, const reco::GenParticle* eNu, const reco::GenParticle* tquark, const reco::GenParticle* antitquark, int binNumber, const pat::Muon* muon, const pat::Electron* electron, const pat::Jet* bjet1, const pat::Jet* jet1, math::XYZTLorentzVector combinedJets, const pat::MET Met, double weight) {

ROOT::Math::Boost boostJets;
boostJets.SetComponents(combinedJets.BoostToCM());


math::XYZTLorentzVector muonP4 = muon->p4();
math::XYZTLorentzVector electronP4 = electron->p4();
math::XYZTLorentzVector bjet1P4 = bjet1->p4();
math::XYZTLorentzVector jet1P4=jet1->p4();
math::XYZTLorentzVector genMuonP4 = genMuon->p4();
math::XYZTLorentzVector genElectronP4 = genElectron->p4();
math::XYZTLorentzVector muNuP4 = muNu->p4();
math::XYZTLorentzVector eNuP4 = eNu->p4();
math::XYZTLorentzVector tquarkP4 = tquark->p4();
math::XYZTLorentzVector antitquarkP4 = antitquark->p4();


std::ofstream myfile;
myfile.open("neuralNetDataTT_3.csv",std::ios_base::app);
myfile << muonP4.Px() << ", "
		   << muonP4.Py() << ", "
       << muonP4.Pz() << ", "
       << muonP4.E()*abs(genMuon->pdgId())/genMuon->pdgId() << ", "
       << electronP4.Px() << ", "
       << electronP4.Py() << ", "
       << electronP4.Pz() << ", "
       << electronP4.E()*abs(genElectron->pdgId())/genElectron->pdgId()  << ", "
       << bjet1P4.Px() << ", "
       << bjet1P4.Py() << ", "
       << bjet1P4.Pz() << ", "
       << bjet1P4.E()  << ", "
       << jet1P4.Px() << ", "
       << jet1P4.Py() << ", "
       << jet1P4.Pz() << ", "
       << jet1P4.E()  << ", "
       << combinedJets.Px() <<", "
       << combinedJets.Py() <<", "
       << combinedJets.Pz() <<", "
       << combinedJets.E()  <<", "
       << Met.pt() << ", "
       << Met.phi() <<", "
       << weight << ", "
       << binNumber << ", "
       << genMuonP4.Px() << ", "
       << genMuonP4.Py() << ", "
       << genMuonP4.Pz() << ", "
       << genMuonP4.E()  << ", "
       << genElectronP4.Px() << ", "
       << genElectronP4.Py() << ", "
       << genElectronP4.Pz() << ", "
       << genElectronP4.E()  << ", "
       << muNuP4.Px() << ", "
       << muNuP4.Py() << ", "
       << muNuP4.Pz() << ", "
       << muNuP4.E()  << ", "
       << eNuP4.Px() << ", "
       << eNuP4.Py() << ", "
       << eNuP4.Pz() << ", "
       << eNuP4.E()  << ", "
       << antitquarkP4.Px() << ", "
       << antitquarkP4.Py() << ", "
       << antitquarkP4.Pz() << ", "
       << antitquarkP4.E() << ", "
       << tquarkP4.Px() << ", "
       << tquarkP4.Py() << ", "
       << tquarkP4.Pz() << ", "
       << tquarkP4.E() << "\n ";


myfile.close();

}


// ------------ method called once each job just before starting event loop  ------------
void
NNstudies::beginJob() {

	std::ofstream myfile;

	// myfile.open("neuralNetDataTT_1.csv",std::ios_base::app);
	// myfile<<"muonP1,muonP2,muonP3,muonP4,electronP1,electronP2,electronP3,electronP4,bjetP1,bjetP2,bjetP3,bjetP4,jetP1,jetP2,jetP3,jetP4,combinedJetsP1,combinedJetsP2,combinedJetsP3,combinedJetsP4,METpt,METphi,eventWeight,binNumber,genMuonP1,genMuonP2,genMuonP3,genMuonP4,genElectronP1,genElectronP2,genElectronP3,genElectronP4,muNuP1,muNuP2,muNuP3,muNuP4,eNuP1,eNuP2,eNuP3,eNuP4,antitquarkP1,antitquarkP2,antitquarkP3,antitquarkP4,tquarkP1,tquarkP2,tquarkP3,tquarkP4\n";
	// myfile.close();

	edm::Service<TFileService> fs;

	TFileDirectory countFolder = fs->mkdir("event_count");
	
	m_histoMaker.book(fs->mkdir("Analysis"));  //2jets

	m_eventsWeight = {countFolder.make<TH1D>("eventsWeight","number of events weighted", 1, 0.0, 1)};

	
}


// ------------ method called once each job just after ending the event loop  ------------
void
NNstudies::endJob() {
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NNstudies::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(NNstudies);
