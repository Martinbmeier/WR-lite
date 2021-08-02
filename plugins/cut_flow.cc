
// -*- C++ -*-
//
// Package:    ExoAnalysis/cut_flow
// Class:      cut_flow
//
/**\class cut_flow cut_flow.cc ExoAnalysis/cut_flow/plugins/cut_flow.cc
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
#include "DataFormats/Candidate/interface/Candidate.h"
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

class cut_flow : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit cut_flow(const edm::ParameterSet&);
		~cut_flow();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		double dR2(double eta1, double eta2, double phi1, double phi2);
		double dPhi(double phi1, double phi2);
		bool tWfinder(const edm::Event&, const reco::GenParticle* );
		bool passElectronTrig(const edm::Event&);
		//double transverseSphericity(math::XYZTLorentzVector p1, math::XYZTLorentzVector p2, math::XYZTLorentzVector p3);
		//void saveElectronData(eventBits2 * iBit, double matched1Mass, double matched2Mass);
		//void saveMuonData(eventBits2 * iBit, double matched1Mass, double matched2Mass);
		
		
		cutFlowHistos m_histoMaker;
		TH1D m_eventsWeight;


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
cut_flow::cut_flow(const edm::ParameterSet& iConfig)
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


cut_flow::~cut_flow()
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
cut_flow::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	bool electronTrigger=false;  //cut1
	bool oneHeepElectron=false;  //cut2
	bool oneMuonHighpT=false;    //cut3
	bool twoJets=false;				//cut4
	bool angularSeparation=false;//cut5
	bool dileptonMass=false;		//cut6
	bool oneBTag=false;				//cut7
	bool twoBTag=false;				//cut8
	bool muonIsolation1=false;	//cut9
	bool muonIsolation2=false;	//cut10

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
	
  
	float eventCount = eventInfo->weight()/fabs(eventInfo->weight());
	double eventWeight = eventInfo->weight();
	
	
	edm::Handle<std::vector<reco::Vertex>> vertices;
	iEvent.getByToken(m_offlineVerticesToken, vertices);
	if(!vertices.isValid()) {
		throw cms::Exception("Vertex collection not valid!");
	}

	iBit.hasPVertex = myEvent.PVselection(vertices);






//start reco

	edm::Handle<std::vector<pat::Jet>> recoJetsAK4;  
	iEvent.getByToken(m_AK4recoCHSJetsToken, recoJetsAK4);  


	//const pat::Jet* leadJet = 0;
	//const pat::Jet* subleadJet = 0;
	int jetCount = 0;

	int btagcount = 0;

	//bool btagged;

	const pat::Jet* Jet1;
	const pat::Jet* Jet2;


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
			if(BJP > 0.4184){ btagcount++; }		
			if(jetCount==0){Jet1=&(*iJet);}
			if(jetCount>0){Jet2=&(*iJet);}
			jetCount++;
		}

		if(jetCount>1){twoJets=true;}
		if(btagcount>0){oneBTag=true;}
		if(btagcount>1){twoBTag=true;}




	const pat::Muon* leadMuon;

	const pat::Electron* leadElectron;




//cut flow things

//check electron trigger
if (passElectronTrig(iEvent)){ electronTrigger=true; }

//muon/electron reconstruction

	//muon reco

		double leadMuonpT = -1000;
		double newLeadMuonpT;
		//int muCount=0;

   	for(std::vector<pat::Muon>::const_iterator iMuon = highMuons->begin(); iMuon != highMuons->end(); iMuon++){

   		//if(fabs(iMuon->eta()) > 2.4 || iMuon->tunePMuonBestTrack()->pt() < 10 || !(iMuon->isHighPtMuon(*myEvent.PVertex)) || (iMuon->isolationR03().sumPt/iMuon->pt() > .1)){ continue;} //preliminary cut
   		if(iMuon->isHighPtMuon(*myEvent.PVertex)){oneMuonHighpT=true;}

   		
   		newLeadMuonpT=iMuon->pt();

   		if(newLeadMuonpT>leadMuonpT){leadMuonpT=newLeadMuonpT; leadMuon=&(*iMuon);}

		}

		if(leadMuon->passed(reco::Muon::TkIsoLoose)){muonIsolation1=true;}
		if(leadMuon->passed(reco::Muon::TkIsoTight)){muonIsolation1=true; muonIsolation2=true; }
		

	//electron reco

		double leadElectronpT = -1000;
		double newLeadElectronpT;

			//for all reco electrons, loop through gen electrons to find spatial matches
			for(std::vector<pat::Electron>::const_iterator iElectron = highElectrons->begin(); iElectron != highElectrons->end(); iElectron++){	
				//if(fabs(iElectron->eta()) > 2.4) {continue;}
				//if(iElectron->pt() < 10 ) {continue;}
				oneHeepElectron=true;

				newLeadElectronpT=iElectron->pt();
   			//newLeadElectronp4=iElectron->p4();

   			if(newLeadElectronpT>leadElectronpT){leadElectronpT=newLeadElectronpT; leadElectron=&(*iElectron);}

			}


			if(oneHeepElectron && oneMuonHighpT && twoJets){
				double dileptonSeparation=sqrt(dR2(leadMuon->eta(), leadElectron->eta(), leadMuon->phi(), leadElectron->phi()));
			   double muonJet1Sep=sqrt(dR2(Jet1->eta(), leadMuon->eta(), Jet1->phi(), leadMuon->phi()));
			   double muonJet2Sep=sqrt(dR2(Jet2->eta(), leadMuon->eta(), Jet2->phi(), leadMuon->phi()));
			   double electronJet1Sep=sqrt(dR2(Jet1->eta(), leadElectron->eta(), Jet1->phi(), leadElectron->phi()));
				double electronJet2Sep=sqrt(dR2(Jet2->eta(), leadElectron->eta(), Jet2->phi(), leadElectron->phi()));
				double jetSeparation=sqrt(dR2(Jet2->eta(), Jet1->eta(), Jet2->phi(), Jet1->phi()));
				if(dileptonSeparation>0.4 && muonJet1Sep>0.4 && muonJet2Sep>0.4 && electronJet1Sep > 0.4 && electronJet2Sep>0.4 && jetSeparation>0.4){angularSeparation=true;} //check for lepton separation
			} 
			if((leadElectron->p4()+leadMuon->p4()).mass()>150){dileptonMass=true;} //check for dilepton mass

		
	m_eventsWeight.Fill(0.5, eventCount);

	m_histoMaker.fill(leadMuonpT,0,eventWeight);

	if(electronTrigger){
		m_histoMaker.fill(leadMuonpT,1,eventWeight);
		if(oneHeepElectron){
			m_histoMaker.fill(leadMuonpT,2,eventWeight);
			if(oneMuonHighpT){
				m_histoMaker.fill(leadMuonpT,3,eventWeight);
				if(twoJets){
					m_histoMaker.fill(leadMuonpT,4,eventWeight);
					if(angularSeparation){
						m_histoMaker.fill(leadMuonpT,5,eventWeight);
						if(dileptonMass){
							m_histoMaker.fill(leadMuonpT,6,eventWeight);
							if(oneBTag){
								m_histoMaker.fill(leadMuonpT,7,eventWeight);
								if(twoBTag){
									m_histoMaker.fill(leadMuonpT,8,eventWeight);
									if(muonIsolation1){
										m_histoMaker.fill(leadMuonpT,9,eventWeight);
										if(muonIsolation2)
											m_histoMaker.fill(leadMuonpT,10,eventWeight);
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
double cut_flow::dR2(double eta1, double eta2, double phi1, double phi2) {
    double deta = eta1 - eta2;
    double dphi = dPhi(phi1, phi2);
    return deta*deta + dphi*dphi;
}
double cut_flow::dPhi(double phi1, double phi2) {
    double raw_dphi = phi1 - phi2;
    if (fabs(raw_dphi) < ROOT::Math::Pi()) return raw_dphi;
    double region = std::round(raw_dphi / (2.*ROOT::Math::Pi()));
    return raw_dphi - 2.*ROOT::Math::Pi()*region;
}

//check if a lepton can be traced back to a W boson and then to a top
bool cut_flow::tWfinder(const edm::Event& iEvent, const reco::GenParticle* lepton) {

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

bool cut_flow::passElectronTrig(const edm::Event& iEvent) {
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


// ------------ method called once each job just before starting event loop  ------------
void
cut_flow::beginJob() {

	edm::Service<TFileService> fs; 

	
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

	m_eventsWeight = {m_histoFolder.make<TH1D>("eventsWeight","number of events weighted", 1, 0.0, 1)};



	
}


// ------------ method called once each job just after ending the event loop  ------------
void
cut_flow::endJob() {
	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
cut_flow::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(cut_flow);
