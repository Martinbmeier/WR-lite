#ifndef eventHistos2_h
#define eventHistos2_h


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

//ROOT CLASSES
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
//C++ CLASSES
#include <iostream>
//LOCAL CLASSES
#include "eventHistos2.h"
#include "eventBits2.h"


class eventHistos2 {
	public:
	eventHistos2();
	//void book(TFileDirectory histoFolder);
	void book(TFileDirectory histoFolder, int nCut, int genType);
	void fill(eventBits2& event, int cutNumber);
	//eventHistos2() : m_histoFolder(5) {}



	private:
	

	//std::vector<TFileDirectory> m_histoFolder = std::vector<TFileDirectory>(5);
	TFileDirectory m_histoFolder;
	//edm::Service<TFileService> fs; 

	//General stats histos
	
	
	//1D

	std::vector<std::vector<TH1D*>> m_matchStatus{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_unmatchedRecoPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_unmatchedRecoEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_unmatchedRecoPt{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_unmatchedGenPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_unmatchedGenEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_unmatchedGenPt{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_unmatchedRecoDR{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_numRecoMuons{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_numRecoElectrons{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_numGenMuons{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_numGenElectrons{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_genMinusRecoMuons{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genMinusRecoElectrons{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH1D*>> m_genRecoMuonPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genMuonPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genMuonPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genMuonEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoMuonPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoMuonPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoMuonEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignGenRecoMuonPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignGenRecoMuonPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoMuonPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoMuonPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoMuonEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoMuonEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoMuonPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoMuonPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_bestMuondR{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_TWgenMuon{4,std::vector<TH1D*>(4,0)};


	std::vector<std::vector<TH1D*>> m_genRecoElectronPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genElectronPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genElectronPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_genElectronEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoElectronPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoElectronPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_recoElectronEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignGenRecoElectronPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignGenRecoElectronPtRatio{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoElectronPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoElectronPhi{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoElectronEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoElectronEta{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_WrongSignRecoElectronPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_RightSignRecoElectronPt{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_bestElectrondR{4,std::vector<TH1D*>(4,0)};
	std::vector<std::vector<TH1D*>> m_TWgenElectron{4,std::vector<TH1D*>(4,0)};

	std::vector<std::vector<TH2D*>> m_phiEta{4,std::vector<TH2D*>(4,0)};

	/*
	TH1D * m_genRecoMuonPtRatio;
	TH1D * m_genMuonPt;
	TH1D * m_genMuonPhi;
	TH1D * m_recoMuonPt;
	TH1D * m_recoMuonPhi;
	TH1D * m_recoMuonEta;
	TH1D * m_WrongSignGenRecoMuonPtRatio;
	TH1D * m_RightSignGenRecoMuonPtRatio;
	TH1D * m_WrongSignRecoMuonPhi;
	TH1D * m_RightSignRecoMuonPhi;
	TH1D * m_WrongSignRecoMuonPt;
	TH1D * m_RightSignRecoMuonPt;
	TH1D * m_bestMuondR;



	TH1D * m_genRecoElectronPtRatio;
	TH1D * m_genElectronPt;
	TH1D * m_genElectronPhi;
	TH1D * m_recoElectronPt;
	TH1D * m_recoElectronPhi;
	TH1D * m_recoElectronEta;
	TH1D * m_WrongSignGenRecoElectronPtRatio;
	TH1D * m_RightSignGenRecoElectronPtRatio;
	TH1D * m_WrongSignRecoElectronPhi;
	TH1D * m_RightSignRecoElectronPhi;
	TH1D * m_WrongSignRecoElectronPt;
	TH1D * m_RightSignRecoElectronPt;
	TH1D * m_bestElectrondR;
	*/

	//TFileDirectory * top;
	//TFileDirectory * allEventsTD;
	//TFileDirectory * twoMuonsTD ;
	//TFileDirectory * twoElectronsTD; 
	//TFileDirectory * muonElectronTD ;
	//TFileDirectory * muonTauTD ;
	//TFileDirectory * electronTauTD;


	//TH1D * m_numElectrons;
	//TH1D * m_numMuons;

	//TH1D * m_notTTbar;



   
};

#endif