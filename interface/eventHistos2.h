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
	void book(TFileDirectory histoFolder, int nDirectory);
	void fill(eventBits2& event);
	//eventHistos2() : m_histoFolder(5) {}



	private:
	

	//std::vector<TFileDirectory> m_histoFolder = std::vector<TFileDirectory>(5);
	TFileDirectory m_histoFolder;
	//edm::Service<TFileService> fs; 

	//General stats histos
	
	
	//1D

	std::vector<TH1D*> m_matchStatus = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_unmatchedRecoPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_unmatchedRecoEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_unmatchedRecoPt = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_unmatchedGenPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_unmatchedGenEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_unmatchedGenPt = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_unmatchedDR = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_numRecoMuons = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_numRecoElectrons = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_numGenMuons = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_numGenElectrons = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_genMinusRecoMuons = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genMinusRecoElectrons = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_genRecoMuonPtRatio = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genMuonPt = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genMuonPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genMuonEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoMuonPt = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoMuonPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoMuonEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignGenRecoMuonPtRatio = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignGenRecoMuonPtRatio = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoMuonPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoMuonPhi = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoMuonEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoMuonEta = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoMuonPt = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoMuonPt = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_bestMuondR = std::vector<TH1D*>(5);
	std::vector<TH1D*> m_TWgenMuon = std::vector<TH1D*>(5);

	std::vector<TH1D*> m_genRecoElectronPtRatio		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genElectronPt		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genElectronPhi		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_genElectronEta		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoElectronPt		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoElectronPhi		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_recoElectronEta		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignGenRecoElectronPtRatio	=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignGenRecoElectronPtRatio	=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoElectronPhi		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoElectronPhi		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoElectronEta		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoElectronEta		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_WrongSignRecoElectronPt		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_RightSignRecoElectronPt		=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_bestElectrondR				=std::vector<TH1D*>(5);
	std::vector<TH1D*> m_TWgenElectron = std::vector<TH1D*>(5);

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

	int nDirectory;

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