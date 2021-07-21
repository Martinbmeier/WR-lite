//C++ CLASSES
#include <iostream>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
//LOCAL CLASSES
#include "ExoAnalysis/WR_lite/interface/eventBits2.h"
#include "ExoAnalysis/WR_lite/interface/eventHistos2.h"

#include "TH1D.h"
#include "TStyle.h"



eventHistos2::eventHistos2 () {}
  
//void eventHistos2::book(TFileDirectory histoFolder) {
void eventHistos2::book(TFileDirectory histoFolder, int nCut, int genType) {

	gStyle->SetOptStat("omen");
	m_histoFolder=histoFolder;

	//Unmatchded lepton histos
	
	m_matchStatus[nCut][genType] = {m_histoFolder.make<TH1D>("matchStatus","were there successful gen-reco matches in this event?", 2,1,2)};
	m_matchStatus[nCut][genType]->GetXaxis()->SetBinLabel(1,"matched");
	m_matchStatus[nCut][genType]->GetXaxis()->SetBinLabel(2,"unmatched");

	m_numRecoMuons[nCut][genType] = {m_histoFolder.make<TH1D>("numRecoMuons","num reco", 10,0,10)};
	m_numRecoElectrons[nCut][genType] = {m_histoFolder.make<TH1D>("numRecoElectrons","num reco", 10,0,10)};

	m_numGenMuons[nCut][genType] = {m_histoFolder.make<TH1D>("numGenMuons","num gen", 10,0,10)};
	m_numGenElectrons[nCut][genType] = {m_histoFolder.make<TH1D>("numGenElectrons","num gen", 10,0,10)};

	m_genMinusRecoMuons[nCut][genType] = {m_histoFolder.make<TH1D>("genMinusRecoMuons","gen-reco", 10,-5,5)};
	m_genMinusRecoElectrons[nCut][genType] = {m_histoFolder.make<TH1D>("genMinusRecoElectrons","gen-reco", 10,-5,5)};

	m_unmatchedRecoPt[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedRecoPt","Pt of unmatched lepton", 100,0,500)};
	m_unmatchedRecoPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	m_unmatchedRecoEta[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedRecoEta","Eta of unmatched lepton", 100,-3.4,3.4)};
	m_unmatchedRecoEta[nCut][genType]->GetXaxis()-> SetTitle("eta (radians)");
	m_unmatchedRecoPhi[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedRecoPhi","Phi of unmatched lepton", 100,-3.4,3.4)};
	m_unmatchedRecoPhi[nCut][genType]->GetXaxis()-> SetTitle("phi (radians)");

	m_unmatchedGenPt[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedGenPt","Pt of unmatched gen lepton", 100,0,500)};
	m_unmatchedGenPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	m_unmatchedGenEta[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedGenEta","Eta of unmatched gen lepton", 100,-3.4,3.4)};
	m_unmatchedGenEta[nCut][genType]->GetXaxis()-> SetTitle("eta (radians)");
	m_unmatchedGenPhi[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedGenPhi","Phi of unmatched gen lepton", 100,-3.4,3.4)};
	m_unmatchedGenPhi[nCut][genType]->GetXaxis()-> SetTitle("phi (radians)");

	m_unmatchedRecoDR[nCut][genType] = {m_histoFolder.make<TH1D>("unmatchedDR","smallest delta R of unmatched lepton", 100,0,5)};

	gStyle->SetOptStat("omen");

	//Muon Histos

	m_genRecoMuonPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("genRecoMuonPtRatio" , "ratio of gen pt to reco pt for muons" , 50 , 0 , 5 )};
	gStyle->SetOptStat("omen");

	m_genMuonPt[nCut][genType] = {m_histoFolder.make<TH1D>("genMuonPt","Pt for gen muons",100,0,500)};
	m_genMuonPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_genMuonPhi[nCut][genType] = {m_histoFolder.make<TH1D>("genMuonPhi","Phi for gen muons",100,-3.4,3.4)};
	m_genMuonPhi[nCut][genType]->GetXaxis()-> SetTitle("Phi");

	m_genMuonEta[nCut][genType] = {m_histoFolder.make<TH1D>("genMuonEta","Eta for gen muons",100,-3.4,3.4)};
	m_genMuonEta[nCut][genType]->GetXaxis()-> SetTitle("Eta");

	m_recoMuonPt[nCut][genType] =  {m_histoFolder.make<TH1D>("recoMuonPt","Pt for reco muons",100,0,500)};
	m_recoMuonPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");

	m_recoMuonPhi[nCut][genType] = {m_histoFolder.make<TH1D>("recoMuonPhi","Phi for reco muons",100,-3.4,3.4)};
	m_recoMuonPhi[nCut][genType]->GetXaxis()-> SetTitle("Phi");

	m_recoMuonEta[nCut][genType] = {m_histoFolder.make<TH1D>("recoMuonEta","Eta for reco muons",100,-3.4,3.4)};
	m_recoMuonEta[nCut][genType]->GetXaxis()-> SetTitle("Eta");
	gStyle->SetOptStat("omen");
	
	m_WrongSignGenRecoMuonPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignGenRecoMuonPtRatio" , "ratio of gen pt to reco pt (wrong charge matched muons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_WrongSignRecoMuonPt[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonPt" , "Pt for wrong sign reco muons" , 100 , 0 , 500 )};
	m_WrongSignRecoMuonPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_WrongSignRecoMuonPhi[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonPhi" , "phi for wrong sign muons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoMuonPhi[nCut][genType]->GetXaxis()->SetTitle("phi");

	m_WrongSignRecoMuonEta[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoMuonEta" , "Eta for wrong sign muons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoMuonEta[nCut][genType]->GetXaxis()->SetTitle("eta");

	m_RightSignGenRecoMuonPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignGenRecoMuonPtRatio" , "ratio of gen pt to reco pt (right charge matched muons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_RightSignRecoMuonPt[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoMuonPt" , "Pt for right sign reco muons" , 100 , 0 , 500 )};
	m_RightSignRecoMuonPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	
	m_RightSignRecoMuonPhi[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoMuonPhi" , "phi for right sign muons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoMuonPhi[nCut][genType]->GetXaxis()->SetTitle("phi");

	m_RightSignRecoMuonEta[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoMuonEta" , "eta for right sign muons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoMuonEta[nCut][genType]->GetXaxis()->SetTitle("eta");

	
	m_bestMuondR[nCut][genType] = {m_histoFolder.make<TH1D>("bestMuondR" , "smallest deltaR (muons)" , 50 , 0 , 0.03 )};
	m_bestMuondR[nCut][genType]->GetXaxis()-> SetTitle("dR");
	gStyle->SetOptStat("omen");




	//electron histograms
	m_genRecoElectronPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("genRecoElectronPtRatio" , "ratio of gen pt to reco pt for electrons" , 50 , 0 , 5 )};


	m_genElectronPt[nCut][genType] = {m_histoFolder.make<TH1D>("genElectronPt","Pt for gen electrons",100,0,500)};
	m_genElectronPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_genElectronPhi[nCut][genType] = {m_histoFolder.make<TH1D>("genElectronPhi","Phi for gen electrons",100,-3.4,3.4)};
	m_genElectronPhi[nCut][genType]->GetXaxis()-> SetTitle("Phi");

	m_genElectronEta[nCut][genType] = {m_histoFolder.make<TH1D>("genElectronEta","Eta for gen electrons",100,-3.4,3.4)};
	m_genElectronEta[nCut][genType]->GetXaxis()-> SetTitle("Eta");

	m_recoElectronPt[nCut][genType] = {m_histoFolder.make<TH1D>("recoElectronPt","Pt for reco electrons",100,0,500)};
	m_recoElectronPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_recoElectronPhi[nCut][genType] = {m_histoFolder.make<TH1D>("recoElectronPhi","Phi for reco electrons",100,-3.4,3.4)};
	m_recoElectronPhi[nCut][genType]->GetXaxis()-> SetTitle("Phi");

	m_recoElectronEta[nCut][genType] = {m_histoFolder.make<TH1D>("recoElectronEta","Eta for reco electrons",100,-3.4,3.4)};
	m_recoElectronEta[nCut][genType]->GetXaxis()-> SetTitle("Eta");


	m_WrongSignGenRecoElectronPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignGenRecoElectronPtRatio" , "ratio of gen pt to reco pt (wrong charge matched electrons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_WrongSignRecoElectronPt[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronPt" , "Pt for wrong sign reco electrons" , 100 , 0 , 500 )};
	m_WrongSignRecoElectronPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_WrongSignRecoElectronPhi[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronPhi" , "phi for wrong sign electrons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoElectronPhi[nCut][genType]->GetXaxis()->SetTitle("phi");

	m_WrongSignRecoElectronEta[nCut][genType] = {m_histoFolder.make<TH1D>("WrongSignRecoElectronEta" , "Eta for wrong sign electrons" , 100 , -3.4 , 3.4 )};
	m_WrongSignRecoElectronEta[nCut][genType]->GetXaxis()->SetTitle("eta");

	m_RightSignGenRecoElectronPtRatio[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignGenRecoElectronPtRatio" , "ratio of gen pt to reco pt (right charge matched electrons)" , 100 , 0 , 10 )};
	gStyle->SetOptStat("omen");

	m_RightSignRecoElectronPt[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoElectronPt" , "Pt for right sign reco electrons" , 100 , 0 , 500 )};
	m_RightSignRecoElectronPt[nCut][genType]->GetXaxis()-> SetTitle("Pt (GeV)");
	gStyle->SetOptStat("omen");

	m_RightSignRecoElectronPhi[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoElectronPhi" , "phi for right sign electrons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoElectronPhi[nCut][genType]->GetXaxis()->SetTitle("phi");

	m_RightSignRecoElectronEta[nCut][genType] = {m_histoFolder.make<TH1D>("RightSignRecoElectronEta" , "Eta for right sign electrons" , 100 , -3.4 , 3.4 )};
	m_RightSignRecoElectronEta[nCut][genType]->GetXaxis()->SetTitle("eta");
	
	m_bestElectrondR[nCut][genType] = {m_histoFolder.make<TH1D>("bestElectrondR" , "smallest deltaR (electrons)" , 50 , 0 , 0.03 )};
	m_bestElectrondR[nCut][genType]->GetXaxis()-> SetTitle("dR");
	gStyle->SetOptStat("omen");


	m_TWgenMuon[nCut][genType] = {m_histoFolder.make<TH1D>("TWgenMuon" , "t->W status for muons" , 2 , 1 , 2 )};
	m_TWgenMuon[nCut][genType]->GetXaxis()->SetBinLabel(1,"traced to tW");
	m_TWgenMuon[nCut][genType]->GetXaxis()->SetBinLabel(2,"not traced to tW");
	m_TWgenElectron[nCut][genType] = {m_histoFolder.make<TH1D>("TWgenElectron" , "t->W status for electrons" , 2 , 1 , 2 )};
	m_TWgenElectron[nCut][genType]->GetXaxis()->SetBinLabel(1,"traced to tW");
	m_TWgenElectron[nCut][genType]->GetXaxis()->SetBinLabel(2,"not traced to tW");

	m_phiEta[nCut][genType] = {m_histoFolder.make<TH2D>("phiEta" , "phi vs. eta for unmatched gen leptons" , 100 , -3.4 , 3.4, 100, -3.4 , 3.4 )};

	

}

//General histogram filling
void eventHistos2::fill(eventBits2& event, int cutNumber) {

int iSize;
int n;

if(event.twoMuons){
	for(j=0;j<2;j++){
	if(j==0){n=0;}
	if(j==1){n=3;}

		m_numRecoMuons[cutNumber][n]->Fill(event.muonRecoCount);
		m_numGenMuons[cutNumber][n]->Fill(event.muonGenCount);
		m_numRecoElectrons[cutNumber][n]->Fill(event.electronRecoCount);
		m_numGenElectrons[cutNumber][n]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[cutNumber][n]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[cutNumber][n]->Fill(event.electronGenCount-event.electronRecoCount);

	if(!event.failedMatch){
		m_matchStatus[cutNumber][n]->Fill("matched",1);


		m_genMuonPt[cutNumber][n]->Fill(event.lepton1Pt);
		m_genMuonPt[cutNumber][n]->Fill(event.lepton2Pt);

		m_genMuonPhi[cutNumber][n]->Fill(event.lepton1Phi);
		m_genMuonPhi[cutNumber][n]->Fill(event.lepton2Phi);

		m_genMuonEta[cutNumber][n]->Fill(event.lepton1Eta);
		m_genMuonEta[cutNumber][n]->Fill(event.lepton2Eta);

		m_recoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
		m_recoMuonPt[cutNumber][n]->Fill(event.Muon2Pt);

		m_recoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
		m_recoMuonPhi[cutNumber][n]->Fill(event.Muon2Phi);

		m_recoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		m_recoMuonEta[cutNumber][n]->Fill(event.Muon2Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);

			//std::cout<<"ptRatio: ";
			//std::cout<<event.Muon1PtRatio<<std::endl;
			//std::cout<<"event weight: ";
			//std::cout<<event.eventWeight<<std::endl;
			//std::cout<<"--------------------"<<std::endl;
		}

		if(event.Muon2chargeMatch){
			m_RightSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon2PtRatio); 
			m_RightSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon2Phi);
			m_RightSignRecoMuonPt[cutNumber][n]->Fill(event.Muon2Pt);
			m_RightSignRecoMuonEta[cutNumber][n]->Fill(event.Muon2Eta);
		}

		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		}

		if(!event.Muon2chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon2PtRatio); 
			m_WrongSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon2Phi);
			m_WrongSignRecoMuonPt[cutNumber][n]->Fill(event.Muon2Pt);
			m_WrongSignRecoMuonEta[cutNumber][n]->Fill(event.Muon2Eta);
		}

		if(event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("traced to tW",1);
		}
		if(event.Muon2TW){
			m_TWgenMuon[cutNumber][n]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("not traced to tW",1);
		}
		if(!event.Muon2TW){
			m_TWgenMuon[cutNumber][n]->Fill("not traced to tW",1);
		}


		m_bestMuondR[cutNumber][n]->Fill(event.Muon1dR);
		m_bestMuondR[cutNumber][n]->Fill(event.Muon2dR);
		m_genRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio);
		m_genRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon2PtRatio);

	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[cutNumber][n]->Fill("unmatched",1);

		if(!event.lepton1Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton1Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton1Phi,event.lepton1Eta);
		}

		if(!event.lepton2Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton2Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton2Phi,event.lepton2Eta);
		}
		
		

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPhi[cutNumber][n]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoEta[cutNumber][n]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPt[cutNumber][n]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoDR[cutNumber][n]->Fill(event.unmatchedDR[iValue]);
		}
	}
}
} 
/*	
if(event.twoElectrons){
	n=1

		m_numRecoMuons[cutNumber][n]->Fill(event.muonRecoCount);
		m_numGenMuons[cutNumber][n]->Fill(event.muonGenCount);
		m_numRecoElectrons[cutNumber][n]->Fill(event.electronRecoCount);
		m_numGenElectrons[cutNumber][n]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[cutNumber][n]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[cutNumber][n]->Fill(event.electronGenCount-event.electronRecoCount);

	if(!event.failedMatch){
		m_matchStatus[cutNumber][n]->Fill("matched",1);


		m_genElectronPt[cutNumber][n]->Fill(event.lepton1Pt);
		m_genElectronPt[cutNumber][n]->Fill(event.lepton2Pt);

		m_genElectronPhi[cutNumber][n]->Fill(event.lepton1Phi);
		m_genElectronPhi[cutNumber][n]->Fill(event.lepton2Phi);

		m_genElectronEta[cutNumber][n]->Fill(event.lepton1Eta);
		m_genElectronEta[cutNumber][n]->Fill(event.lepton2Eta);

		m_recoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
		m_recoElectronPt[cutNumber][n]->Fill(event.Electron2Pt);

		m_recoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
		m_recoElectronPhi[cutNumber][n]->Fill(event.Electron2Phi);

		m_recoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		m_recoElectronEta[cutNumber][n]->Fill(event.Electron2Eta);

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}

		if(event.Electron2chargeMatch){
			m_RightSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron2PtRatio); 
			m_RightSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron2Phi);
			m_RightSignRecoElectronPt[cutNumber][n]->Fill(event.Electron2Pt);
			m_RightSignRecoElectronEta[cutNumber][n]->Fill(event.Electron2Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}

		if(!event.Electron2chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron2PtRatio); 
			m_WrongSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron2Phi);
			m_WrongSignRecoElectronPt[cutNumber][n]->Fill(event.Electron2Pt);
			m_WrongSignRecoElectronEta[cutNumber][n]->Fill(event.Electron2Eta);
		}
		if(event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("traced to tW",1);
		}
		if(event.Electron2TW){
			m_TWgenElectron[cutNumber][n]->Fill("traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("not traced to tW",1);
		}
		if(!event.Electron2TW){
			m_TWgenElectron[cutNumber][n]->Fill("not traced to tW",1);
		}
		m_bestElectrondR[cutNumber][n]->Fill(event.Electron1dR);
		m_bestElectrondR[cutNumber][n]->Fill(event.Electron2dR);

		m_genRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio);
		m_genRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron2PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[cutNumber][n]->Fill("unmatched",1);

		if(!event.lepton1Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton1Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton1Phi,event.lepton1Eta);
		}

		if(!event.lepton2Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton2Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton2Phi,event.lepton2Eta);
		}

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPhi[cutNumber][n]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoEta[cutNumber][n]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPt[cutNumber][n]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoDR[cutNumber][n]->Fill(event.unmatchedDR[iValue]);
		}


	} 
}
*/

	
if(event.muonElectron){
for(j=0;j<2;j++){

	if(j==0){n=1;}
	if(j==1){n=3;}
		m_numRecoMuons[cutNumber][n]->Fill(event.muonRecoCount);
		m_numGenMuons[cutNumber][n]->Fill(event.muonGenCount);
		m_numRecoElectrons[cutNumber][n]->Fill(event.electronRecoCount);
		m_numGenElectrons[cutNumber][n]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[cutNumber][n]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[cutNumber][n]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[cutNumber][n]->Fill("matched",1);



		if(abs(event.lepton1Id)==13){
			m_genMuonPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_genMuonPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_genMuonEta[cutNumber][n]->Fill(event.lepton1Eta);

			m_genElectronPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_genElectronPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_genElectronEta[cutNumber][n]->Fill(event.lepton2Eta);
		}
		else{ 		
			m_genMuonPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_genMuonPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_genMuonEta[cutNumber][n]->Fill(event.lepton2Eta);

			m_genElectronPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_genElectronPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_genElectronEta[cutNumber][n]->Fill(event.lepton1Eta);
		}

		m_recoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
		m_recoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
		m_recoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);

		m_recoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
		m_recoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
		m_recoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		}

		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		}

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}

		if(event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("traced to tW",1);
		}
		if(event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("not traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("not traced to tW",1);
		}

		m_bestMuondR[cutNumber][n]->Fill(event.Muon1dR);
		m_bestElectrondR[cutNumber][n]->Fill(event.Electron1dR);
		m_genRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio);
		m_genRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[cutNumber][n]->Fill("unmatched",1);

		if(!event.lepton1Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton1Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton1Phi,event.lepton1Eta);
		}

		if(!event.lepton2Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton2Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton2Phi,event.lepton2Eta);
		}

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPhi[cutNumber][n]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoEta[cutNumber][n]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPt[cutNumber][n]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoDR[cutNumber][n]->Fill(event.unmatchedDR[iValue]);
		}

	}
}
}

if(event.muonTau){
	for(j=0;j<2;j++){

	if(j==0){n=2;}
	if(j==1){n=3;}

		m_numRecoMuons[cutNumber][n]->Fill(event.muonRecoCount);
		m_numGenMuons[cutNumber][n]->Fill(event.muonGenCount);
		m_numRecoElectrons[cutNumber][n]->Fill(event.electronRecoCount);
		m_numGenElectrons[cutNumber][n]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[cutNumber][n]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[cutNumber][n]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[cutNumber][n]->Fill("matched",1);


		if(abs(event.lepton1Id)==13){
			m_genMuonPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_genMuonPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_genMuonEta[cutNumber][n]->Fill(event.lepton1Eta);
		}
		else{ 		
			m_genMuonPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_genMuonPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_genMuonEta[cutNumber][n]->Fill(event.lepton2Eta);
		}

		m_recoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
		m_recoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
		m_recoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);

		if(event.Muon1chargeMatch){
			m_RightSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_RightSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_RightSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_RightSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		}


		if(!event.Muon1chargeMatch){
			m_WrongSignGenRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio); 
			m_WrongSignRecoMuonPhi[cutNumber][n]->Fill(event.Muon1Phi);
			m_WrongSignRecoMuonPt[cutNumber][n]->Fill(event.Muon1Pt);
			m_WrongSignRecoMuonEta[cutNumber][n]->Fill(event.Muon1Eta);
		}
		if(event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("traced to tW",1);
		}
		if(!event.Muon1TW){
			m_TWgenMuon[cutNumber][n]->Fill("not traced to tW",1);
		}

		m_bestMuondR[cutNumber][n]->Fill(event.Muon1dR);
		m_genRecoMuonPtRatio[cutNumber][n]->Fill(event.Muon1PtRatio);

	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[cutNumber][n]->Fill("unmatched",1);

		if(!event.lepton1Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton1Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton1Phi,event.lepton1Eta);
		}

		if(!event.lepton2Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton2Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton2Phi,event.lepton2Eta);
		}

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPhi[cutNumber][n]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoEta[cutNumber][n]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPt[cutNumber][n]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoDR[cutNumber][n]->Fill(event.unmatchedDR[iValue]);
		}

	}
}
}

/*
if(event.electronTau){
	n=4

		m_numRecoMuons[cutNumber][n]->Fill(event.muonRecoCount);
		m_numGenMuons[cutNumber][n]->Fill(event.muonGenCount);
		m_numRecoElectrons[cutNumber][n]->Fill(event.electronRecoCount);
		m_numGenElectrons[cutNumber][n]->Fill(event.electronGenCount);
		m_genMinusRecoMuons[cutNumber][n]->Fill(event.muonGenCount-event.muonRecoCount);
		m_genMinusRecoElectrons[cutNumber][n]->Fill(event.electronGenCount-event.electronRecoCount);
	if(!event.failedMatch){
		m_matchStatus[cutNumber][n]->Fill("matched",1);


		if(abs(event.lepton1Id)==11){
			m_genElectronPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_genElectronPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_genElectronEta[cutNumber][n]->Fill(event.lepton1Eta);
		}
		else{ 		
			m_genElectronPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_genElectronPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_genElectronEta[cutNumber][n]->Fill(event.lepton2Eta);
		}

		m_recoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
		m_recoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
		m_recoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);

		if(event.Electron1chargeMatch){
			m_RightSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_RightSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_RightSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_RightSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}

		if(!event.Electron1chargeMatch){
			m_WrongSignGenRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio); 
			m_WrongSignRecoElectronPhi[cutNumber][n]->Fill(event.Electron1Phi);
			m_WrongSignRecoElectronPt[cutNumber][n]->Fill(event.Electron1Pt);
			m_WrongSignRecoElectronEta[cutNumber][n]->Fill(event.Electron1Eta);
		}
				

		if(event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("traced to tW",1);
		}
		if(!event.Electron1TW){
			m_TWgenElectron[cutNumber][n]->Fill("not traced to tW",1);
		}

		m_bestElectrondR[cutNumber][n]->Fill(event.Electron1dR);
		m_genRecoElectronPtRatio[cutNumber][n]->Fill(event.Electron1PtRatio);
	}
	if(event.failedMatch && !event.failedGenPtEta){
		m_matchStatus[cutNumber][n]->Fill("unmatched",1);

		if(!event.lepton1Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton1Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton1Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton1Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton1Phi,event.lepton1Eta);
		}

		if(!event.lepton2Matched){
			m_unmatchedGenPt[cutNumber][n]->Fill(event.lepton2Pt);
			m_unmatchedGenEta[cutNumber][n]->Fill(event.lepton2Eta);
			m_unmatchedGenPhi[cutNumber][n]->Fill(event.lepton2Phi);
			m_phiEta[cutNumber][n]->Fill(event.lepton2Phi,event.lepton2Eta);
		}

		iSize=event.unmatchedPhi.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPhi[cutNumber][n]->Fill(event.unmatchedPhi[iValue]);
		}

		iSize=event.unmatchedEta.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoEta[cutNumber][n]->Fill(event.unmatchedEta[iValue]);
		}

		iSize=event.unmatchedPt.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoPt[cutNumber][n]->Fill(event.unmatchedPt[iValue]);
		}

		iSize=event.unmatchedDR.size();
		for(int iValue=0; iValue<iSize; iValue++){
			m_unmatchedRecoDR[cutNumber][n]->Fill(event.unmatchedDR[iValue]);
		}	
	}
}
*/

//top->cd();

}


	
