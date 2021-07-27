#ifndef cutFlowHistos_h
#define cutFlowHistos_h


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
#include "cutFlowHistos.h"
#include "eventBits2.h"


class cutFlowHistos {
	public:
	cutFlowHistos();
	//void book(TFileDirectory histoFolder);
	void book(TFileDirectory histoFolder, int nCut);
	void fill(double pT, int cutNumber);
	//cutFlowHistos() : m_histoFolder(5) {}



	private:
	
	TFileDirectory m_histoFolder;

	std::vector<TH1D*> m_recoMuonPt{10};

   
};

#endif