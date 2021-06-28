// system include files
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "ExoAnalysis/WR_lite/interface/eventBits.h"
#include "ExoAnalysis/WR_lite/interface/eventInfo.h"
#include "ExoAnalysis/WR_lite/interface/eventHistos.h"
#include "ExoAnalysis/WR_lite/interface/neuralNet.h"

