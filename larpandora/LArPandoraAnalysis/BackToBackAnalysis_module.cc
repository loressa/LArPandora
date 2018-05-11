/**
 *  @file   larpandora/LArPandoraAnalysis/ConsolidatedPFParticleAnalysisTemplate_module.cc
 *
 *  @brief  A template analysis module for using the Pandora consolidated output
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "TTree.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  ConsolidatedPFParticleAnalysisTemplate class
 */
class BackToBackAnalysis : public art::EDAnalyzer
{
public:

	typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
	typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
	typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
	typedef std::map< art::Ptr<recob::PFParticle>, PFParticleVector > PFParticlePrimaryMap;
	typedef std::vector< art::Ptr<recob::Track> > TrackVector;
	typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;

	typedef std::map< art::Ptr<recob::PFParticle>, art::Ptr<recob::Track> > PFParticleToTrackMap;
	typedef std::map< art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower> > PFParticleToShowerMap;
  
    typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
	typedef std::set< art::Ptr<simb::MCParticle> > MCParticleSet;
	typedef std::set< art::Ptr<simb::MCTruth> > MCTruthSet;
	typedef std::map< art::Ptr<simb::MCParticle>,  std::vector< art::Ptr<recob::PFParticle> > > MCParticlesToPFParticleVector;
	typedef std::map< art::Ptr<recob::PFParticle>,  std::vector< art::Ptr<simb::MCParticle> > > PFParticleToMCParticleVector;
  
    typedef std::map< art::Ptr<recob::PFParticle>, HitVector> PFParticleHitsMap;
    typedef std::map< art::Ptr<simb::MCParticle>, HitVector> MCParticleHitsMap;
	typedef std::map< art::Ptr<simb::MCParticle>, PFParticleHitsMap > MCParticlePFParticleHitVector;
	typedef std::map< art::Ptr<recob::PFParticle>, MCParticleHitsMap > PFParticleMCParticleHitVector;
  
  
    /**
     *  @brief  Constructor
     *
     *  @param  pset the set of input fhicl parameters
     */
    BackToBackAnalysis(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief  Configure memeber variables using FHiCL parameters
     *
     *  @param  pset the set of input fhicl parameters
     */
    void reconfigure(fhicl::ParameterSet const &pset);
	
	void beginJob();
    void endJob();

    /**
     *  @brief  Analyze an event!
     *
     *  @param  evt the art event to analyze
     */
    void analyze(const art::Event &evt);

private:


	bool IsNeutrinoInduced(const art::Ptr<simb::MCParticle> pMCParticle, const MCParticlesToMCTruth &artMCParticlesToMCTruth) const;

    /**
     *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
     *
     *  @param  pfParticleHandle the handle for the PFParticle collection
     *  @param  pfParticleMap the mapping from ID to PFParticle
     */
    void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

    /**
     *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
     *
     *  @param  pfParticleMap the mapping from ID to PFParticle
     *  @param  crParticles a vector to hold the top-level PFParticles reconstructed under the cosmic hypothesis
     *  @param  nuParticles a vector to hold the final-states of the reconstruced neutrino
     */
    void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles);

	/**
	*  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
	*
	*  @param  pfParticleMap the mapping from ID to PFParticle
	*  @param  crParticles a vector to hold the top-level PFParticles reconstructed under the cosmic hypothesis
	*  @param  nuParticles a vector to hold the final-states of the reconstruced neutrino
	*/
	void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticlePrimaryMap &finalParticles);

    /**
     *  @brief  Collect associated tracks and showers to particles in an input particle vector
     *
     *  @param  particles a vector holding PFParticles from which to find the associated tracks and showers
     *  @param  pfParticleHandle the handle for the PFParticle collection
     *  @param  evt the art event to analyze
     *  @param  tracks a vector to hold the associated tracks
     *  @param  showers a vector to hold the associated showers
     */
	void CollectTracksAndShowers(const PFParticleVector &particles,const PFParticleHandle &pfParticleHandle, const art::Event &evt, PFParticleToTrackMap &tracks, PFParticleToShowerMap &showers, std::string &trackLabel, std::string &showerLabel);


	void FindBackToBackTracks(const PFParticleToTrackMap &tracks, const PFParticlesToVertices &particlesToVertices, int &pairs, int &pairsNu, int &pairsReal, int &pairsRealNu,
		const PFParticleMap &recoParticleMap,PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits);

	bool DaughtersOfSameRecoNeutrino(const art::Ptr<recob::PFParticle> &pParticle1, const art::Ptr<recob::PFParticle> &pParticle2, const PFParticleMap &recoParticleMap);
			       
	void GetRecoToTrueMatches(const PFParticlesToHits &particlesToHits, const HitsToMCParticles &hitsToMCParticles,
		PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits) const;	

	void GetRecoToTrueMatchesDifferentCollections(const PFParticlesToHits &particlesToHits, const HitsToMCParticles &hitsToMCParticles,
		PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits) const;		

	bool SameMCParticleMatch(const art::Ptr<recob::PFParticle> &pParticle1, const art::Ptr<recob::PFParticle> &pParticle2, PFParticleToMCParticleVector &matchedParticles, 
		PFParticleMCParticleHitVector &matchedHits);
		   

	std::string m_pandoraLabel;         ///< The label for the pandora producer
	std::string m_trackLabel;           ///< The label for the track producer from PFParticles
	std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
	std::string m_pandoraNuLabel;         ///< The label for the pandora producer
	std::string m_trackPandoraNuLabel;           ///< The label for the track producer from PFParticles
	std::string m_showerPandoraNuLabel;          ///< The label for the shower producer from PFParticles
  
	std::string  m_hitfinderLabel;         ///<   
	std::string  m_hitfinderLabel2;         ///<   
	std::string  m_geantModuleLabel;       ///<
	std::string  m_backtrackerLabel;       ///<
  
  
    TTree       *m_pRecoTree;             ///<
   
    int          m_run;                   ///<
    int          m_event;                 ///<
    int          m_index;                 ///<
	int          m_nuHits;
	int    		 m_crRemovedHits;
	int 		 m_gaushits;
	int		 	 m_backToBackTracks;
	int 		 m_backToBackTracksNu;
	int 		 m_backToBackTracksPandoraNu;
	int 		 m_backToBackTracksNuPandoraNu;
	int 		 m_backToBackTracksRealBroken;
	int 		 m_backToBackTracksRealBrokenPandoraNu;
	int 		 m_backToBackTracksRealBrokenNu;
	int 		 m_backToBackTracksRealBrokenNuPandoraNu;
	float 		 m_maxVertexDistance;
	float		 m_maxAngleTolerance;
  
	bool         m_useDaughterPFParticles; ///<
	bool         m_useDaughterMCParticles; ///<
	bool         m_addDaughterPFParticles; ///<
	bool         m_addDaughterMCParticles; ///<
	bool         m_recursiveMatching;      ///< 

};

DEFINE_ART_MODULE(BackToBackAnalysis)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Pandora/PdgTable.h"

#include <iostream>

namespace lar_pandora
{

BackToBackAnalysis::BackToBackAnalysis(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BackToBackAnalysis::reconfigure(fhicl::ParameterSet const &pset)
{
  m_pandoraLabel = pset.get<std::string>("PFParticleModule","pandora");
  m_trackLabel = pset.get<std::string>("TrackModule","pandora");
  m_showerLabel = pset.get<std::string>("ShowerModule","pandora");
  m_pandoraNuLabel = pset.get<std::string>("PFParticleModule","pandoraNu");
  m_trackPandoraNuLabel = pset.get<std::string>("TrackModule","pandoraNu");
  m_showerPandoraNuLabel = pset.get<std::string>("ShowerModule","pandoraNu");
  
  m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
  m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
  m_hitfinderLabel2 = pset.get<std::string>("HitFinderModule2","pandoraCosmicHitRemoval");
  m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
  
  m_maxVertexDistance = pset.get<float>("MaxVertexDistance",10.f);
  m_maxAngleTolerance = pset.get<float>("MaxAngleTolerance",0.2f);
  
  m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles",false);
  m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles",true);
  m_addDaughterPFParticles = pset.get<bool>("AddDaughterPFParticles",true);
  m_addDaughterMCParticles = pset.get<bool>("AddDaughterMCParticles",true);
  m_recursiveMatching = pset.get<bool>("RecursiveMatching",false);

}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void BackToBackAnalysis::beginJob()
{
    mf::LogDebug("LArPandora") << " *** BackToBackAnalysis::beginJob() *** " << std::endl;
    std::cout << " *** BackToBackAnalysis::beginJob() *** " << std::endl;
    
    //
    art::ServiceHandle<art::TFileService> tfs;
    
    m_pRecoTree = tfs->make<TTree>("pandora", "Back to back Tracks analysis");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("nuHits", &m_nuHits, "nuHits/I");
    m_pRecoTree->Branch("crRemovedHits", &m_crRemovedHits, "crRemovedHits/I");
    m_pRecoTree->Branch("gaushits", &m_gaushits, "gaushits/I");
    //    m_pRecoTree->Branch("index", &m_index, "index/I");
    m_pRecoTree->Branch("nPairs", &m_backToBackTracks, "nPairs/I");
    m_pRecoTree->Branch("nPairsPandoraNu", &m_backToBackTracksPandoraNu, "nPairsPandoraNu/I");
    m_pRecoTree->Branch("nPairsNu", &m_backToBackTracksNu, "nPairsNu/I");
    m_pRecoTree->Branch("nPairsNuPandoraNu", &m_backToBackTracksNuPandoraNu, "nPairsNuPandoraNu/I");
	m_pRecoTree->Branch("nPairsRealBroken", &m_backToBackTracksRealBroken, "nPairsRealBroken/I");
	m_pRecoTree->Branch("nPairsRealBrokenNu", &m_backToBackTracksRealBrokenNu, "nPairsRealBrokenNu/I");
	m_pRecoTree->Branch("nPairsRealBrokenPandoraNu", &m_backToBackTracksRealBrokenPandoraNu, "nPairsRealBrokenPandoraNu/I");
	m_pRecoTree->Branch("nPairsRealBrokenNuPandoraNu", &m_backToBackTracksRealBrokenNuPandoraNu, "nPairsRealBrokenNuPandoraNu/I");
	
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
void BackToBackAnalysis::endJob()
{
}


//------------------------------------------------------------------------------------------------------------------------------------------

void BackToBackAnalysis::analyze(const art::Event &evt)
{
    // Collect the PFParticles from the event
	m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;
	
	PFParticleHandle pfParticleHandle;
	PFParticleHandle pfParticlePandoraNuHandle;
	evt.getByLabel(m_pandoraLabel, pfParticleHandle);
	evt.getByLabel(m_pandoraNuLabel, pfParticlePandoraNuHandle);
  
	if (!pfParticleHandle.isValid())
    {
        mf::LogDebug("BackToBackAnalysis") << "  Failed to find the PFParticles." << std::endl;
        return;
    }


	// Collect Hits                                                                                                                           
	// ============                                                                                                                           
	HitVector gausHits,crRemovedHits;
	LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, gausHits);
	LArPandoraHelper::CollectHits(evt, m_hitfinderLabel2, crRemovedHits);

	m_gaushits = gausHits.size();
	m_crRemovedHits = crRemovedHits.size();

    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
	PFParticleIdMap pfParticleMap, pfParticleMapPandoraNu;
	this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
	this->GetPFParticleIdMap(pfParticlePandoraNuHandle, pfParticleMapPandoraNu);
    
	PFParticleMap recoParticleMap, recoParticleMapPandoraNu;
	PFParticleVector recoParticleVector, recoParticleVectorPandoraNu;

	LArPandoraHelper::CollectPFParticles(evt, m_pandoraLabel, recoParticleVector);
	LArPandoraHelper::CollectPFParticles(evt, m_pandoraNuLabel, recoParticleVectorPandoraNu);

	LArPandoraHelper::BuildPFParticleMap(recoParticleVector, recoParticleMap);
	LArPandoraHelper::BuildPFParticleMap(recoParticleVectorPandoraNu, recoParticleMapPandoraNu);

    // Produce two PFParticle vectors containing final-state particles:
    // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
    // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
    std::vector< art::Ptr<recob::PFParticle> > crParticles;
    std::vector< art::Ptr<recob::PFParticle> > nuParticles;
    this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, nuParticles);
    
    // Use as required!
    // -----------------------------
    //   What follows is an example showing how one might access the reconstructed neutrino final-state tracks and showers
    
    // These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
    PFParticleToTrackMap tracks, tracksPandoraNu;
    PFParticleToShowerMap showers, showersPandoraNu;
    this->CollectTracksAndShowers(recoParticleVector, pfParticleHandle, evt, tracks, showers, m_trackLabel, m_showerLabel);
    this->CollectTracksAndShowers(recoParticleVectorPandoraNu, pfParticlePandoraNuHandle, evt, tracksPandoraNu, showersPandoraNu, m_trackPandoraNuLabel, m_showerPandoraNuLabel);

    // Collect Vertices and PFParticle <-> Vertex Associations                                        
    // =======================================================                                                                          
    VertexVector recoVertexVector, recoVertexVectorPandoraNu;
    PFParticlesToVertices recoParticlesToVertices, recoParticlesToVerticesPandoraNu;
    LArPandoraHelper::CollectVertices(evt, m_pandoraNuLabel, recoVertexVectorPandoraNu, recoParticlesToVerticesPandoraNu);
    LArPandoraHelper::CollectVertices(evt, m_pandoraLabel, recoVertexVector, recoParticlesToVertices);
	
	 // Collect MCParticles and match True Particles to Hits
    // ====================================================
    MCParticleVector trueParticleVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;
    
    LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, trueParticleVector);
    LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
	
	//map of MCParticles to hits
    // ==========================
    LArPandoraHelper::BuildMCParticleHitMaps(evt, m_geantModuleLabel, m_hitfinderLabel, m_backtrackerLabel,trueParticlesToHits, trueHitsToParticles,
         (m_useDaughterMCParticles ? (m_addDaughterMCParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));
	
	//find neutrino hits
	m_nuHits = 0;
	for (MCParticlesToHits::const_iterator iter1 = trueParticlesToHits.begin(), iterEnd1 = trueParticlesToHits.end();
         iter1 != iterEnd1; ++iter1)
      {
        const art::Ptr<simb::MCParticle> trueParticle = iter1->first;
		if (!this->IsNeutrinoInduced(trueParticle, particlesToTruth))
			continue;
        const HitVector &hitVector = iter1->second;
		m_nuHits += hitVector.size();
	  }
	std::cout << " Neutrino hits = " << m_nuHits << std::endl;
		 
	// maps of PFParticles to hits
    // ===========================
    PFParticlesToHits particlesToHits, particlesToHitsPandoraNu;
    HitsToPFParticles hitsToParticles, hitsToParticlesPandoraNu;
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, particlesToHits, hitsToParticles,
         (m_useDaughterPFParticles ? (m_addDaughterPFParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraNuLabel, particlesToHitsPandoraNu, hitsToParticlesPandoraNu,
         (m_useDaughterPFParticles ? (m_addDaughterPFParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));

	
	 // Match Reco Particles to True Particles                                                                   
    // ======================================        
    PFParticleToMCParticleVector matchedParticles, matchedParticlesPandoraNu;
    PFParticleMCParticleHitVector matchedHits, matchedHitsPandoraNu;
    
    this->GetRecoToTrueMatches(particlesToHits, trueHitsToParticles, matchedParticles, matchedHits);
    this->GetRecoToTrueMatchesDifferentCollections(particlesToHitsPandoraNu, trueHitsToParticles, matchedParticlesPandoraNu, matchedHitsPandoraNu);
    
    std::cout << "pandora:" << std::endl;
    this->FindBackToBackTracks(tracks, recoParticlesToVertices, m_backToBackTracks, m_backToBackTracksNu, m_backToBackTracksRealBroken, m_backToBackTracksRealBrokenNu,  
								recoParticleMap, matchedParticles, matchedHits);
    std::cout << "pandoraNu:" << std::endl;
    this->FindBackToBackTracks(tracksPandoraNu, recoParticlesToVerticesPandoraNu, m_backToBackTracksPandoraNu, m_backToBackTracksNuPandoraNu, m_backToBackTracksRealBrokenPandoraNu, 
	                              m_backToBackTracksRealBrokenNuPandoraNu, recoParticleMapPandoraNu, matchedParticlesPandoraNu, matchedHitsPandoraNu);
   
    //  FillTree
	m_pRecoTree->Fill();

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BackToBackAnalysis::IsNeutrinoInduced(const art::Ptr<simb::MCParticle> pMCParticle, const MCParticlesToMCTruth &artMCParticlesToMCTruth) const
{
    MCParticlesToMCTruth::const_iterator iter = artMCParticlesToMCTruth.find(pMCParticle);

    if (artMCParticlesToMCTruth.end() == iter)
        return false;

    const art::Ptr<simb::MCTruth> pMCTruth = iter->second;
    return (simb::kBeamNeutrino == pMCTruth->Origin());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BackToBackAnalysis::GetRecoToTrueMatches(const PFParticlesToHits &particlesToHits, const HitsToMCParticles &hitsToMCParticles,
											  PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits) const
    
  {
      
    for (PFParticlesToHits::const_iterator iter1 = particlesToHits.begin(), iterEnd1 = particlesToHits.end();
         iter1 != iterEnd1; ++iter1)
      {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;

        MCParticleSet trueMatches;
        MCParticleHitsMap trueParticleHits;
        const HitVector &hitVector = iter1->second;
        
        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {
            const art::Ptr<recob::Hit> hit = *iter2;
            
            HitsToMCParticles::const_iterator iter3 = hitsToMCParticles.find(hit);
            if (hitsToMCParticles.end() == iter3)
                continue;
        
            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

            trueParticleHits[trueParticle].push_back(hit);
            if (trueMatches.count(trueParticle) > 0)
              continue;
              
            trueMatches.insert(trueParticle);
            matchedParticles[recoParticle].push_back(trueParticle);
          }
        matchedHits[recoParticle] = trueParticleHits;
      }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

void BackToBackAnalysis::GetRecoToTrueMatchesDifferentCollections(const PFParticlesToHits &particlesToHits, const HitsToMCParticles &hitsToMCParticles,
											                      PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits) const
    
  {
      
    for (PFParticlesToHits::const_iterator iter1 = particlesToHits.begin(), iterEnd1 = particlesToHits.end();
         iter1 != iterEnd1; ++iter1)
      {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;

        MCParticleSet trueMatches;
        MCParticleHitsMap trueParticleHits;
        const HitVector &hitVector = iter1->second;
        
        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {            
			  //ATTN: known issue with crHitRemovalTruthMatch, can't be used as in GetRecoToTrueMatches, instead, comapre wire ids and times
			  const art::Ptr<recob::Hit> hit = *iter2;
              const double hitTime(hit->PeakTime());
              const geo::WireID wire(hit->WireID());
              
              for (HitsToMCParticles::const_iterator iter3 = hitsToMCParticles.begin(), iter3End = hitsToMCParticles.end(); iter3 != iter3End; ++iter3)
              {
                  const art::Ptr<recob::Hit> hit2 = iter3->first;
                  const double hitTime2(hit2->PeakTime());
                  const geo::WireID wire2(hit2->WireID());
                  
                  if ((hitTime == hitTime2) && (wire == wire2))
				  {
        
						const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

						trueParticleHits[trueParticle].push_back(hit);
						if (trueMatches.count(trueParticle) > 0)
							continue;
              
						trueMatches.insert(trueParticle);
						matchedParticles[recoParticle].push_back(trueParticle);
				  }
			}
		}
        matchedHits[recoParticle] = trueParticleHits;
      }
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void BackToBackAnalysis::FindBackToBackTracks(const PFParticleToTrackMap &tracks, const PFParticlesToVertices &particlesToVertices,
										        int &pairs, int &pairsNu, int &pairsReal, int &pairsRealNu,const PFParticleMap &recoParticleMap, 
												PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits)
  {
    
	pairs = 0; //pairs of back to back tracks
    pairsNu = 0; //paris of back to back tracks coming from a reco neutrino vertex
	pairsReal = 0; //this are pairs of back to back related to real broken tracks (both pieces matching the same MC particle)  
	pairsRealNu = 0; //same as above with the addition of coming from a reco neutrino vertex
	
	PFParticleSet usedPFParticles;
	for (PFParticleToTrackMap::const_iterator it = tracks.begin(); it != tracks.end(); ++it)
    {
		
		const art::Ptr<recob::PFParticle> pParticle(it->first);
		const art::Ptr<recob::Track> pTrack(it->second);
		
		//check we haven't used it yet
		if (usedPFParticles.count(pParticle) > 0)
			continue;
		usedPFParticles.insert(pParticle);
		
		//Track->PFParticle->Vertex navigation
		PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(pParticle);
		if (particlesToVertices.end() == vIter)
			continue;
		
		const VertexVector &vertexVector = vIter->second;
		if (vertexVector.empty())
			continue;
			
		if (vertexVector.size() !=1)
			continue;
				
		const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
		double v1[3] = {0.0, 0.0, 0.0} ;
		recoVertex->XYZ(v1);
		
		//Track->direction (at vertex) navigation		
		TVector3 direction1 = pTrack->VertexDirection(); 
		
		//now loop over the rest of tracks 
		for (PFParticleToTrackMap::const_iterator it2 = tracks.begin(); it2 != tracks.end(); ++it2)
		{
		
			const art::Ptr<recob::PFParticle> pParticle2(it2->first);
			const art::Ptr<recob::Track> pTrack2(it2->second);
		
			//check we haven't used it yet - this will now skip pParticle itself
			if (usedPFParticles.count(pParticle2) > 0)
				continue;
		
			//Track->PFParticle->Vertex navigation
			PFParticlesToVertices::const_iterator vIter2 = particlesToVertices.find(pParticle2);
			if (particlesToVertices.end() == vIter2)
				continue;
		
			const VertexVector &vertexVector2 = vIter2->second;
			if (vertexVector2.empty())
				continue;
			
			if (vertexVector2.size() !=1)
				continue;
				
			const art::Ptr<recob::Vertex> recoVertex2 = *(vertexVector2.begin());
			double v2[3] = {0.0, 0.0, 0.0} ;
			recoVertex2->XYZ(v2);
		
			//Track->direction (at vertex) navigation		
			TVector3 direction2 = pTrack2->VertexDirection(); 
			
			//check proximity of vertices, using pandora::CartesianVector for convenience 
			pandora::CartesianVector vertex1(v1[0], v1[1], v1[2]);
			pandora::CartesianVector vertex2(v2[0], v2[1], v2[2]);
			if ((vertex1-vertex2).GetMagnitude() > m_maxVertexDistance)
				continue;
		
			//check angle approximately 180 degrees (+/-pi)
			double angle(direction1.Angle(direction2));	
			if ((std::abs(angle-3.14159) < m_maxAngleTolerance) || (std::abs(angle+3.14159) < m_maxAngleTolerance))
			{
				std::cout << " vertex1 = (" << vertex1.GetX() <<  ", "<< vertex1.GetY() << " ," <<  vertex1.GetZ() 
						  << "), vertex2 = (" << vertex2.GetX() <<  ", "<< vertex2.GetY() << " ," <<  vertex2.GetZ() << "), and angle " << angle << std::endl;
				std::cout << "pParticle = " << pParticle->Self() << " and pParticle2  = " << pParticle2->Self() << std::endl;

				++pairs;
				
				bool realBroken(this->SameMCParticleMatch(pParticle, pParticle2, matchedParticles, matchedHits));
				if (realBroken)
					++pairsReal;
					
				if (this->DaughtersOfSameRecoNeutrino(pParticle, pParticle2,recoParticleMap))
					{
						++pairsNu;
						if (realBroken)
							++pairsRealNu;
                                          
					}
					
				usedPFParticles.insert(pParticle2);
			}
		}
	}

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BackToBackAnalysis::SameMCParticleMatch(const art::Ptr<recob::PFParticle> &pParticle1, const art::Ptr<recob::PFParticle> &pParticle2, 
											PFParticleToMCParticleVector &matchedParticles, PFParticleMCParticleHitVector &matchedHits)
{

	PFParticleMCParticleHitVector::const_iterator iter1 = matchedHits.find(pParticle1);
	PFParticleMCParticleHitVector::const_iterator iter2 = matchedHits.find(pParticle2);
	
	PFParticleToMCParticleVector::const_iterator pIter1 = matchedParticles.find(pParticle1);
	PFParticleToMCParticleVector::const_iterator pIter2 = matchedParticles.find(pParticle2);
	
	if ((matchedParticles.end() == pIter1) || (matchedParticles.end() == pIter2))
		return false;
	
	const MCParticleVector &MCParticleVector1 = pIter1->second;	
	const MCParticleVector &MCParticleVector2 = pIter2->second;	
	
	if ((MCParticleVector1.empty()) || (MCParticleVector2.empty()))
		return false;
		
	//find the main MC particle matching the reco particle (i.e. that with most hits)	
	unsigned int maxHits1(0), maxHits2(0);	
	art::Ptr<simb::MCParticle> bestTrueParticleMatch1, bestTrueParticleMatch2;
	for (const art::Ptr<simb::MCParticle> trueParticle : MCParticleVector1)	
	{
		const MCParticleHitsMap &trueHitsMap = iter1->second;
		MCParticleHitsMap::const_iterator trueHitsIter = trueHitsMap.find(trueParticle);
		const HitVector &hitVector = trueHitsIter->second;
		if (hitVector.size() > maxHits1)
		{
			maxHits1 = hitVector.size();
			bestTrueParticleMatch1 = trueParticle;
		}
	}
	for (const art::Ptr<simb::MCParticle> trueParticle : MCParticleVector2)	
	{
		const MCParticleHitsMap &trueHitsMap = iter2->second;
		MCParticleHitsMap::const_iterator trueHitsIter = trueHitsMap.find(trueParticle);
		const HitVector &hitVector = trueHitsIter->second;
		if (hitVector.size() > maxHits2)
		{
			maxHits2 = hitVector.size();
			bestTrueParticleMatch2 = trueParticle;
		}
	}
 
	if (bestTrueParticleMatch1 == bestTrueParticleMatch2)
	{
		//std::cout << "bestTrueParticleMatch1 = " << bestTrueParticleMatch1->TrackId() << " and bestTrueParticleMatch2 = " << bestTrueParticleMatch2->TrackId() << std::endl;
		//std::cout << " MC PDG = " << bestTrueParticleMatch1->PdgCode() << std::endl;
		std::cout << " maxHits1 = " << maxHits1 << " and maxHits2 =" << maxHits2 << std::endl;
		return true;
	}
	//	return (bestTrueParticleMatch1 == bestTrueParticleMatch2);
	else 
		return false;
 
}

//------------------------------------------------------------------------------------------------------------------------------------------

  bool BackToBackAnalysis::DaughtersOfSameRecoNeutrino(const art::Ptr<recob::PFParticle> &pParticle1, const art::Ptr<recob::PFParticle> &pParticle2, 
														const PFParticleMap &recoParticleMap)
{

  const art::Ptr<recob::PFParticle> parentParticle1 = LArPandoraHelper::GetParentPFParticle(recoParticleMap, pParticle1);
  const art::Ptr<recob::PFParticle> parentParticle2 = LArPandoraHelper::GetParentPFParticle(recoParticleMap, pParticle2);
  
  if (parentParticle1->Self() == pParticle1->Self())
    return false;

  if (parentParticle2->Self() == pParticle2->Self())
    return false;
					
  if (parentParticle1->Self() != parentParticle2->Self())
    return false;
  
  if ((parentParticle1->PdgCode() != 14) && (parentParticle1->PdgCode() != 12))
    return false;

 // std::cout << " parentParticle1->PdgCode() = " << parentParticle1->PdgCode() << std::endl;

  return true;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BackToBackAnalysis::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
{
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
        if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
        {
            throw cet::exception("BackToBackAnalysis") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void BackToBackAnalysis::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles)
{
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
    {
        const art::Ptr<recob::PFParticle> pParticle(it->second);

        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;

        // Check if this particle is identified as the neutrino
        const int pdg(pParticle->PdgCode());
        const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

        // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
        if (!isNeutrino)
        {
            crParticles.push_back(pParticle);
            continue;
        }

        // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
        //       If this is not the case please handle accordingly
        if (!nuParticles.empty())
        {
            throw cet::exception("BackToBackAnalysis") << "  This event contains multiple reconstructed neutrinos!";
        }

        // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
        for (const size_t daughterId : pParticle->Daughters())
        {
            if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                throw cet::exception("BackToBackAnalysis") << "  Invalid PFParticle collection!";

            nuParticles.push_back(pfParticleMap.at(daughterId));
        }
    }
}

  //------------------------------------------------------------------------------------------------------------------------------------------

  void BackToBackAnalysis::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticlePrimaryMap &finalParticles)
                                                                              
  {
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
      {
        const art::Ptr<recob::PFParticle> pParticle(it->second);

        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;

        PFParticleVector daughters;
        // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
        for (const size_t daughterId : pParticle->Daughters())
          {
            if (pfParticleMap.find(daughterId) == pfParticleMap.end())
              throw cet::exception("BackToBackAnalysis") << "  Invalid PFParticle collection!";

            daughters.push_back(pfParticleMap.at(daughterId));
          }
        finalParticles.insert(PFParticlePrimaryMap::value_type(pParticle, daughters));
      }
  }


//------------------------------------------------------------------------------------------------------------------------------------------
    
  void BackToBackAnalysis::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, PFParticleToTrackMap &tracks, 
  PFParticleToShowerMap &showers, std::string &trackLabel, std::string &showerLabel)
{
    // Get the associations between PFParticles and tracks/showers from the event
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, trackLabel);
    art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, evt, showerLabel);
   
    for (const art::Ptr<recob::PFParticle> &pParticle : particles)
    {
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
        const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
        const unsigned int nTracks(associatedTracks.size());
        const unsigned int nShowers(associatedShowers.size());

        // Check if the PFParticle has no associated tracks or showers
        if (nTracks == 0 && nShowers == 0)
        {
            mf::LogDebug("BackToBackAnalysis") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
            continue;
        }

        // Check if there is an associated track
        if (nTracks == 1 && nShowers == 0)
        {
          const art::Ptr<recob::Track> track(associatedTracks.front());
          //            tracks.push_back(associatedTracks.front());
          tracks[pParticle] = track;
            continue;
        }

        // Check if there is an associated shower
        if (nTracks == 0 && nShowers == 1)
        {
          const art::Ptr<recob::Shower> shower(associatedShowers.front());
          //            showers.push_back(associatedShowers.front());
          showers[pParticle] = shower;
            continue;
        }

        throw cet::exception("BackToBackAnalysis") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
    }
}

} //namespace lar_pandora
