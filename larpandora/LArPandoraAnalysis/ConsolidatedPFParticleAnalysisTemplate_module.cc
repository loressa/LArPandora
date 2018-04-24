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

#include "TTree.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  ConsolidatedPFParticleAnalysisTemplate class
 */
class ConsolidatedPFParticleAnalysisTemplate : public art::EDAnalyzer
{
public:
    typedef art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap;
    typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
    typedef std::map< art::Ptr<recob::PFParticle>, PFParticleVector > PFParticlePrimaryMap;
    typedef std::vector< art::Ptr<recob::Track> > TrackVector;
    typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
  typedef std::map< art::Ptr<recob::PFParticle>, HitVector> PFParticleHitsMap;
  typedef std::map< art::Ptr<simb::MCParticle>, PFParticleHitsMap > MCParticlePFParticleHitVector; 

    /**
     *  @brief  Constructor
     *
     *  @param  pset the set of input fhicl parameters
     */
    ConsolidatedPFParticleAnalysisTemplate(fhicl::ParameterSet const &pset);
    
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

  unsigned int m_matchingMinSharedHits;
    float m_matchingMinCompleteness;
    float m_matchingMinPurity;

private:

    TTree       *m_pRecoTree;             ///<
    //LORENA - I need: number of true CR muons, number of neutrino-induced hits,
    // number of pfparticles + tracks matching true CR in pandoraNu, and consolidated
    // consolidated isCR matches
    // think about neutrino position, etc... - is there a vertex in the middle?
    
    // one entry per true CR muons - loop over MC particles, only go here isTrueCR
    int          m_run;                   ///<
    int          m_event;                 ///<
    int          m_index;                 ///<
    int          m_nuHits;                ///< number of neutrino-induced hits
    
    int          m_nRecoMatchesPandoraNu;
    int          m_nTrackMatchesPandoraNu;
    int          m_nTrackMatchesSameVertexPandoraNu; ///number of tracks matching from same vertex
    int          m_nRecoMatches;
    int          m_nTrackMatches;
    bool         m_isConsolidatedCR;
    int         m_isConsolidatedCRInt;
   
    
    
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> > MCParticleSet;
  typedef std::set< art::Ptr<simb::MCTruth> > MCTruthSet;
  typedef std::map< art::Ptr<simb::MCParticle>,  std::vector< art::Ptr<recob::PFParticle> > > MCParticlesToPFParticleVector;
  
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
    void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers);
    
    
    
    
    void GetTrueToRecoMatches(const MCParticlesToHits &trueParticlesToHits, const HitsToPFParticles &hitsToParticles,
                              MCParticlesToPFParticleVector &matchedParticles, MCParticlePFParticleHitVector &matchedHits) const;
    
    void GetTrueToRecoMatchesDifferentCollections(const MCParticlesToHits &trueParticlesToHits, const HitsToPFParticles &hitsToParticles,
                                                  MCParticlesToPFParticleVector &matchedParticles, MCParticlePFParticleHitVector &matchedHits) const;
    
    
 
  /**
   *  @brief Perform comparison between true reco matches and print out
   *
   */
  void CompareTrueRecoParticles(const MCParticlesToHits &trueParticlesToHits, const MCParticlesToMCTruth &particlesToTruth, const MCParticleMap &trueParticleMap, const MCParticlesToPFParticleVector &matchedParticlesPandoraNu, const MCParticlesToPFParticleVector &matchedParticles, const PFParticlesToVertices &recoParticlesToVertices, const PFParticlesToVertices &recoParticlesToVerticesPandoraNu, const PFParticleMap &recoParticleMap, const PFParticleMap &recoParticleMapPandoraNu, const PFParticlesToHits &particlesToHits, const PFParticlesToHits &particlesToHitsPandoraNu, const MCParticlePFParticleHitVector &matchedHits,const MCParticlePFParticleHitVector &matchedHitsPandoraNu);



    std::string m_pandoraLabel;         ///< The label for the pandora producer
    std::string m_trackLabel;           ///< The label for the track producer from PFParticles
    std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
    std::string m_pandoraNuLabel;         ///< The label for the pandora producer
    std::string m_trackPandoraNuLabel;           ///< The label for the track producer from PFParticles
    std::string m_showerPandoraNuLabel;          ///< The label for the shower producer from PFParticles
  
  std::string  m_hitfinderLabel;         ///<
  std::string  m_backtrackerLabel;       ///<
  std::string  m_geantModuleLabel;       ///<
  
  bool         m_useDaughterPFParticles; ///<
  bool         m_useDaughterMCParticles; ///<
  bool         m_addDaughterPFParticles; ///<
  bool         m_addDaughterMCParticles; ///<
  bool         m_recursiveMatching;      ///<    

};

DEFINE_ART_MODULE(ConsolidatedPFParticleAnalysisTemplate)

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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "Pandora/PdgTable.h"

#include <iostream>

namespace lar_pandora
{

ConsolidatedPFParticleAnalysisTemplate::ConsolidatedPFParticleAnalysisTemplate(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConsolidatedPFParticleAnalysisTemplate::reconfigure(fhicl::ParameterSet const &pset)
{

  m_matchingMinSharedHits = 5;
  m_matchingMinCompleteness = 0.1f;
  m_matchingMinPurity = 0.5f;

    m_pandoraLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_trackLabel = pset.get<std::string>("TrackModule","pandora");
    m_showerLabel = pset.get<std::string>("ShowerModule","pandora");
    m_pandoraNuLabel = pset.get<std::string>("PFParticleModule","pandoraNu");
    m_trackPandoraNuLabel = pset.get<std::string>("TrackModule","pandoraNu");
    m_showerPandoraNuLabel = pset.get<std::string>("ShowerModule","pandoraNu");

    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
    m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
    //    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","pandoraCosmicHitRemoval");
    // m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","crHitRemovalTruthMatch");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");

    m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles",false);
    m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles",true);
    m_addDaughterPFParticles = pset.get<bool>("AddDaughterPFParticles",true);
    m_addDaughterMCParticles = pset.get<bool>("AddDaughterMCParticles",true);
    m_recursiveMatching = pset.get<bool>("RecursiveMatching",false);
}
    
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::beginJob()
{
    mf::LogDebug("LArPandora") << " *** ConsolidatedPFParticleAnalysisTemplate::beginJob() *** " << std::endl;
    std::cout << " *** ConsolidatedPFParticleAnalysisTemplate::beginJob() *** " << std::endl;
    
    //
    art::ServiceHandle<art::TFileService> tfs;
    
    m_pRecoTree = tfs->make<TTree>("pandora", "LAr PFParticles");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("index", &m_index, "index/I");
    m_pRecoTree->Branch("nuHits", &m_nuHits, "nuHits/I");
    m_pRecoTree->Branch("nRecoMatchesPandoraNu", &m_nRecoMatchesPandoraNu, "nRecoMatchesPandoraNu/I");
    m_pRecoTree->Branch("nTrackMatchesPandoraNu", &m_nTrackMatchesPandoraNu, "nTrackMatchesPandoraNu/I");
    m_pRecoTree->Branch("nTrackMatchesSameVertexPandoraNu", &m_nTrackMatchesSameVertexPandoraNu, "nTrackMatchesSameVertexPandoraNu/I");
    m_pRecoTree->Branch("nRecoMatches", &m_nRecoMatches, "nRecoMatches/I");
    m_pRecoTree->Branch("nTrackMatches", &m_nTrackMatches, "nTrackMatches/I");
    m_pRecoTree->Branch("isConsolidatedCR", &m_isConsolidatedCRInt, "isConsolidatedCR/I");
    
    
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConsolidatedPFParticleAnalysisTemplate::analyze(const art::Event &evt)
{
    
    m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;
    
    // Collect the PFParticles from the event
    // ====================================================
    PFParticleHandle pfParticleHandle;
    PFParticleHandle pfParticlePandoraNuHandle;
    evt.getByLabel(m_pandoraLabel, pfParticleHandle);
    evt.getByLabel(m_pandoraNuLabel, pfParticlePandoraNuHandle);
   
    if (!pfParticleHandle.isValid())
    {
        mf::LogDebug("ConsolidatedPFParticleAnalysisTemplate") << "  Failed to find the PFParticles." << std::endl;
        return;
    }
    if (!pfParticlePandoraNuHandle.isValid())
    {
        mf::LogDebug("ConsolidatedPFParticleAnalysisTemplate") << "  Failed to find the PFParticles." << std::endl;
        return;
    }
    
    // Collect MCParticles and match True Particles to Hits
    // ====================================================
    MCParticleVector trueParticleVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;
    
    LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, trueParticleVector);
    LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
    
    // Collect Vertices and PFParticle <-> Vertex Associations
    // =======================================================
    VertexVector recoVertexVector, recoVertexVectorPandoraNu;
    PFParticlesToVertices recoParticlesToVertices, recoParticlesToVerticesPandoraNu;
    LArPandoraHelper::CollectVertices(evt, m_pandoraNuLabel, recoVertexVectorPandoraNu, recoParticlesToVerticesPandoraNu);
    LArPandoraHelper::CollectVertices(evt, m_pandoraLabel, recoVertexVector, recoParticlesToVertices);
    
    // Collect Hits
    // ============
    HitVector hitVector;
    LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);
    
    // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
    // =============================================================================
    PFParticleIdMap pfParticleMap, pfParticleMapPandoraNu;
    this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);
    this->GetPFParticleIdMap(pfParticlePandoraNuHandle, pfParticleMapPandoraNu);
 
    // Produce  PFParticle vectors containing final-state particles:
    // =============================================================
    std::vector< art::Ptr<recob::PFParticle> > crParticles;
    std::vector< art::Ptr<recob::PFParticle> > nuParticles;
    this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, nuParticles);

    PFParticlePrimaryMap finalParticles;
    this->GetFinalStatePFParticleVectors(pfParticleMapPandoraNu,finalParticles);

    // Vectors for tracks and showers for the final-states particles
    // =============================================================
    std::vector< art::Ptr<recob::Track> > tracks, trackscr, tracksPandoraNu;
    std::vector< art::Ptr<recob::Shower> > showers, showerscr, showersPandoraNu;
    this->CollectTracksAndShowers(nuParticles, pfParticleHandle, evt, tracks, showers);
    this->CollectTracksAndShowers(crParticles, pfParticleHandle, evt, trackscr, showerscr);
    this->CollectTracksAndShowers(nuParticles, pfParticlePandoraNuHandle, evt, tracksPandoraNu, showersPandoraNu);

    // maps of PFParticles to hits
    // ===========================
    PFParticlesToHits particlesToHits, particlesToHitsPandoraNu;
    HitsToPFParticles hitsToParticles, hitsToParticlesPandoraNu;
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, particlesToHits, hitsToParticles,
         (m_useDaughterPFParticles ? (m_addDaughterPFParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraNuLabel, particlesToHitsPandoraNu, hitsToParticlesPandoraNu,
         (m_useDaughterPFParticles ? (m_addDaughterPFParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));
 
    //map of MCParticles to hits
    // ==========================
    LArPandoraHelper::BuildMCParticleHitMaps(evt, m_geantModuleLabel, m_hitfinderLabel, m_backtrackerLabel,trueParticlesToHits, trueHitsToParticles,
         (m_useDaughterMCParticles ? (m_addDaughterMCParticles ? LArPandoraHelper::kAddDaughters : LArPandoraHelper::kUseDaughters) : LArPandoraHelper::kIgnoreDaughters));

    // Print a summary of the consolidated event
    std::cout << "Consolidated event summary:" << std::endl;
    std::cout << "  - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << std::endl;
    std::cout << "    ... of which are track-like   : " << trackscr.size() << std::endl;
    std::cout << "    ... of which are showers-like : " << showerscr.size() << std::endl;
    std::cout << "  - Number of neutrino final-state PFParticles : " << nuParticles.size() << std::endl;
    std::cout << "    ... of which are track-like   : " << tracks.size() << std::endl;
    std::cout << "    ... of which are showers-like : " << showers.size() << std::endl;
    
    std::cout << "PandoraNu event summary:" << std::endl;
    std::cout << "  - Found 'neutrinos': " << std::endl;
    for (PFParticlePrimaryMap::const_iterator it = finalParticles.begin(); it != finalParticles.end(); ++it)
    {
        const art::Ptr<recob::PFParticle> neutrino(it->first);
        const PFParticleVector daughters(it->second);
        std::cout << "    - PFParticle ID: " << neutrino->Self() << " is 'neutrino' with " << daughters.size() << " daughters" << std::endl;
    }

    std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;
    std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;
    std::cout << "  TrueParticles with hits: " << trueParticlesToHits.size() << std::endl;

    // Build Reco and True Particle Maps (for Parent/Daughter Navigation)                                                                                                  
    // =================================================================                                                                                                   
    MCParticleMap trueParticleMap;
    PFParticleMap recoParticleMap, recoParticleMapPandoraNu;
    PFParticleVector recoParticleVector, recoParticleVectorPandoraNu;

    LArPandoraHelper::CollectPFParticles(evt, m_pandoraLabel, recoParticleVector);
    LArPandoraHelper::CollectPFParticles(evt, m_pandoraNuLabel, recoParticleVectorPandoraNu);

    LArPandoraHelper::BuildMCParticleMap(trueParticleVector, trueParticleMap);
    LArPandoraHelper::BuildPFParticleMap(recoParticleVector, recoParticleMap);
    LArPandoraHelper::BuildPFParticleMap(recoParticleVectorPandoraNu, recoParticleMapPandoraNu);

    // Match Reco Particles to True Particles                                                                   
    // ======================================        
    MCParticlesToPFParticleVector matchedParticles, matchedParticlesPandoraNu;
    MCParticlePFParticleHitVector matchedHits, matchedHitsPandoraNu;
    
    this->GetTrueToRecoMatches(trueParticlesToHits, hitsToParticles, matchedParticles, matchedHits);
    this->GetTrueToRecoMatchesDifferentCollections(trueParticlesToHits, hitsToParticlesPandoraNu, matchedParticlesPandoraNu, matchedHitsPandoraNu);

    // Compare true and reconstructed particles
    // ========================================
    this->CompareTrueRecoParticles(trueParticlesToHits, particlesToTruth, trueParticleMap, matchedParticlesPandoraNu, matchedParticles, recoParticlesToVertices, recoParticlesToVerticesPandoraNu, recoParticleMap, recoParticleMapPandoraNu, particlesToHits, particlesToHitsPandoraNu, matchedHits, matchedHitsPandoraNu);

    
    //loop over MC particles, check whether true cosmic ray, and fill tree         m_pRecoTree->Fill();


}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConsolidatedPFParticleAnalysisTemplate::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
{
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
        if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
        {
            throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticleVector &crParticles, PFParticleVector &nuParticles)
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
            throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  This event contains multiple reconstructed neutrinos!";
        }

        // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
        for (const size_t daughterId : pParticle->Daughters())
        {
            if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  Invalid PFParticle collection!";

            nuParticles.push_back(pfParticleMap.at(daughterId));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, PFParticlePrimaryMap &finalParticles)
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
                throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  Invalid PFParticle collection!";

            daughters.push_back(pfParticleMap.at(daughterId));
        }
        finalParticles.insert(PFParticlePrimaryMap::value_type(pParticle, daughters));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers)
{
    // Get the associations between PFParticles and tracks/showers from the event
    art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
    art::FindManyP< recob::Shower > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);
   
    for (const art::Ptr<recob::PFParticle> &pParticle : particles)
    {
        const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
        const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
        const unsigned int nTracks(associatedTracks.size());
        const unsigned int nShowers(associatedShowers.size());

        // Check if the PFParticle has no associated tracks or showers
        if (nTracks == 0 && nShowers == 0)
        {
            mf::LogDebug("ConsolidatedPFParticleAnalysisTemplate") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
            continue;
        }

        // Check if there is an associated track
        if (nTracks == 1 && nShowers == 0)
        {
            tracks.push_back(associatedTracks.front());
            continue;
        }

        // Check if there is an associated shower
        if (nTracks == 0 && nShowers == 1)
        {
            showers.push_back(associatedShowers.front());
            continue;
        }

        throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------
void ConsolidatedPFParticleAnalysisTemplate::GetTrueToRecoMatches(const MCParticlesToHits &trueParticlesToHits, const HitsToPFParticles &hitsToParticles,
                                                                  MCParticlesToPFParticleVector &matchedParticles, MCParticlePFParticleHitVector &matchedHits) const
    
  {
      
    for (MCParticlesToHits::const_iterator iter1 = trueParticlesToHits.begin(), iterEnd1 = trueParticlesToHits.end();
         iter1 != iterEnd1; ++iter1)
      {
        const art::Ptr<simb::MCParticle> trueParticle = iter1->first;

        PFParticleSet recoMatches;
        PFParticleHitsMap recoParticleHits;
        const HitVector &hitVector = iter1->second;
        
        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {
            const art::Ptr<recob::Hit> hit = *iter2;
            
            HitsToPFParticles::const_iterator iter3 = hitsToParticles.find(hit);
            if (hitsToParticles.end() == iter3)
                continue;
        
            const art::Ptr<recob::PFParticle> recoParticle = iter3->second;

            recoParticleHits[recoParticle].push_back(hit);
            if (recoMatches.count(recoParticle) > 0)
              continue;
              
            recoMatches.insert(recoParticle);
            matchedParticles[trueParticle].push_back(recoParticle);
          }
        matchedHits[trueParticle] = recoParticleHits;
      }
  }


//------------------------------------------------------------------------------------------------------------------------------------------
  void ConsolidatedPFParticleAnalysisTemplate::GetTrueToRecoMatchesDifferentCollections(const MCParticlesToHits &trueParticlesToHits, const HitsToPFParticles &hitsToParticles, MCParticlesToPFParticleVector &matchedParticles, MCParticlePFParticleHitVector &matchedHits) const  
                                        
  {
      
      
    for (MCParticlesToHits::const_iterator iter1 = trueParticlesToHits.begin(), iterEnd1 = trueParticlesToHits.end();
           iter1 != iterEnd1; ++iter1)
      {
          const art::Ptr<simb::MCParticle> trueParticle = iter1->first;
          
          PFParticleSet recoMatches;
          PFParticleHitsMap recoParticleHits;
          const HitVector &hitVector = iter1->second;
          
          for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
          {
              const art::Ptr<recob::Hit> hit = *iter2;
              const double hitTime(hit->PeakTime());
              const geo::WireID wire(hit->WireID());
              
              for (HitsToPFParticles::const_iterator iter3 = hitsToParticles.begin(), iter3End = hitsToParticles.end(); iter3 != iter3End; ++iter3)
              {
                  const art::Ptr<recob::Hit> hit2 = iter3->first;
                  const double hitTime2(hit2->PeakTime());
                  const geo::WireID wire2(hit2->WireID());
                  
                  if ((hitTime == hitTime2) && (wire == wire2))
                  {
                        const art::Ptr<recob::PFParticle> recoParticle = iter3->second;
                        recoParticleHits[recoParticle].push_back(hit);

                        if (recoMatches.count(recoParticle) > 0)
                          continue;
                      
                        recoMatches.insert(recoParticle);
                        matchedParticles[trueParticle].push_back(recoParticle);
                  }
              }
          }
          matchedHits[trueParticle] = recoParticleHits;
      }
  }

//------------------------------------------------------------------------------------------------------------------------------------------
    
void ConsolidatedPFParticleAnalysisTemplate::CompareTrueRecoParticles(const MCParticlesToHits &trueParticlesToHits, const MCParticlesToMCTruth &particlesToTruth, const MCParticleMap &trueParticleMap, const MCParticlesToPFParticleVector &matchedParticlesPandoraNu, const MCParticlesToPFParticleVector &matchedParticles, const PFParticlesToVertices &recoParticlesToVertices, const PFParticlesToVertices &recoParticlesToVerticesPandoraNu, const PFParticleMap &recoParticleMap, const PFParticleMap &recoParticleMapPandoraNu, const PFParticlesToHits &particlesToHits, const PFParticlesToHits &particlesToHitsPandoraNu, const MCParticlePFParticleHitVector &matchedHits, const MCParticlePFParticleHitVector &matchedHitsPandoraNu)
{
 
    for (MCParticlesToHits::const_iterator iter = trueParticlesToHits.begin(), iterEnd = trueParticlesToHits.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> trueParticle = iter->first;
        const HitVector &trueHitVector = iter->second;
        
        if (trueHitVector.empty())
            continue;
        
        int m_mcPdg = trueParticle->PdgCode();
        int m_nMCHits = trueHitVector.size();
        std::cout << " ===============================================================" << std::endl;
        std::cout << " MC particle with PDG code " << m_mcPdg << " and " << m_nMCHits << " hits" << std::endl;
        // Get the true parent neutrino
        //check this
        MCParticlesToMCTruth::const_iterator nuIter = particlesToTruth.find(trueParticle);
        MCParticlePFParticleHitVector::const_iterator pandoraIter = matchedHits.find(trueParticle);
        MCParticlePFParticleHitVector::const_iterator pandoraNuIter = matchedHitsPandoraNu.find(trueParticle);
        std::map< art::Ptr<recob::PFParticle>, HitVector> recoParticleHits;
        const art::Ptr<simb::MCTruth> trueEvent = nuIter->second;
        
        if (trueEvent->NeutrinoSet())
        {
            const simb::MCNeutrino neutrino(trueEvent->GetNeutrino());
            int m_mcNuPdg = neutrino.Nu().PdgCode();
            bool m_mcIsCC = ((simb::kCC == neutrino.CCNC()) ? 1 : 0);
            std::cout <<  "   (corresponding neutrino: " << m_mcNuPdg  << " and is CC? " << m_mcIsCC << " ) " << std::endl;
        }
        try
        {
            const art::Ptr<simb::MCParticle> parentParticle(LArPandoraHelper::GetParentMCParticle(trueParticleMap, trueParticle));
            const art::Ptr<simb::MCParticle> primaryParticle(LArPandoraHelper::GetFinalStateMCParticle(trueParticleMap, trueParticle));
            int m_mcParentPdg = ((parentParticle != trueParticle) ? parentParticle->PdgCode() : 0);
            int m_mcPrimaryPdg = primaryParticle->PdgCode();
            bool m_mcIsPrimary = (primaryParticle == trueParticle);
            bool m_mcIsDecay = ("Decay" == trueParticle->Process());
            std::cout <<  "   whose parent is " << m_mcParentPdg << " and final state = " << m_mcPrimaryPdg << " is primary?" << m_mcIsPrimary << " is decay?" << m_mcIsDecay << std::endl;
            
        }
        catch (cet::exception &e){}
        
        //reco matches for pandoraNu
        std::cout << "   *** pandoraNu matching *** " << std::endl;
        MCParticlesToPFParticleVector::const_iterator pIter1 = matchedParticlesPandoraNu.find(trueParticle);
        if (matchedParticlesPandoraNu.end() != pIter1)
        {
            const PFParticleVector &particleVector = pIter1->second;
            if (!particleVector.empty())
            {
                for (const art::Ptr<recob::PFParticle> recoParticle : particleVector)
                {
              
                    PFParticlesToVertices::const_iterator pIter4 = recoParticlesToVerticesPandoraNu.find(recoParticle);
                    const PFParticleHitsMap &recoHitsMap = pandoraNuIter->second;
                    PFParticleHitsMap::const_iterator pandoraHitsIter = recoHitsMap.find(recoParticle);
                    const HitVector &hitVector = pandoraHitsIter->second;
                    //cuts on purity and completeness to consider a match
                    PFParticlesToHits::const_iterator recoHitsIter = particlesToHitsPandoraNu.find(recoParticle);
                    const HitVector &recoHits = recoHitsIter->second;
                    const float purity((recoHits.size() > 0) ? static_cast<float>(hitVector.size()) / static_cast<float>(recoHits.size()) : 0.f);
                    const float completeness((m_nMCHits > 0) ? static_cast<float>(hitVector.size()) / static_cast<float>(m_nMCHits) : 0.f);
                    if ((hitVector.size() < m_matchingMinSharedHits) || (purity < m_matchingMinPurity) || (completeness < m_matchingMinCompleteness))
                      continue;

                    std::cout << "     - reco match PFParticle " << recoParticle->Self() ;
                    std::cout << ((recoParticle->PdgCode()==13) ? ", track-like," : ", shower-like,");
                    std::cout << " with " << hitVector.size() << " matched hits";
                    
                    if (recoParticlesToVerticesPandoraNu.end() != pIter4)
                    {
                        const VertexVector &vertexVector = pIter4->second;
                        if (!vertexVector.empty())
                        {
                            if (vertexVector.size() !=1)
                                std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;
                    
                            const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
                            double xyz[3] = {0.0, 0.0, 0.0} ;
                            recoVertex->XYZ(xyz);
                            std::cout << " and vertex : (" << xyz[0] << " , " << xyz[1] << " , " << xyz[2] << ") " <<  std::endl;
                        }
                        else
                        {
                            std::cout << " " << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << " " << std::endl;
                    }
            
                    const art::Ptr<recob::PFParticle> parentParticle = LArPandoraHelper::GetParentPFParticle(recoParticleMapPandoraNu, recoParticle);
                    int m_pfoParentPdg = parentParticle->PdgCode();
                    if (parentParticle->Self() != recoParticle->Self())
                      {
                        std::cout << "        ... daughter of PFParticle " << parentParticle->Self() << " with pdg = " << m_pfoParentPdg ; //<< std::endl;
                        PFParticlesToVertices::const_iterator pIter5 = recoParticlesToVerticesPandoraNu.find(parentParticle);
                        if (recoParticlesToVerticesPandoraNu.end() != pIter5)
                          {
                            const VertexVector &vertexVector = pIter5->second;
                            if (!vertexVector.empty())
                              {
                                if (vertexVector.size() !=1)
                                  std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;
                                
                                const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
                                double xyz[3] = {0.0, 0.0, 0.0} ;
                                recoVertex->XYZ(xyz);
                                std::cout << " and vertex : ("<< xyz[0] << " , " << xyz[1] <<" , " << xyz[2] << ") " <<  std::endl;
                              }
                          }
                      }
                }
            }
            
        }
        // and consolidated reco matches
        std::cout << "   *** consolidated matching *** " << std::endl;
        MCParticlesToPFParticleVector::const_iterator pIter6 = matchedParticles.find(trueParticle);
        if (matchedParticles.end() != pIter6)
        {
            const PFParticleVector &particleVector = pIter6->second;
            if (!particleVector.empty())
            {
                for (const art::Ptr<recob::PFParticle> recoParticle : particleVector)
                {

                  const PFParticleHitsMap &recoHitsMap = pandoraIter->second;
                  PFParticleHitsMap::const_iterator pandoraHitsIter = recoHitsMap.find(recoParticle);
                  const HitVector &hitVector =pandoraHitsIter->second;

                  //cuts on purity and completeness to consider a match
                  PFParticlesToHits::const_iterator recoHitsIter = particlesToHits.find(recoParticle);
                  const HitVector &recoHits = recoHitsIter->second;
                  const float purity((recoHits.size() > 0) ? static_cast<float>(hitVector.size()) / static_cast<float>(recoHits.size()) : 0.f);
                  const float completeness((m_nMCHits > 0) ? static_cast<float>(hitVector.size()) / static_cast<float>(m_nMCHits) : 0.f);
                  if ((hitVector.size() < m_matchingMinSharedHits) || (purity < m_matchingMinPurity) || (completeness < m_matchingMinCompleteness))
                    continue;
                  
                  std::cout << "     - reco match PFParticle " << recoParticle->Self() ;
                  std::cout << ((recoParticle->PdgCode()==13) ? ", track-like," : ", shower-like,");
                  std::cout << " with " << hitVector.size() << " matched hits" ;

                    PFParticlesToVertices::const_iterator pIter7 = recoParticlesToVertices.find(recoParticle);
                    //                    std::cout << "     - reco match PFParticle " << recoParticle->Self() << " with pdg = " << recoParticle->PdgCode() ;
                    if (recoParticlesToVertices.end() != pIter7)
                    {
                        const VertexVector &vertexVector = pIter7->second;
                        if (!vertexVector.empty())
                        {
                            if (vertexVector.size() !=1)
                                std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;
                            
                            const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
                            double xyz[3] = {0.0, 0.0, 0.0} ;
                            recoVertex->XYZ(xyz);
                            std::cout << " and vertex : (" << xyz[0] << " , " << xyz[1] << " , " << xyz[2] << ") " <<  std::endl;
                        }
                        else
                        {
                            std::cout << " " << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << " " << std::endl;
                    }
                
            
                    const art::Ptr<recob::PFParticle> parentParticle = LArPandoraHelper::GetParentPFParticle(recoParticleMap, recoParticle);
                    int m_pfoParentPdg = parentParticle->PdgCode();
                    if (parentParticle->Self() != recoParticle->Self())
                      {
                        std::cout << "        ... daugther of PFParticle " << parentParticle->Self() << " with pdg = " << m_pfoParentPdg ; // << std::endl;
                        PFParticlesToVertices::const_iterator pIter8 = recoParticlesToVertices.find(parentParticle);
                        if (recoParticlesToVertices.end() != pIter8)
                          {
                            const VertexVector &vertexVector = pIter8->second;
                            if (!vertexVector.empty())
                              {
                                if (vertexVector.size() !=1)
                                  std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;
                                
                                const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
                                double xyz[3] = {0.0, 0.0, 0.0} ;
                                recoVertex->XYZ(xyz);
                                std::cout << " and vertex : ("<< xyz[0] << " , " << xyz[1] <<" , " << xyz[2] << ") " <<  std::endl;
                              }
                          }
                      }
            
            
                }
            }
    
        }
    }
}

} //namespace lar_pandora
