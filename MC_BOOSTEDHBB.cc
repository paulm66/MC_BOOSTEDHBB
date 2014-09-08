// -*- C++ -*-
#ifndef RIVET_MC_BOOSTEDHBB_HH
#define RIVET_MC_BOOSTEDHBB_HH

#include <iostream>
#include <map>
#include <string>

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"

using std::map;
using std::string;

namespace Rivet {


    class MC_BOOSTEDHBB : public Analysis {
        public:

            /// Constructor
            MC_BOOSTEDHBB()
                : Analysis("MC_BOOSTEDHBB") {

                    return;
                }


        public:

            /// @name Analysis methods
            //@{

            /// Book histograms and initialise projections before the run
            void init() {

                ChargedLeptons clfs(FinalState(-2.5, 2.5, 25*GeV));
                addProjection(clfs, "ChargedLeptons");

                MissingMomentum mmfs(FinalState(-4.2, 4.2, 0*GeV));
                addProjection(mmfs, "MissingMomentum");

                // calo jets
                // don't include high-pt neutrinos or leptons in jets
                // TODO
                // include electrons?
                IdentifiedFinalState nufs(FinalState(-2.5, 2.5, 25*GeV));
                nufs.acceptNeutrinos();

                MergedFinalState leptonsAndNeutrinos(clfs, nufs);
                addProjection(leptonsAndNeutrinos, "LeptonsAndNeutrinos");

                VetoedFinalState caloParts(FinalState(-4.2, 4.2));
                caloParts.addVetoOnThisFinalState(leptonsAndNeutrinos);

                addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.2), "AntiKt02CaloJets");
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), "AntiKt04CaloJets");
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.6), "AntiKt06CaloJets");
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.8), "AntiKt08CaloJets");
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 1.0), "AntiKt10CaloJets");


                // "track" jets
                ChargedFinalState trackParts(-2.5, 2.5, 0.5*GeV);

                addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.2), "AntiKt02TrackJets");
                addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.4), "AntiKt04TrackJets");
                addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.6), "AntiKt06TrackJets");
                addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.8), "AntiKt08TrackJets");
                addProjection(FastJets(trackParts, FastJets::ANTIKT, 1.0), "AntiKt10TrackJets");


                addJetCollection("AntiKt02CaloJets");
                addJetCollection("AntiKt04CaloJets");
                addJetCollection("AntiKt06CaloJets");
                addJetCollection("AntiKt08CaloJets");
                addJetCollection("AntiKt10CaloJets");

                addJetCollection("AntiKt02TrackJets");
                addJetCollection("AntiKt04TrackJets");
                addJetCollection("AntiKt06TrackJets");
                addJetCollection("AntiKt08TrackJets");
                addJetCollection("AntiKt10TrackJets");

                // TODO

                minJetPtCut["AntiKt02CaloJets"] = 25*GeV;
                minJetPtCut["AntiKt04CaloJets"] = 25*GeV;
                minJetPtCut["AntiKt06CaloJets"] = 25*GeV;
                minJetPtCut["AntiKt08CaloJets"] = 25*GeV;
                minJetPtCut["AntiKt10CaloJets"] = 25*GeV;

                minJetPtCut["AntiKt02TrackJets"] = 25*GeV;
                minJetPtCut["AntiKt04TrackJets"] = 25*GeV;
                minJetPtCut["AntiKt06TrackJets"] = 25*GeV;
                minJetPtCut["AntiKt08TrackJets"] = 25*GeV;
                minJetPtCut["AntiKt10TrackJets"] = 25*GeV;

                
                return;
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                // leptons
                // const Particles &leptons =
                    // applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

                fillJetCollection("AntiKt02CaloJets", event);
                fillJetCollection("AntiKt04CaloJets", event);
                fillJetCollection("AntiKt06CaloJets", event);
                fillJetCollection("AntiKt08CaloJets", event);
                fillJetCollection("AntiKt10CaloJets", event);

                fillJetCollection("AntiKt02TrackJets", event);
                fillJetCollection("AntiKt04TrackJets", event);
                fillJetCollection("AntiKt06TrackJets", event);
                fillJetCollection("AntiKt08TrackJets", event);
                fillJetCollection("AntiKt10TrackJets", event);

            }


            /// Normalise histograms etc., after the run
            void finalize() {

                for (map<string, map<string, Histo1DPtr> >::iterator p = jet1DHistos.begin(); p != jet1DHistos.end(); ++p) {
                    for (map<string, Histo1DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                        scale(q->second, crossSection()/sumOfWeights()); // norm to cross section
                    }
                }

                return;
            }

            //@}


        private:

            map<string, map<string, Histo1DPtr> > jet1DHistos;
            map<string, map<string, Profile1DPtr> > jet1DProfiles;
            map<string, double> minJetPtCut;

            void addJetCollection(const string &name) {
                jet1DHistos[name] = map<string, Histo1DPtr>();
                jet1DProfiles[name] = map<string, Profile1DPtr>();

                MSG_DEBUG("Adding jet collection " << name);

                map<string, Histo1DPtr> &mh1D = jet1DHistos.at(name);
                map<string, Profile1DPtr> &mp1D = jet1DProfiles.at(name);

                mh1D["_n"] = bookHisto1D(name + "_n", 10, 0, 10);
                mh1D["_pt"] = bookHisto1D(name + "_pt", 50, 0, 2000*GeV);
                mh1D["_eta"] = bookHisto1D(name + "_eta", 50, -5, 5);
                mh1D["_phi"] = bookHisto1D(name + "_phi", 50, 0, 2*PI);
                mh1D["_E"] = bookHisto1D(name + "_E", 50, 0, 2000*GeV);
                mh1D["_m"] = bookHisto1D(name + "_m", 50, 0, 500*GeV);

                mh1D["0_pt"] = bookHisto1D(name + "0_pt", 50, 0, 2000*GeV);
                mh1D["0_eta"] = bookHisto1D(name + "0_eta", 50, -5, 5);
                mh1D["0_phi"] = bookHisto1D(name + "0_phi", 50, 0, 2*PI);
                mh1D["0_E"] = bookHisto1D(name + "0_E", 50, 0, 2000*GeV);
                mh1D["0_m"] = bookHisto1D(name + "0_m", 50, 0, 500*GeV);

                mh1D["1_pt"] = bookHisto1D(name + "1_pt", 50, 0, 2000*GeV);
                mh1D["1_eta"] = bookHisto1D(name + "1_eta", 50, -5, 5);
                mh1D["1_phi"] = bookHisto1D(name + "1_phi", 50, 0, 2*PI);
                mh1D["1_E"] = bookHisto1D(name + "1_E", 50, 0, 2000*GeV);
                mh1D["1_m"] = bookHisto1D(name + "1_m", 50, 0, 500*GeV);


                mh1D["01_dr"] = bookHisto1D(name + "01_dr", 50, 0, 8);
                mh1D["01_deta"] = bookHisto1D(name + "01_deta", 50, 0, 8);
                mh1D["01_dphi"] = bookHisto1D(name + "01_dphi", 50, 0, PI);

                mh1D["01_pt"] = bookHisto1D(name + "01_pt", 50, 0, 2000*GeV);
                mh1D["01_eta"] = bookHisto1D(name + "01_eta", 50, -5, 5);
                mh1D["01_phi"] = bookHisto1D(name + "01_phi", 50, 0, 2*PI);
                mh1D["01_E"] = bookHisto1D(name + "01_E", 50, 0, 2000*GeV);
                mh1D["01_m"] = bookHisto1D(name + "01_m", 50, 0, 500*GeV);

                mp1D["_pt_m"] = bookProfile1D(name + "_pt_m", 50, 0, 2000*GeV);
                mp1D["_m_pt"] = bookProfile1D(name + "_m_pt", 50, 0, 500*GeV);

                mp1D["0_pt_m"] = bookProfile1D(name + "0_pt_m", 50, 0, 2000*GeV);
                mp1D["0_m_pt"] = bookProfile1D(name + "0_m_pt", 50, 0, 500*GeV);

                mp1D["1_pt_m"] = bookProfile1D(name + "1_pt_m", 50, 0, 2000*GeV);
                mp1D["1_m_pt"] = bookProfile1D(name + "1_m_pt", 50, 0, 500*GeV);

                mp1D["0_pt_1_pt"] = bookProfile1D(name + "0_pt_1_pt", 50, 0, 2000*GeV);
                mp1D["1_pt_0_pt"] = bookProfile1D(name + "1_pt_0_pt", 50, 0, 2000*GeV);

                mp1D["01_pt_01_dr"] = bookProfile1D(name + "01_pt_01_dr", 50, 0, 2000*GeV);
                mp1D["01_dr_01_pt"] = bookProfile1D(name + "01_dr_01_pt", 50, 0, 5);

                mp1D["01_pt_01_deta"] = bookProfile1D(name + "01_pt_01_deta", 50, 0, 2000*GeV);
                mp1D["01_deta_01_pt"] = bookProfile1D(name + "01_deta_01_pt", 50, 0, 5);

                mp1D["01_pt_01_m"] = bookProfile1D(name + "01_pt_01_m", 50, 0, 2000*GeV);
                mp1D["01_m_01_pt"] = bookProfile1D(name + "01_m_01_pt", 50, 0, 500*GeV);

                return;
            }


            void fillJetCollection(const string &name, const Event &event) {
                MSG_DEBUG("Filling jet collection " << name);

                double weight = event.weight();
                const FastJets &fj =
                    applyProjection<FastJets>(event, name);
                const Jets &jets = fj.jetsByPt(minJetPtCut[name]);

                map<string, Histo1DPtr> &mh1D = jet1DHistos[name];
                map<string, Profile1DPtr> &mp1D = jet1DProfiles[name];

                mh1D["_n"]->fill(jets.size(), weight);

                foreach (const Jet &jet, jets) {
                    mh1D["_pt"]->fill(jet.pT(), weight);
                    mh1D["_eta"]->fill(jet.eta(), weight);
                    mh1D["_phi"]->fill(jet.phi(), weight);
                    mh1D["_E"]->fill(jet.E(), weight);
                    mh1D["_m"]->fill(jet.mass(), weight);

                    mp1D["_pt_m"]->fill(jet.pT(), jet.mass(), weight);
                    mp1D["_m_pt"]->fill(jet.mass(), jet.pT(), weight);
                }

                if (jets.size()) {
                    mh1D["0_pt"]->fill(jets[0].pT(), weight);
                    mh1D["0_eta"]->fill(jets[0].eta(), weight);
                    mh1D["0_phi"]->fill(jets[0].phi(), weight);
                    mh1D["0_E"]->fill(jets[0].energy(), weight);
                    mh1D["0_m"]->fill(jets[0].mass(), weight);

                    mp1D["0_pt_m"]->fill(jets[0].pT(), jets[0].mass(), weight);
                    mp1D["0_m_pt"]->fill(jets[0].mass(), jets[0].pT(), weight);
                }

                if (jets.size() > 1) {
                    mh1D["1_pt"]->fill(jets[1].pT(), weight);
                    mh1D["1_eta"]->fill(jets[1].eta(), weight);
                    mh1D["1_phi"]->fill(jets[1].phi(), weight);
                    mh1D["1_E"]->fill(jets[1].energy(), weight);
                    mh1D["1_m"]->fill(jets[1].mass(), weight);

                    mp1D["1_pt_m"]->fill(jets[1].pT(), jets[1].mass(), weight);
                    mp1D["1_m_pt"]->fill(jets[1].mass(), jets[1].pT(), weight);


                    mh1D["01_dr"]->fill(Rivet::deltaR(jets[0], jets[1]), weight);
                    mh1D["01_deta"]->fill(Rivet::deltaEta(jets[0], jets[1]), weight);
                    mh1D["01_dphi"]->fill(Rivet::deltaPhi(jets[0], jets[1]), weight);

                    FourMomentum p = jets[0].momentum() + jets[1].momentum();
                    mh1D["01_pt"]->fill(p.pT(), weight);
                    mh1D["01_eta"]->fill(p.eta(), weight);
                    mh1D["01_phi"]->fill(p.phi(), weight);
                    mh1D["01_E"]->fill(p.E(), weight);
                    mh1D["01_m"]->fill(p.mass(), weight);


                    mp1D["0_pt_1_pt"]->fill(jets[0].pT(), jets[1].pT(), weight);
                    mp1D["1_pt_0_pt"]->fill(jets[1].pT(), jets[0].pT(), weight);

                    mp1D["01_pt_01_m"]->fill(p.pT(), p.mass(), weight);
                    mp1D["01_m_01_pt"]->fill(p.mass(), p.pT(), weight);

                    mp1D["01_pt_01_dr"]->fill(p.pT(), Rivet::deltaR(jets[0], jets[1]), weight);
                    mp1D["01_dr_01_pt"]->fill(Rivet::deltaR(jets[0], jets[1]), p.pT(), weight);

                    mp1D["01_pt_01_deta"]->fill(p.pT(), Rivet::deltaEta(jets[0], jets[1]), weight);
                    mp1D["01_deta_01_pt"]->fill(Rivet::deltaEta(jets[0], jets[1]), p.pT(), weight);
                }

                return;
            }

    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);

}

#endif
