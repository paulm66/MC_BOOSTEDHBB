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
                : Analysis("MC_BOOSTEDHBB")
            {    }


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
                IdentifiedFinalState nufs(FinalState(-2.5, 2.5, 25*GeV));
                nufs.acceptNeutrinos();

                IdentifiedFinalState mufs(FinalState(-2.5, 2.5, 25*GeV));
                mufs.acceptIdPair(PID::MUON);
                addProjection(mufs, "Muons");

                // don't include high-pt neutrinos or muons in jets
                MergedFinalState muonsAndNeutrinos(mufs, nufs);
                addProjection(muonsAndNeutrinos, "Muons and Neutrinos");

                VetoedFinalState caloParts(FinalState(-4.2, 4.2));
                caloParts.addVetoOnThisFinalState(muonsAndNeutrinos);

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

                return;
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                /// @todo Do the event by event analysis here

                // leptons
                const Particles &leptons =
                    applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

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

                /// @todo Normalise, scale and otherwise manipulate histograms here

                for (map<string, map<string, Histo1DPtr> >::iterator p = jetHists.begin(); p != jetHists.end(); ++p) {
                    for (map<string, Histo1DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                        scale(q->second, crossSection()/sumOfWeights()); // norm to cross section
                    }
                }

                return;
            }

            //@}


        private:

            map<string, map<string, Histo1DPtr> > jetHists;

            void addJetCollection(const string &name) {
                jetHists[name] = map<string, Histo1DPtr>();

                MSG_DEBUG("Adding jet collection " << name);

                map<string, Histo1DPtr> &m = jetHists.at(name);

                m["_n"] = bookHisto1D(name + "_n", 5, 0, 5);
                m["_pt"] = bookHisto1D(name + "_pt", 50, 0, 2000*GeV);
                m["_eta"] = bookHisto1D(name + "_eta", 50, -4.2, 4.2);
                m["_phi"] = bookHisto1D(name + "_phi", 50, -PI, PI);
                m["_E"] = bookHisto1D(name + "_E", 50, 0, 2000*GeV);
                m["_m"] = bookHisto1D(name + "_m", 50, 0, 2000*GeV);

                m["0_pt"] = bookHisto1D(name + "0_pt", 50, 0, 2000*GeV);
                m["0_eta"] = bookHisto1D(name + "0_eta", 50, -4.2, 4.2);
                m["0_phi"] = bookHisto1D(name + "0_phi", 50, -PI, PI);
                m["0_E"] = bookHisto1D(name + "0_E", 50, 0, 2000*GeV);
                m["0_m"] = bookHisto1D(name + "0_m", 50, 0, 2000*GeV);

                m["1_pt"] = bookHisto1D(name + "1_pt", 50, 0, 2000*GeV);
                m["1_eta"] = bookHisto1D(name + "1_eta", 50, -4.2, 4.2);
                m["1_phi"] = bookHisto1D(name + "1_phi", 50, -PI, PI);
                m["1_E"] = bookHisto1D(name + "1_E", 50, 0, 2000*GeV);
                m["1_m"] = bookHisto1D(name + "1_m", 50, 0, 2000*GeV);


                m["_dr01"] = bookHisto1D(name + "_dr01", 50, 0, 8);
                m["_deta01"] = bookHisto1D(name + "_deta01", 50, 0, 8);
                m["_dphi01"] = bookHisto1D(name + "_dphi01", 50, 0, 8);
                m["_m01"] = bookHisto1D(name + "_m01", 50, 0, 2000*GeV);

                return;
            }


            void fillJetCollection(const string &name, const Event &event) {
                MSG_DEBUG("Filling jet collection " << name);

                double weight = event.weight();
                const FastJets &fj =
                    applyProjection<FastJets>(event, name);
                const Jets &jets = fj.jetsByPt();

                map<string, Histo1DPtr> &m = jetHists[name];

                m["_n"]->fill(jets.size(), weight);

                foreach (const Jet &jet, jets) {
                    m["_pt"]->fill(jet.pT(), weight);
                    m["_eta"]->fill(jet.eta(), weight);
                    m["_phi"]->fill(jet.phi(), weight);
                    m["_E"]->fill(jet.E(), weight);
                    m["_m"]->fill(jet.mass(), weight);
                }

                if (jets.size()) {
                    m["0_pt"]->fill(jets[0].pT(), weight);
                    m["0_eta"]->fill(jets[0].eta(), weight);
                    m["0_phi"]->fill(jets[0].phi(), weight);
                    m["0_E"]->fill(jets[0].energy(), weight);
                    m["0_m"]->fill(jets[0].mass(), weight);
                }

                if (jets.size() > 1) {
                    m["1_pt"]->fill(jets[1].pT(), weight);
                    m["1_eta"]->fill(jets[1].eta(), weight);
                    m["1_phi"]->fill(jets[1].phi(), weight);
                    m["1_E"]->fill(jets[1].energy(), weight);
                    m["1_m"]->fill(jets[1].mass(), weight);

                    m["_dr01"]->fill(Rivet::deltaR(jets[0], jets[1]), weight);
                    m["_deta01"]->fill(Rivet::deltaEta(jets[0], jets[1]), weight);
                    m["_dphi01"]->fill(Rivet::deltaPhi(jets[0], jets[1]), weight);
                    m["_m01"]->fill((jets[0].momentum() + jets[1].momentum()).mass(), weight);
                }

                return;
            }


            // Data members like post-cuts event weight counters go here


        private:

            /// @name Histograms
            //@{
            //@}


    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);

}

#endif
