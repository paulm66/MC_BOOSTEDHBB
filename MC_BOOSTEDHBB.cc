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

                // TODO
                // minimum pt cutoff?
                MissingMomentum mmfs(FinalState(-4.2, 4.2, 0*GeV));
                addProjection(mmfs, "MissingMomentum");

                // calo jets constituents
                // TODO
                // don't include high-pt neutrinos or leptons in jets
                // include electrons?
                IdentifiedFinalState nufs(FinalState(-2.5, 2.5, 25*GeV));
                nufs.acceptNeutrinos();

                MergedFinalState leptonsAndNeutrinos(clfs, nufs);
                addProjection(leptonsAndNeutrinos, "LeptonsAndNeutrinos");

                VetoedFinalState caloParts(FinalState(-4.2, 4.2));
                caloParts.addVetoOnThisFinalState(leptonsAndNeutrinos);

                // "track" jets constituents
                VetoedFinalState trackParts(-2.5, 2.5, 0.5*GeV);
                trackParts.addVetoOnThisFinalState(leptonsAndNeutrinos);

                // register 0.4 calo jets
                string s = "AntiKt04CaloJets";
                jetCollections.push_back(s);
                bookPartCollection(s);
                minJetPtCut[s] = 25*GeV;
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), s);


                // and corresponding b-tagged jets
                s = "AntiKt04CaloJetsB";
                bookPartCollection(s);

                // register 1.0 calo jets
                s = "AntiKt10CaloJets";
                jetCollections.push_back(s);
                bookPartCollection(s);
                minJetPtCut[s] = 250*GeV;
                addProjection(FastJets(caloParts, FastJets::ANTIKT, 1.0), s);

                // and corresponding b-tagged jets
                s = "AntiKt10CaloJetsB";
                bookPartCollection(s);


                // register 0.3 track jets
                s = "AntiKt03TrackJets";
                jetCollections.push_back(s);
                bookPartCollection(s);
                minJetPtCut[s] = 25*GeV;
                addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), s);

                // and corresponding b-tagged jets
                s = "AntiKt03TrackJetsB";
                bookPartCollection(s);

                bookPartCollection("Leptons");
                bookPart("MissingMomentum");

                bookPart("Higgs");
                bookPartCollection("HiggsTrackJets");

                // other histograms
                minLeptonBJet_m = bookHisto1D("MinLeptonBJet_m", 50, 0, 1000*GeV);
                minLeptonBJet_dr = bookHisto1D("MinLeptonBJet_dr", 50, 0, 5);

                return;
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                // leptons
                // TODO
                // isolation?
                const Particles &leptons =
                    applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

                if (!leptons.size())
                    vetoEvent;

                const Particle& mm =
                    Particle(0, -applyProjection<MissingMomentum>(event, "MissingMomentum").visibleMomentum());

                // TODO
                // mass window cuts?
                if (leptons.size() >= 2) {
                    const FourMomentum &zmom = leptons[0].momentum() + leptons[1].momentum();
                    if (abs(zmom.mass() - 90*GeV) > 15*GeV)
                        vetoEvent;
                }


                fillPartCollection("Leptons", leptons, weight);
                fillPart("MissingMomentum", mm, weight);

                foreach (const string &name, jetCollections) {
                    const FastJets &fj =
                        applyProjection<FastJets>(event, name);
                    const Jets &jets = fj.jetsByPt(minJetPtCut[name]);
                    fillPartCollection(name, jets, weight);

                    Jets bjets;
                    foreach (const Jet& jet, jets)
                        if (jet.bTags().size()) bjets.push_back(jet);

                    fillPartCollection(name + "B", bjets, weight);
                }


                // search for higgs: highest-pt akt10 jet matched to two b-tagged track jets
                Jets antiKt10CaloJets =
                    applyProjection<FastJets>(event, "AntiKt10CaloJets").jetsByPt(250*GeV);
                Jets antiKt03TrackJets =
                    applyProjection<FastJets>(event, "AntiKt03TrackJets").jetsByPt(25*GeV);

                Jet higgs;
                Jets matchedTrackJets;
                foreach (const Jet &calojet, antiKt10CaloJets) {
                    matchedTrackJets.clear();

                    // TODO
                    // mass window optimized?
                    if (abs(125*GeV - calojet.mass()) > 30*GeV)
                        continue;

                    foreach (const Jet &trackjet, antiKt03TrackJets) {
                        // is it near the calojet?
                        if (Rivet::deltaR(calojet, trackjet) > 1.0)
                            continue;

                        // is it b-tagged?
                        if (!trackjet.bTags().size())
                            continue;

                        matchedTrackJets.push_back(trackjet);
                    }

                    if (matchedTrackJets.size() >= 2) {
                        higgs = calojet;
                        break;
                    }
                }

                double drMin = -1, massMin = -1;
                if (higgs.pT()) {
                    fillPart("Higgs", higgs, weight);
                    fillPartCollection("HiggsTrackJets", matchedTrackJets, weight);

                    foreach(Jet &trackjet, matchedTrackJets) {
                        if (drMin < 0) {
                            drMin = Rivet::deltaR(trackjet, leptons[0]);
                            massMin = (trackjet.momentum() + leptons[0].momentum()).mass();
                            continue;
                        }

                        double dr = Rivet::deltaR(trackjet, leptons[0]);
                        double mass = (trackjet.momentum() + leptons[0].momentum()).mass();

                        drMin = dr < drMin ? dr : drMin;
                        massMin = mass < massMin ? mass : massMin;
                    }

                    minLeptonBJet_m->fill(massMin, weight);
                    minLeptonBJet_dr->fill(drMin, weight);
                }

                return;
            }


            /// Normalise histograms etc., after the run
            void finalize() {

                double norm = crossSection()/sumOfWeights();
                for (map<string, map<string, Histo1DPtr> >::iterator p = histos1D.begin(); p != histos1D.end(); ++p) {
                    for (map<string, Histo1DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                        q->second->scaleW(norm); // norm to cross section
                    }
                }

                for (map<string, map<string, Histo2DPtr> >::iterator p = histos2D.begin(); p != histos2D.end(); ++p) {
                    for (map<string, Histo2DPtr>::iterator q = p->second.begin(); q != p->second.end(); ++q) {
                        q->second->scaleW(norm); // norm to cross section
                    }
                }

                minLeptonBJet_m->scaleW(norm);
                minLeptonBJet_dr->scaleW(norm);

                return;
            }

            //@}


        private:

            vector<string> jetCollections;
            map<string, map<string, Histo1DPtr> > histos1D;
            map<string, map<string, Histo2DPtr> > histos2D;
            map<string, double> minJetPtCut;

            Histo1DPtr minLeptonBJet_m;
            Histo1DPtr minLeptonBJet_dr;

            void bookPart(const string &name) {
                histos1D[name] = map<string, Histo1DPtr>();
                histos2D[name] = map<string, Histo2DPtr>();

                MSG_DEBUG("Adding particle " << name);

                map<string, Histo1DPtr> &mh1D = histos1D.at(name);
                map<string, Histo2DPtr> &mh2D = histos2D.at(name);

                mh1D["_pt"] = bookHisto1D(name + "_pt", 50, 0, 2000*GeV);
                mh1D["_eta"] = bookHisto1D(name + "_eta", 50, -5, 5);
                mh1D["_phi"] = bookHisto1D(name + "_phi", 50, 0, 2*PI);
                mh1D["_E"] = bookHisto1D(name + "_E", 50, 0, 2000*GeV);
                mh1D["_m"] = bookHisto1D(name + "_m", 50, 0, 500*GeV);

                mh2D["_m_pt"] = bookHisto2D(name + "_m_pt", 50, 0, 2000*GeV, 50, 0, 500*GeV);

                return;
            }

            template <class T>
            void fillPart(const string &name, const T &part, double weight) {
                MSG_DEBUG("Filling particle " << name);

                map<string, Histo1DPtr> &mh1D = histos1D[name];
                map<string, Histo2DPtr> &mh2D = histos2D[name];

                mh1D["_pt"]->fill(part.pT(), weight);
                mh1D["_eta"]->fill(part.eta(), weight);
                mh1D["_phi"]->fill(part.phi(), weight);
                mh1D["_E"]->fill(part.E(), weight);
                mh1D["_m"]->fill(part.mass(), weight);

                mh2D["_m_pt"]->fill(part.pT(), part.mass(), weight);

                return;
            }

            void bookPartCollection(const string &name) {
                histos1D[name] = map<string, Histo1DPtr>();
                histos2D[name] = map<string, Histo2DPtr>();

                MSG_DEBUG("Adding particle collection " << name);

                map<string, Histo1DPtr> &mh1D = histos1D.at(name);
                map<string, Histo2DPtr> &mh2D = histos2D.at(name);

                mh1D["_n"] = bookHisto1D(name + "_n", 10, 0, 10);
                mh1D["_pt"] = bookHisto1D(name + "_pt", 50, 0, 2000*GeV);
                mh1D["_eta"] = bookHisto1D(name + "_eta", 50, -5, 5);
                mh1D["_phi"] = bookHisto1D(name + "_phi", 50, 0, 2*PI);
                mh1D["_E"] = bookHisto1D(name + "_E", 50, 0, 2000*GeV);
                mh1D["_m"] = bookHisto1D(name + "_m", 50, 0, 500*GeV);
                mh1D["_mbypt"] = bookHisto1D(name + "_mbypt", 50, 0, 5);

                mh1D["0_pt"] = bookHisto1D(name + "0_pt", 50, 0, 2000*GeV);
                mh1D["0_eta"] = bookHisto1D(name + "0_eta", 50, -5, 5);
                mh1D["0_phi"] = bookHisto1D(name + "0_phi", 50, 0, 2*PI);
                mh1D["0_E"] = bookHisto1D(name + "0_E", 50, 0, 2000*GeV);
                mh1D["0_m"] = bookHisto1D(name + "0_m", 50, 0, 500*GeV);
                mh1D["0_mbypt"] = bookHisto1D(name + "0_mbypt", 50, 0, 5);

                mh1D["1_pt"] = bookHisto1D(name + "1_pt", 50, 0, 2000*GeV);
                mh1D["1_eta"] = bookHisto1D(name + "1_eta", 50, -5, 5);
                mh1D["1_phi"] = bookHisto1D(name + "1_phi", 50, 0, 2*PI);
                mh1D["1_E"] = bookHisto1D(name + "1_E", 50, 0, 2000*GeV);
                mh1D["1_m"] = bookHisto1D(name + "1_m", 50, 0, 500*GeV);
                mh1D["1_mbypt"] = bookHisto1D(name + "1_mbypt", 50, 0, 5);


                mh1D["01_dr"] = bookHisto1D(name + "01_dr", 50, 0, 8);
                mh1D["01_deta"] = bookHisto1D(name + "01_deta", 50, 0, 8);
                mh1D["01_dphi"] = bookHisto1D(name + "01_dphi", 50, 0, PI);

                mh1D["01_pt"] = bookHisto1D(name + "01_pt", 50, 0, 2000*GeV);
                mh1D["01_eta"] = bookHisto1D(name + "01_eta", 50, -5, 5);
                mh1D["01_phi"] = bookHisto1D(name + "01_phi", 50, 0, 2*PI);
                mh1D["01_E"] = bookHisto1D(name + "01_E", 50, 0, 2000*GeV);
                mh1D["01_m"] = bookHisto1D(name + "01_m", 50, 0, 500*GeV);
                mh1D["01_mbypt"] = bookHisto1D(name + "01_mbypt", 50, 0, 5);

                mh2D["_m_pt"] = bookHisto2D(name + "_m_pt", 50, 0, 2000*GeV, 50, 0, 500*GeV);
                mh2D["0_m_pt"] = bookHisto2D(name + "0_m_pt", 50, 0, 2000*GeV, 50, 0, 500*GeV);
                mh2D["1_m_pt"] = bookHisto2D(name + "1_m_pt", 50, 0, 2000*GeV, 50, 0, 500*GeV);
                mh2D["0_pt_1_pt"] = bookHisto2D(name + "0_pt_1_pt", 50, 0, 2000*GeV, 50, 0, 2000*GeV);
                mh2D["01_dr_01_pt"] = bookHisto2D(name + "01_dr_01_pt", 50, 0, 2000*GeV, 50, 0, 5);
                mh2D["01_deta_01_pt"] = bookHisto2D(name + "01_deta_01_pt", 50, 0, 2000*GeV, 50, 0, 5);
                mh2D["01_m_01_pt"] = bookHisto2D(name + "01_m_01_pt", 50, 0, 2000*GeV, 50, 0, 500*GeV);

                return;
            }


            template <class T>
            void fillPartCollection(const string &name, const vector<T> &parts, double weight) {
                MSG_DEBUG("Filling collection " << name);

                map<string, Histo1DPtr> &mh1D = histos1D[name];
                map<string, Histo2DPtr> &mh2D = histos2D[name];

                mh1D["_n"]->fill(parts.size(), weight);

                foreach (const T &part, parts) {
                    mh1D["_pt"]->fill(part.pT(), weight);
                    mh1D["_eta"]->fill(part.eta(), weight);
                    mh1D["_phi"]->fill(part.phi(), weight);
                    mh1D["_E"]->fill(part.E(), weight);
                    mh1D["_m"]->fill(part.mass(), weight);
                    mh1D["_mbypt"]->fill(part.mass()/part.pT(), weight);

                    mh2D["_m_pt"]->fill(part.pT(), part.mass(), weight);
                }

                if (parts.size()) {
                    mh1D["0_pt"]->fill(parts[0].pT(), weight);
                    mh1D["0_eta"]->fill(parts[0].eta(), weight);
                    mh1D["0_phi"]->fill(parts[0].phi(), weight);
                    mh1D["0_E"]->fill(parts[0].energy(), weight);
                    mh1D["0_m"]->fill(parts[0].mass(), weight);
                    mh1D["0_mbypt"]->fill(parts[0].mass()/parts[0].pT(), weight);

                    mh2D["0_m_pt"]->fill(parts[0].pT(), parts[0].mass(), weight);
                }

                if (parts.size() > 1) {
                    mh1D["1_pt"]->fill(parts[1].pT(), weight);
                    mh1D["1_eta"]->fill(parts[1].eta(), weight);
                    mh1D["1_phi"]->fill(parts[1].phi(), weight);
                    mh1D["1_E"]->fill(parts[1].energy(), weight);
                    mh1D["1_m"]->fill(parts[1].mass(), weight);
                    mh1D["1_mbypt"]->fill(parts[1].mass()/parts[1].pT(), weight);

                    mh2D["1_m_pt"]->fill(parts[1].pT(), parts[1].mass(), weight);


                    mh1D["01_dr"]->fill(Rivet::deltaR(parts[0], parts[1]), weight);
                    mh1D["01_deta"]->fill(Rivet::deltaEta(parts[0], parts[1]), weight);
                    mh1D["01_dphi"]->fill(Rivet::deltaPhi(parts[0], parts[1]), weight);

                    FourMomentum p = parts[0].momentum() + parts[1].momentum();
                    mh1D["01_pt"]->fill(p.pT(), weight);
                    mh1D["01_eta"]->fill(p.eta(), weight);
                    mh1D["01_phi"]->fill(p.phi(), weight);
                    mh1D["01_E"]->fill(p.E(), weight);
                    mh1D["01_m"]->fill(p.mass(), weight);
                    mh1D["01_mbypt"]->fill(p.mass()/p.pT(), weight);


                    mh2D["0_pt_1_pt"]->fill(parts[0].pT(), parts[1].pT(), weight);
                    mh2D["01_m_01_pt"]->fill(p.pT(), p.mass(), weight);
                    mh2D["01_dr_01_pt"]->fill(p.pT(), Rivet::deltaR(parts[0], parts[1]), weight);
                    mh2D["01_deta_01_pt"]->fill(p.pT(), Rivet::deltaEta(parts[0], parts[1]), weight);
                }

                return;
            }

    };



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);

}

#endif
