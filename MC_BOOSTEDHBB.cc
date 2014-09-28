// -*- C++ -*-
#include "MC_BOOSTEDHBB.hh"

#include <iostream>
#include <map>
#include <string>

#include "Rivet/Tools/Logging.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"

using std::map;
using std::string;


// TODO
// global variables. ick
const string ptlab = "$p_T$ / GeV";
const string mlab = "mass / GeV";
const string drlab = "$\\Delta R";
const string etalab = "$\\eta";
const string philab = "$\\phi";

namespace Rivet {

/// @name Analysis methods
//@{

/// Book histograms and initialise projections before the run
void MC_BOOSTEDHBB::init() {
    getLog().setLevel(Log::DEBUG);

    ChargedLeptons clfs(FinalState(-2.5, 2.5, 25*GeV));
    addProjection(clfs, "ChargedLeptons");

    // TODO
    // minimum pt cutoff?
    MissingMomentum mmfs(FinalState(-4.2, 4.2, 0*GeV));
    addProjection(mmfs, "MissingMomentum");

    addProjection(HeavyHadrons(), "HeavyHadrons");

    addProjection(ZFinder(-2.5, 2.5, 25*GeV, 11, 75*GeV, 105*GeV), "ZeeFinder");
    addProjection(ZFinder(-2.5, 2.5, 25*GeV, 13, 75*GeV, 105*GeV), "ZmumuFinder");

    addProjection(WFinder(-2.5, 2.5, 25*GeV, 11, 65*GeV, 95*GeV, 25*GeV), "WenuFinder");
    addProjection(WFinder(-2.5, 2.5, 25*GeV, 13, 65*GeV, 95*GeV, 25*GeV), "WmunuFinder");

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
    VetoedFinalState trackParts(ChargedFinalState(-2.5, 2.5, 0.5*GeV));
    trackParts.addVetoOnThisFinalState(leptonsAndNeutrinos);


    // register 0.4 calo jets
    string s = "AntiKt04CaloJets";
    jetColls.push_back(s);
    bookFourMomColl(s);
    minJetPtCut[s] = 25*GeV;
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), s);

    // and corresponding b-tagged jets
    bookFourMomColl("AntiKt04CaloJetsB");


    // register 1.0 calo jets
    s = "AntiKt10CaloJets";
    jetColls.push_back(s);
    bookFourMomColl(s);
    minJetPtCut[s] = 250*GeV;
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 1.0), s);

    // and corresponding b-tagged jets
    bookFourMomColl("AntiKt10CaloJetsB");

    // register 0.3 track jets
    s = "AntiKt03TrackJets";
    jetColls.push_back(s);
    bookFourMomColl(s);
    minJetPtCut[s] = 25*GeV;
    addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), s);

    // and corresponding b-tagged jets
    bookFourMomColl("AntiKt03TrackJetsB");


    // register Z and W bosons
    bookFourMom("ZBoson");
    bookFourMom("WBoson");

    // register leptons and met
    bookFourMomColl("Leptons");
    bookFourMom("MissingMomentum");


    // register special collections
    bookFourMomPair("HiggsTrackJets");
    bookFourMom("HiggsFatJet");

    bookFourMomPair("VBosonHiggs");

    bookFourMomComp("BHadVsNearestAntiKt03TrackJet");
    bookFourMomComp("BHadVsNearestAntiKt04CaloJet");
    bookFourMomComp("BHadVsNearestAntiKt10CaloJet");

    return;
}


/// Perform the per-event analysis
void MC_BOOSTEDHBB::analyze(const Event& event) {
    const double weight = event.weight();

    // find a boson...
    const Particles &zeebosons =
        applyProjection<ZFinder>(event, "ZeeFinder").bosons();
    const Particles &zmumubosons =
        applyProjection<ZFinder>(event, "ZmumuFinder").bosons();

    const Particles &wenubosons =
        applyProjection<WFinder>(event, "WenuFinder").bosons();
    const Particles &wmunubosons =
        applyProjection<WFinder>(event, "WmunuFinder").bosons();

    Particle vboson;
    if (zeebosons.size())
        vboson = zeebosons[0];
    else if (zmumubosons.size())
        vboson = zmumubosons[0];
    else if (wenubosons.size())
        vboson = wenubosons[0];
    else if (wmunubosons.size())
        vboson = wmunubosons[0];
    else
        vetoEvent;

    if (vboson.pid() == 23)
        fillFourMom("ZBoson", vboson, weight);
    else
        fillFourMom("WBoson", vboson, weight);


    // leptons
    // TODO
    // isolation?
    const Particles &leptons =
        applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

    fillFourMomColl("Leptons", leptons, weight);

    const Particle& mm =
        Particle(0, -applyProjection<MissingMomentum>(event, "MissingMomentum").visibleMomentum());

    fillFourMom("MissingMomentum", mm.mom(), weight);


    foreach (const string &name, jetColls) {
        const FastJets &fj =
            applyProjection<FastJets>(event, name);
        const Jets &jets = fj.jetsByPt(minJetPtCut[name]);
        fillFourMomColl(name, jets, weight);

        Jets bjets;
        foreach (const Jet& jet, jets)
            if (jet.bTags().size()) bjets.push_back(jet);

        fillFourMomColl(name + "B", bjets, weight);
    }


    Jets antiKt03TrackJets =
        applyProjection<FastJets>(event, "AntiKt03TrackJets").jetsByPt(25*GeV);

    Jets antiKt04CaloJets =
        applyProjection<FastJets>(event, "AntiKt04CaloJets").jetsByPt(25*GeV);

    Jets antiKt10CaloJets =
        applyProjection<FastJets>(event, "AntiKt10CaloJets").jetsByPt(250*GeV);


    // very simple boosted higgs tagging
    // search for higgs:
    // highest-pt akt10 jet
    // matched to two b-tagged track jets
    // within mass window
    Jet higgs;
    Jets matchedTrackJets;
    bool foundHiggs = false;
    foreach (const Jet &calojet, antiKt10CaloJets) {
        matchedTrackJets.clear();

        if (abs(calojet.mass() - 125*GeV) > 25*GeV)
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
            foundHiggs = true;
            break;
        }
    }

    if (foundHiggs) {
        MSG_DEBUG("Higgs candidate found");
        fillFourMomPair("HiggsTrackJets", matchedTrackJets[0].mom(), matchedTrackJets[1].mom(), weight);
        fillFourMom("HiggsFatJet", higgs.mom(), weight);
        fillFourMomPair("VBosonHiggs", higgs.mom(), vboson.mom(), weight);
    } else
        MSG_DEBUG("No Higgs candidate found");


    // treat leading b-tagged trackjets as higgs
    // book kinematics with leading lepton(s)
    Jets bjets;
    foreach(Jet &trackjet, antiKt03TrackJets) {
        if (trackjet.bTags().size())
            bjets.push_back(trackjet);
    }


    // look at deltaR(b-hadron, jet)
    const Particles& bhads =
        applyProjection<HeavyHadrons>(event, "HeavyHadrons").bHadrons();

    Jet bestjet;
    double drMin = -1;
    foreach (const Particle &bhad, bhads) {

        foreach(Jet &trackjet, antiKt03TrackJets) {
            if (drMin < 0) {
                drMin = Rivet::deltaR(trackjet, bhad);
                bestjet = trackjet;
            }

            bestjet = Rivet::deltaR(bhad, trackjet) < drMin ? trackjet : bestjet;
        }

        if (drMin > 0)
            fillFourMomComp("BHadVsNearestAntiKt03TrackJet", bhad.mom(), bestjet.mom(), weight);

        // reset
        drMin = -1;

        foreach(Jet &calojet, antiKt04CaloJets) {
            if (drMin < 0) {
                drMin = Rivet::deltaR(calojet, bhad);
                bestjet = calojet;
            }

            bestjet = Rivet::deltaR(bhad, calojet) < drMin ? calojet : bestjet;
        }

        if (drMin > 0)
            fillFourMomComp("BHadVsNearestAntiKt04CaloJet", bhad.mom(), bestjet.mom(), weight);

        // reset
        drMin = -1;

        foreach(Jet &calojet, antiKt10CaloJets) {
            if (drMin < 0) {
                drMin = Rivet::deltaR(calojet, bhad);
                bestjet = calojet;
            }

            bestjet = Rivet::deltaR(bhad, calojet) < drMin ? calojet : bestjet;
        }

        fillFourMomComp("BHadVsNearestAntiKt10CaloJet", bhad.mom(), bestjet.mom(), weight);
    }

    return;
}


/// Normalise histograms etc., after the run
void MC_BOOSTEDHBB::finalize() {

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

    return;
}


Histo1DPtr MC_BOOSTEDHBB::bookHisto(const string& name, const string& title,
        const string& xlabel, int nxbins, double xmin, double xmax) {

    double xbinwidth = (xmax - xmin)/nxbins;

    char buff[100];
    sprintf(buff, "events / %.2f", xbinwidth);
    string ylabel = buff;

    return bookHisto1D(name, nxbins, xmin, xmax, title, xlabel, ylabel);
}


Histo2DPtr MC_BOOSTEDHBB::bookHisto(const string& name, const string& title,
        const string& xlabel, int nxbins, double xmin, double xmax,
        const string& ylabel, int nybins, double ymin, double ymax) {

    double xbinwidth = (xmax - xmin)/nxbins;
    double ybinwidth = (ymax - ymin)/nybins;

    char buff[100];
    sprintf(buff, "events / %.2f / %.2f", xbinwidth, ybinwidth);
    string zlabel = buff;

    return bookHisto2D(name, nxbins, xmin, xmax, nybins, ymin, ymax, title, xlabel, ylabel, zlabel);
}


void MC_BOOSTEDHBB::bookFourMom(const string &name) {
    MSG_DEBUG("Booking " << name << " histograms.");

    histos1D[name]["pt"] = bookHisto(name + "_pt", name, ptlab, 50, 0, 2000*GeV);
    histos1D[name]["eta"] = bookHisto(name + "_eta", name, "$\\eta$", 50, -5, 5);
    histos1D[name]["m"] = bookHisto(name + "_m", name, mlab, 50, 0, 500*GeV);

    histos2D[name]["m_vs_pt"] = bookHisto(name + "_m_vs_pt", name,
            ptlab, 50, 0, 2000*GeV,
            mlab, 50, 0, 500*GeV);

    return;
}


void MC_BOOSTEDHBB::bookFourMomPair(const string &name) {
    // pairs of particles also are "particles"
    bookFourMom(name);

    // extra histograms for pairs of particles
    histos1D[name]["dr"] = bookHisto(name + "_dr", drlab, "", 50, 0, 5);

    histos2D[name]["dr_vs_pt"] = bookHisto(name + "_dr_vs_pt", name,
            ptlab, 50, 0, 2000*GeV,
            drlab, 50, 0, 5);

    histos2D[name]["pt1_vs_pt2"] = bookHisto(name + "_pt1_vs_pt2", name,
            ptlab, 50, 0, 2000*GeV,
            ptlab, 50, 0, 2000*GeV);

    return;
}


void MC_BOOSTEDHBB::bookFourMomComp(const string &name) {

    histos1D[name]["dr"] = bookHisto(name + "_dr", drlab, name, 50, 0, 0.5);

    histos1D[name]["pt1_minus_pt2"] = bookHisto(name + "_pt1_minus_pt2", name,
            ptlab, 50, -100*GeV, 100*GeV);

    histos1D[name]["pt1_by_pt2"] = bookHisto(name + "_pt1_by_pt2", name,
            "$p_{T,1}/p_{T,2}$" , 50, 0, 3);

    histos2D[name]["pt1_vs_pt2"] = bookHisto(name + "_pt1_vs_pt2", name,
            ptlab, 50, 0, 2000*GeV,
            ptlab, 50, 0, 2000*GeV);

    return;
}


void MC_BOOSTEDHBB::bookFourMomColl(const string &name) {
    bookFourMom(name);
    bookFourMomPair(name + "LeadingPair");

    histos1D[name]["n"] = bookHisto(name + "_n", "multiplicity", "", 10, 0, 10);

    return;
}


void MC_BOOSTEDHBB::fillFourMom(const string &name, const FourMomentum &part, double weight) {
    MSG_DEBUG("Filling " << name << " histograms");

    histos1D[name]["pt"]->fill(part.pT(), weight);
    histos1D[name]["eta"]->fill(part.eta(), weight);
    histos1D[name]["m"]->fill(part.mass(), weight);
    histos2D[name]["m_vs_pt"]->fill(part.pT(), part.mass(), weight);

    return;
}


void MC_BOOSTEDHBB::fillFourMomPair(const string &name, const FourMomentum &p1, const FourMomentum &p2, double weight) {
    fillFourMom(name, p1 + p2, weight);

    double dr = Rivet::deltaR(p1, p2);
    double pt = (p1 + p2).pT();

    histos1D[name]["dr"]->fill(dr, weight);
    histos2D[name]["dr_vs_pt"]->fill(pt, dr);
    histos2D[name]["pt1_vs_pt2"]->fill(p1.pT(), p2.pT());

    return;
}


void MC_BOOSTEDHBB::fillFourMomComp(const string &name, const FourMomentum &p1, const FourMomentum &p2, double weight) {

    double dr = Rivet::deltaR(p1, p2);

    histos1D[name]["dr"]->fill(dr, weight);
    histos1D[name]["pt1_minus_pt2"]->fill(p1.pT() - p2.pT());
    histos1D[name]["pt1_by_pt2"]->fill(p1.pT() / p2.pT());
    histos2D[name]["pt1_vs_pt2"]->fill(p1.pT(), p2.pT());

    return;
}


template <class T>
void MC_BOOSTEDHBB::fillFourMomColl(const string &name, const vector<T> &parts, double weight) {

    MSG_DEBUG("Filling " << parts.size() << " members of collection " << name);
    histos1D[name]["n"]->fill(parts.size(), weight);

    foreach (const T &part, parts)
        fillFourMom(name, part.mom(), weight);

    if (parts.size() >= 2)
        fillFourMomPair(name + "LeadingPair", parts[0].mom(), parts[1].mom(), weight);

    return;
}

//@}

} // Rivet
