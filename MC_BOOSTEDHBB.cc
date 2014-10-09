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
    // getLog().setLevel(Log::DEBUG);


    bookChannel("0l1b");
    bookChannel("0l2b");
    bookChannel("1l1b");
    bookChannel("1l2b");
    bookChannel("2l1b");
    bookChannel("2l2b");


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


    /*
    // register 0.4 calo jets
    bookFourMomColl("akt04c_j");
    bookFourMomColl("akt04c_b");
    bookFourMomPair("akt04c_jj");
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), "AntiKt04CaloJets");
    */

    // register 1.0 calo jets
    // bookFourMomColl("akt10c_j");
    // bookFourMomColl("akt10c_b");
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 1.0), "AntiKt10CaloJets");

    // register 0.3 track jets
    bookFourMomColl("akt03t_j");
    bookFourMomColl("akt03t_b");
    bookFourMomPair("akt03t_jj");
    addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), "AntiKt03TrackJets");

    // register Z and W bosons
    bookFourMom("vboson");

    // register leptons and met
    // bookFourMomColl("leptons");
    // bookFourMom("met");


    // register special collections
    // bookFourMomPair("higgs_tjs");
    bookFourMom("higgs");
    bookFourMomPair("vboson_higgs");

    // bookFourMomComp("bhad_akt03t");
    // bookFourMomComp("bhad_akt04c");

    return;
}


/// Perform the per-event analysis
void MC_BOOSTEDHBB::analyze(const Event& event) {
    const double weight = event.weight();

    // find a boson...
    const Particles& zeebosons =
        applyProjection<ZFinder>(event, "ZeeFinder").bosons();
    const Particles& zmumubosons =
        applyProjection<ZFinder>(event, "ZmumuFinder").bosons();

    const Particles& wenubosons =
        applyProjection<WFinder>(event, "WenuFinder").bosons();
    const Particles& wmunubosons =
        applyProjection<WFinder>(event, "WmunuFinder").bosons();

    // leptons
    // TODO
    // isolation?
    const Particles& leptons =
        applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

    const Particle& mm =
        Particle(0, -applyProjection<MissingMomentum>(event, "MissingMomentum").visibleMomentum());


    // which lepton channel?
    string channel;
    Particle vboson;
    if (zeebosons.size() && leptons.size() == 2) {
        vboson = zeebosons[0];
        channel = "2l";
    } else if (zmumubosons.size() && leptons.size() == 2) {
        vboson = zmumubosons[0];
        channel = "2l";
    } else if (wenubosons.size() && leptons.size() == 1) {
        vboson = wenubosons[0];
        channel = "1l";
    } else if (wmunubosons.size() && leptons.size() == 1) {
        vboson = wmunubosons[0];
        channel = "1l";
    } else if (leptons.size() == 0 && mm.pT() > 30*GeV) {
        vboson = Particle(23, mm);
        channel = "0l";
    } else
        vetoEvent;


    // which nbtags channel?
    const Jets& akt03tjs =
        applyProjection<FastJets>(event, "AntiKt03TrackJets").jetsByPt(20*GeV);

    const Jets& akt03tbjs = bTagged(akt03tjs);

    if (akt03tbjs.size() == 0)
        vetoEvent;
    else if (akt03tbjs.size() == 1)
        channel += "1b";
    else if (akt03tbjs.size() >= 2)
        channel += "2b";


    fillFourMomColl(channel, "akt03t_j", akt03tjs, weight);
    fillFourMomColl(channel, "akt03t_b", akt03tbjs, weight);

    if (akt03tbjs.size() >= 2)
        fillFourMomPair(channel, "akt03t_jj", akt03tbjs[0].mom(), akt03tbjs[1].mom(), weight);
    else if (akt03tbjs.size() && akt03tjs.size())
        fillFourMomPair(channel, "akt03t_jj", akt03tbjs[0].mom(), akt03tjs[0].mom(), weight);
    else if (akt03tjs.size() >= 2)
        fillFourMomPair(channel, "akt03t_jj", akt03tjs[0].mom(), akt03tjs[1].mom(), weight);


    /*
    const Jets& akt04cjs =
        applyProjection<FastJets>(event, "AntiKt04CaloJets").jetsByPt(25*GeV);

    const Jets& akt04cbjs = bTagged(akt04cjs);

    fillFourMomColl(channel, "akt04c_j", akt04cjs, weight);
    fillFourMomColl(channel, "akt04c_b", akt04cbjs, weight);

    if (akt04cbjs.size() >= 2)
        fillFourMomPair(channel, "akt04c_jj", akt04cbjs[0].mom(), akt04cbjs[1].mom(), weight);
    else if (akt04cbjs.size() && akt04cjs.size())
        fillFourMomPair(channel, "akt04c_jj", akt04cbjs[0].mom(), akt04cjs[0].mom(), weight);
    else if (akt04cjs.size() >= 2)
        fillFourMomPair(channel, "akt04c_jj", akt04cjs[0].mom(), akt04cjs[1].mom(), weight);
    */

    const Jets& akt10cjs =
        applyProjection<FastJets>(event, "AntiKt10CaloJets").jetsByPt(250*GeV);

    // const Jets& akt10cbjs = bTagged(akt10cjs);

    // fillFourMomColl(channel, "akt10c_j", akt10cjs, weight);
    // fillFourMomColl(channel, "akt10c_b", akt10cbjs, weight);


    fillFourMom(channel, "vboson", vboson, weight);
    // fillFourMomColl(channel, "leptons", leptons, weight);
    // fillFourMom(channel, "met", mm.mom(), weight);


    // very simple boosted higgs tagging
    // search for higgs:
    // highest-pt akt10 jet
    // matched to two b-tagged track jets
    Jet higgs;
    Jets matchedTrackJets;
    foreach (const Jet& cj, akt10cjs) {
        matchedTrackJets.clear();

        int btags = 0;
        foreach (const Jet& tj, akt03tjs) {
            // is it near the calojet?
            if (Rivet::deltaR(cj, tj) < 1.0)
                matchedTrackJets.push_back(tj);
        }

        // require at least two jets and all b-tagged track jets to be
        // associated with the fat jet
        if (matchedTrackJets.size() >= 2 && btags == akt03tbjs.size()) {
            higgs = cj;
            break;
        }
    }

    if (higgs.pT()) {
        MSG_DEBUG("Higgs candidate found");
        // fillFourMomPair(channel, "higgs_tjs", matchedTrackJets[0].mom(), matchedTrackJets[1].mom(), weight);
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } else
        MSG_DEBUG("No Higgs candidate found");


    /*
    // look at deltaR(b-hadron, jet)
    const Particles& bhads =
        applyProjection<HeavyHadrons>(event, "HeavyHadrons").bHadrons();

    Jet bestjet;
    double drMin = -1;
    foreach (const Particle& bhad, bhads) {

        foreach(const Jet& tj, akt03tjs) {
            if (drMin < 0) {
                drMin = Rivet::deltaR(tj, bhad);
                bestjet = tj;
                continue;
            }

            bestjet = Rivet::deltaR(bhad, tj) < drMin ? tj : bestjet;
        }

        if (drMin > 0)
            fillFourMomComp(channel, "bhad_akt03t", bhad.mom(), bestjet.mom(), weight);

        // reset
        drMin = -1;

        foreach(const Jet& cj, akt04cjs) {
            if (drMin < 0) {
                drMin = Rivet::deltaR(cj, bhad);
                bestjet = cj;
                continue;
            }

            bestjet = Rivet::deltaR(bhad, cj) < drMin ? cj : bestjet;
        }

        if (drMin > 0)
            fillFourMomComp(channel, "bhad_akt04c", bhad.mom(), bestjet.mom(), weight);
    }
    */

    return;
}


/// Normalise histograms etc., after the run
void MC_BOOSTEDHBB::finalize() {

    double norm = crossSection()/sumOfWeights();
    for (map< string, map<string, map<string, Histo1DPtr> > >::iterator p = histos1D.begin(); p != histos1D.end(); ++p) {
        for (map<string, map<string, Histo1DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
            for (map<string, Histo1DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
                r->second->scaleW(norm); // norm to cross section
            }
        }
    }

    for (map< string, map<string, map<string, Histo2DPtr> > >::iterator p = histos2D.begin(); p != histos2D.end(); ++p) {
        for (map<string, map<string, Histo2DPtr> >::iterator q = p->second.begin(); q != p->second.end(); ++q) {
            for (map<string, Histo2DPtr>::iterator r = q->second.begin(); r != q->second.end(); ++r) {
                r->second->scaleW(norm); // norm to cross section
            }
        }
    }

    return;
}


void MC_BOOSTEDHBB::bookChannel(const string& channel) {

    channels.push_back(channel);
    histos1D[channel] = map<string, map<string, Histo1DPtr> >();
    histos2D[channel] = map<string, map<string, Histo2DPtr> >();

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


void MC_BOOSTEDHBB::bookFourMom(const string& name) {
    MSG_DEBUG("Booking " << name << " histograms.");

    foreach (const string& chan, channels) {
        histos1D[chan][name]["pt"] = bookHisto(chan + "_" + name + "_pt", name, ptlab, 50, 0, 2000*GeV);
        histos1D[chan][name]["eta"] = bookHisto(chan + "_" + name + "_eta", name, "$\\eta$", 50, -5, 5);
        histos1D[chan][name]["m"] = bookHisto(chan + "_" + name + "_m", name, mlab, 200, 0, 2000*GeV);

        histos2D[chan][name]["m_vs_pt"] = bookHisto(chan + "_" + name + "_m_vs_pt", name,
                ptlab, 50, 0, 2000*GeV,
                mlab, 200, 0, 2000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomPair(const string& name) {
    // pairs of particles also are "particles"
    bookFourMom(name);

    foreach (const string& chan, channels) {
        // extra histograms for pairs of particles
        histos1D[chan][name]["dr"] = bookHisto(chan + "_" + name + "_dr", drlab, "", 50, 0, 5);

        histos2D[chan][name]["dr_vs_pt"] = bookHisto(chan + "_" + name + "_dr_vs_pt", name,
                ptlab, 50, 0, 2000*GeV,
                drlab, 50, 0, 5);

        histos2D[chan][name]["pt1_vs_pt2"] = bookHisto(chan + "_" + name + "_pt1_vs_pt2", name,
                ptlab, 50, 0, 2000*GeV,
                ptlab, 50, 0, 2000*GeV);

        // pt balance
        histos1D[chan][name]["pt1_minus_pt2"] = bookHisto(chan + "_" + name + "_pt1_minus_pt2", name,
                ptlab, 50, -1000*GeV, 1000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomComp(const string& name) {

    foreach (const string& chan, channels) {
        histos1D[chan][name]["dr"] = bookHisto(chan + "_" + name + "_dr", drlab, name, 50, 0, 0.5);

        histos1D[chan][name]["pt1_minus_pt2"] = bookHisto(chan + "_" + name + "_pt1_minus_pt2", name,
                ptlab, 50, -100*GeV, 100*GeV);

        histos1D[chan][name]["pt1_by_pt2"] = bookHisto(chan + "_" + name + "_pt1_by_pt2", name,
                "$p_{T,1}/p_{T,2}$" , 50, 0, 3);

        histos2D[chan][name]["dr_vs_dpt"] = bookHisto(chan + "_" + name + "_dr_vs_dpt", name,
                ptlab, 50, -100*GeV, 100*GeV,
                drlab, 50, 0, 0.5);

        histos2D[chan][name]["pt1_vs_pt2"] = bookHisto(chan + "_" + name + "_pt1_vs_pt2", name,
                ptlab, 50, 0, 2000*GeV,
                ptlab, 50, 0, 2000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomColl(const string& name) {
    bookFourMom(name);

    // bookFourMom(name + "0");
    // bookFourMom(name + "1");

    foreach (const string& chan, channels) {
        histos1D[chan][name]["n"] = bookHisto(chan + "_" + name + "_n", "multiplicity", "", 10, 0, 10);
    }

    return;
}


void MC_BOOSTEDHBB::fillFourMom(const string& channel, const string& name, const FourMomentum& p, double weight) {
    MSG_DEBUG("Filling " << name << " histograms");

    histos1D[channel][name]["pt"]->fill(p.pT(), weight);
    histos1D[channel][name]["eta"]->fill(p.eta(), weight);
    histos1D[channel][name]["m"]->fill(p.mass(), weight);
    histos2D[channel][name]["m_vs_pt"]->fill(p.pT(), p.mass(), weight);

    return;
}


void MC_BOOSTEDHBB::fillFourMomPair(const string& channel, const string& name, const FourMomentum& p1, const FourMomentum& p2, double weight) {
    fillFourMom(channel, name, p1 + p2, weight);

    double dr = Rivet::deltaR(p1, p2);
    double pt = (p1 + p2).pT();

    histos1D[channel][name]["dr"]->fill(dr, weight);
    histos2D[channel][name]["dr_vs_pt"]->fill(pt, dr);
    histos2D[channel][name]["pt1_vs_pt2"]->fill(p2.pT(), p1.pT());
    histos1D[channel][name]["pt1_minus_pt2"]->fill(p1.pT() - p2.pT());

    return;
}


void MC_BOOSTEDHBB::fillFourMomComp(const string& channel, const string& name, const FourMomentum& p1, const FourMomentum& p2, double weight) {

    double dr = Rivet::deltaR(p1, p2);

    histos1D[channel][name]["dr"]->fill(dr, weight);
    histos1D[channel][name]["pt1_minus_pt2"]->fill(p1.pT() - p2.pT());
    histos1D[channel][name]["pt1_by_pt2"]->fill(p1.pT() / p2.pT());
    histos2D[channel][name]["dr_vs_dpt"]->fill(p1.pT() - p2.pT(), dr, weight);
    histos2D[channel][name]["pt1_vs_pt2"]->fill(p2.pT(), p1.pT());

    return;
}


template <class T>
void MC_BOOSTEDHBB::fillFourMomColl(const string& channel, const string& name, const vector<T>& ps, double weight) {

    MSG_DEBUG("Filling " << ps.size() << " members of collection " << name);
    histos1D[channel][name]["n"]->fill(ps.size(), weight);

    foreach (const T& p, ps)
        fillFourMom(channel, name, p.mom(), weight);

    /*
    if (ps.size())
        fillFourMom(channel, name + "0", ps[0].mom(), weight);

    if (ps.size() >= 2)
        fillFourMom(channel, name + "1", ps[1].mom(), weight);
    */

    return;
}


Jets MC_BOOSTEDHBB::bTagged(const Jets& js) {
    Jets bjs;
    foreach(const Jet& j, js)
        // TODO
        // update
        if (j.bTags().size()) bjs.push_back(j);

    return bjs;
}

//@}

} // Rivet
