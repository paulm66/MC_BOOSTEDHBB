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

#include "Rivet/Jet.hh"
#include "Rivet/Projections/FastJets.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/VariableR.hh"

using std::map;
using std::string;
using namespace Rivet::Cuts;


// TODO
// global variables. ick
const string ptlab = "$p_T$ / GeV";
const string mlab = "mass / GeV";
const string drlab = "$\\Delta R$";
const string etalab = "$\\eta$";
const string philab = "$\\phi$";

namespace Rivet {

/// @name Analysis methods
//@{

/// Book histograms and initialise projections before the run
void MC_BOOSTEDHBB::init() {
    // getLog().setLevel(Log::DEBUG);

    bookChannel("ZllBoostedHbb");
    bookChannel("ZllBoostedHb");
    bookChannel("ZllResolvedHbb");
    bookChannel("ZllResolvedHb");
    bookChannel("ZllVRHbb");
    bookChannel("ZllVRHb");

    bookChannel("WlnuBoostedHbb");
    bookChannel("WlnuBoostedHb");
    bookChannel("WlnuResolvedHbb");
    bookChannel("WlnuResolvedHb");
    bookChannel("WlnuVRHbb");
    bookChannel("WlnuVRHb");

    bookChannel("ZnunuBoostedHbb");
    bookChannel("ZnunuBoostedHb");
    bookChannel("ZnunuResolvedHbb");
    bookChannel("ZnunuResolvedHb");
    bookChannel("ZnunuVRHbb");
    bookChannel("ZnunuVRHb");


    ChargedLeptons clfs(FinalState(-2.5, 2.5, 25*GeV));
    addProjection(clfs, "ChargedLeptons");

    // TODO
    // minimum pt cutoff?
    MissingMomentum mmfs(FinalState(-4.2, 4.2, 0*GeV));
    addProjection(mmfs, "MissingMomentum");

    FinalState fs;
    addProjection(ZFinder(fs, etaIn(-2.5, 2.5) & (pT >= 25*GeV), PID::ELECTRON, 75*GeV, 105*GeV), "ZeeFinder");
    addProjection(ZFinder(fs, etaIn(-2.5, 2.5) & (pT >= 25*GeV), PID::MUON, 75*GeV, 105*GeV), "ZmumuFinder");

    addProjection(WFinder(fs, etaIn(-2.5, 2.5) & (pT > 25*GeV), PID::ELECTRON, 65*GeV, 95*GeV, 25*GeV), "WenuFinder");
    addProjection(WFinder(fs, etaIn(-2.5, 2.5) & (pT > 25*GeV), PID::MUON, 65*GeV, 95*GeV, 25*GeV), "WmunuFinder");
    
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


    //adding the caloparts so we can project it out as particles
    addProjection(caloParts, "CaloParts");


    // track jets constituents
    VetoedFinalState trackParts(ChargedFinalState(-2.5, 2.5, 0.5*GeV));
    trackParts.addVetoOnThisFinalState(leptonsAndNeutrinos);

    // register jet collections
    addProjection(FastJets(trackParts, FastJets::ANTIKT, 0.3), "AntiKt03TrackJets");
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 0.4), "AntiKt04CaloJets");
    addProjection(FastJets(caloParts, FastJets::ANTIKT, 1.0), "AntiKt10CaloJets");

    // variable-R jets
    fastjet::JetDefinition::Plugin *vrplug =
        new fastjet::contrib::VariableRPlugin(60*GeV /* rho < mH */, 0.2, 0.6, fastjet::contrib::VariableRPlugin::AKTLIKE);
    addProjection(FastJets(caloParts, vrplug), "AntiKtVRCaloJets");

    fastjet::JetDefinition::Plugin *vrtplug =
        new fastjet::contrib::VariableRPlugin(60*GeV /* rho < mH */, 0.2, 0.6, fastjet::contrib::VariableRPlugin::AKTLIKE);
    addProjection(FastJets(trackParts, vrtplug), "AntiKtVRTrackJets");

    // register Z and W bosons
    bookFourMom("vboson");

    // register special collections
    bookFourMom("higgs");
    bookFourMomPair("vboson_higgs");

    cutflow = bookHisto1D("cutflow", CUTSLEN, 0, CUTSLEN, "cutflow", "cut", "entries");
    jetflow = bookHisto1D("jetflow", 10, 0, 9, "track jets matched", "jets", "entries");

    return;
}


/// Perform the per-event analysis
void MC_BOOSTEDHBB::analyze(const Event& event) {
    const double weight = event.weight();

    // reset cut bits.
    for (unsigned int iCut = 0; iCut < CUTSLEN; ++iCut)
        cutBits[iCut] = false;
    cutBits[NONE] = true;

    // find a boson...
    const Particles& zeebosons =
        applyProjection<ZFinder>(event, "ZeeFinder").bosons();
    const Particles& zmumubosons =
        applyProjection<ZFinder>(event, "ZmumuFinder").bosons();

    const Particles& wenubosons =
        applyProjection<WFinder>(event, "WenuFinder").bosons();
    const Particles& wmunubosons =
        applyProjection<WFinder>(event, "WmunuFinder").bosons();


    const Particles& CalParticles=
        applyProjection<VetoedFinalState>(event, "CaloParts").particles();


    // leptons
    // TODO
    // isolation?
    const Particles& leptons =
        applyProjection<ChargedLeptons>(event, "ChargedLeptons").particles();

    const Particle& mm =
        Particle(0, -applyProjection<MissingMomentum>(event, "MissingMomentum").visibleMomentum());


    const Jets& akt03tjs =
        applyProjection<FastJets>(event, "AntiKt03TrackJets").jetsByPt(25*GeV);
    const Jets& akt03tbjs = bTagged(akt03tjs);

    const Jets& akt04cjs =
        applyProjection<FastJets>(event, "AntiKt04CaloJets").jetsByPt(25*GeV);
    const Jets& akt04cbjs = bTagged(akt04cjs);

    const Jets& akt10cjs =
        applyProjection<FastJets>(event, "AntiKt10CaloJets").jetsByPt(250*GeV);

    const Jets& aktvrcjs =
        applyProjection<FastJets>(event, "AntiKtVRCaloJets").jetsByPt(25*GeV);
    const Jets& aktvrcbjs = bTagged(aktvrcjs);

    const Jets& vrtjs =
        applyProjection<FastJets>(event, "AntiKtVRTrackJets").jetsByPt(25*GeV);




    
    vector<pair<Jet, vector<const Jet* > > > MyTestObj = GhostHunter(CalParticles, vrtjs);
    typedef pair<Jet, vector<const Jet* > > jetAssoc;
    foreach(const jetAssoc &TestObj, MyTestObj){
        jetflow->fill(TestObj.second.size(),1);

    }



    // find vboson
    Particle vboson;
    if (zeebosons.size() && leptons.size() == 2) {
        vboson = zeebosons[0];
        cutBits[ZLL] = true;
    } else if (zmumubosons.size() && leptons.size() == 2) {
        vboson = zmumubosons[0];
        cutBits[ZLL] = true;
    } else if (wenubosons.size() && leptons.size() == 1) {
        vboson = wenubosons[0];
        cutBits[WLNU] = true;
    } else if (wmunubosons.size() && leptons.size() == 1) {
        vboson = wmunubosons[0];
        cutBits[WLNU] = true;
    } else if (leptons.size() == 0 && mm.pT() > 30*GeV) {
        vboson = Particle(23, mm);
        cutBits[ZNUNU] = true;
    }

    cutBits[VBOSON] = cutBits[ZLL] || cutBits[WLNU] || cutBits[ZNUNU];


    // find boosted higgs
    Particle boostedhiggs;

    // very simple boosted higgs tagging
    // exactly one akt10 calo jet
    // all track jets in event must be in the calo jet cone

    if (akt10cjs.size())
        cutBits[ONEAKT10JETINC] = true;

    // exactly one akt10 calo jet
    if (akt10cjs.size() == 1) {
        cutBits[ONEAKT10JETEXC] = true;
        cutBits[BOOSTEDHB] = akt03tbjs.size() == 1;
        cutBits[BOOSTEDHBB] = akt03tbjs.size() >= 2;

        boostedhiggs = Particle(25, akt10cjs[0].mom());

        // all track jets in event must be in the calo jet cone
        foreach (const Jet& tj, akt03tjs) {
            if (Rivet::deltaR(boostedhiggs, tj) > 1.0) {
                cutBits[BOOSTEDHB] = false;
                cutBits[BOOSTEDHBB] = false;
                break;
            }
        }
    }


    // find resolved higgs
    Particle resolvedhiggs;

    if (akt04cjs.size() >= 2)
        cutBits[TWOAKT04JETSINC] = true;

    // very simple resolved higgs tagging
    // exactly two akt04 calo jets, one of which must be b-tagged
    if (akt04cjs.size() == 2) {
        cutBits[TWOAKT04JETSEXC] = true;
        cutBits[RESOLVEDHB] = akt04cbjs.size() == 1;
        cutBits[RESOLVEDHBB] = akt04cbjs.size() == 2;

        resolvedhiggs = Particle(25, akt04cjs[0].mom() + akt04cjs[1].mom());
    }


    Particle vrhiggs;
    if (aktvrcjs.size() >= 2)
        cutBits[TWOAKTVRJETSINC] = true;

    if (aktvrcjs.size() == 2) {
        cutBits[TWOAKTVRJETSEXC] = true;
        cutBits[VRHB] = aktvrcbjs.size() == 1;
        cutBits[VRHBB] = aktvrcbjs.size() == 2;

        vrhiggs = Particle(25, aktvrcjs[0].mom() + aktvrcjs[1].mom());
    }


    if (cutBits[ZLL]) {
        cutBits[ZLLBOOSTEDHB] = cutBits[BOOSTEDHB];
        cutBits[ZLLBOOSTEDHBB] = cutBits[BOOSTEDHBB];
        cutBits[ZLLRESOLVEDHB] = cutBits[RESOLVEDHB];
        cutBits[ZLLRESOLVEDHBB] = cutBits[RESOLVEDHBB];
        cutBits[ZLLVRHB] = cutBits[VRHB];
        cutBits[ZLLVRHBB] = cutBits[VRHBB];
    } else if (cutBits[WLNU]) {
        cutBits[WLNUBOOSTEDHB] = cutBits[BOOSTEDHB];
        cutBits[WLNUBOOSTEDHBB] = cutBits[BOOSTEDHBB];
        cutBits[WLNURESOLVEDHB] = cutBits[RESOLVEDHB];
        cutBits[WLNURESOLVEDHBB] = cutBits[RESOLVEDHBB];
        cutBits[WLNUVRHB] = cutBits[VRHB];
        cutBits[WLNUVRHBB] = cutBits[VRHBB];
    } else if (cutBits[ZNUNU]) {
        cutBits[ZNUNUBOOSTEDHB] = cutBits[BOOSTEDHB];
        cutBits[ZNUNUBOOSTEDHBB] = cutBits[BOOSTEDHBB];
        cutBits[ZNUNURESOLVEDHB] = cutBits[RESOLVEDHB];
        cutBits[ZNUNURESOLVEDHBB] = cutBits[RESOLVEDHBB];
        cutBits[ZNUNUVRHB] = cutBits[VRHB];
        cutBits[ZNUNUVRHBB] = cutBits[VRHBB];
    }


    // fill cuts
    for (int iCut = 0; iCut < CUTSLEN; ++iCut)
        if (cutBits[iCut])
            cutflow->fill(iCut, weight);


    // find channel
    string lepchan;
    if (cutBits[ZLL])
        lepchan = "Zll";
    else if (cutBits[WLNU])
        lepchan = "Wlnu";
    else if (cutBits[ZNUNU])
        lepchan = "Znunu";
    else
        vetoEvent;


    string channel;
    Particle higgs;
    if (cutBits[BOOSTEDHBB]) {
        channel = lepchan + "BoostedHbb";
        higgs = boostedhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } else if (cutBits[BOOSTEDHB]) {
        channel = lepchan + "BoostedHb";
        higgs = boostedhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } 

    if (cutBits[RESOLVEDHBB]) {
        channel = lepchan + "ResolvedHbb";
        higgs = resolvedhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } else if (cutBits[RESOLVEDHB]) {
        channel = lepchan + "ResolvedHb";
        higgs = resolvedhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } 

    if (cutBits[VRHBB]) {
        channel = lepchan + "VRHbb";
        higgs = vrhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } else if (cutBits[VRHB]) {
        channel = lepchan + "VRHb";
        higgs = vrhiggs;
        fillFourMom(channel, "higgs", higgs.mom(), weight);
        fillFourMom(channel, "vboson", vboson, weight);
        fillFourMomPair(channel, "vboson_higgs", vboson.mom(), higgs.mom(), weight);
    } 

    return;
}

//returns a vector of jets that have been constructed from particles and ghost jets 
//the pair gives us an association between each new jet and all of the ghost jets that have been clustered into it
vector<pair<Jet, vector<const Jet* > > >MC_BOOSTEDHBB::GhostHunter(const Particles& parts, const Jets& jets){//ghost busters?

    vector<pair<Jet, vector<const Jet*> > > FinalObject; //final object
    pjs.clear();
    //std::map<int, const Jet*> JetMap;
    vector<const Jet*> JetVec;

    foreach ( const Particle &p, parts){
      fastjet::PseudoJet pj(p.px(), p.py(), p.pz(), p.E()); 
      pj.set_user_index(0);//todo
      pjs.push_back(pj);
    }//add this shit to the header
    int counter=0;
    foreach (const Jet &jet, jets) {
        //so we dont use zero as a user index
        counter++;
        //make jet momentum ~=0 (ghost) and add them to pseudojet in preparation for clustering
        //track them using set_user_index

        JetVec.push_back(&jet);
        //JetMap[counter]= &jet;

        const FourMomentum fv = 1e-20 * jet.momentum();//making a ghost
        fastjet::PseudoJet pj(fv.px(), fv.py(), fv.pz(), fv.E()); 
        pj.set_user_index(-counter);//todo
        pjs.push_back(pj);
    }
    MSG_DEBUG("done adding jets");


    //cluster jet
    //make algorithm type an input to the function?
    //make these variables exist in the header?
    //fastjet::JetDefinition jet_def_small(fastjet::cambridge_algorithm, 0.2);//maybe use cambridge aachen?
    fastjet::JetDefinition jet_def_small(fastjet::cambridge_algorithm, 0.4);//maybe use cambridge aachen?
    fastjet::ClusterSequence clust_seq_small1(pjs, jet_def_small);
    vector<fastjet::PseudoJet> smallJets = fastjet::sorted_by_pt(clust_seq_small1.inclusive_jets(20*GeV));
   

    //now build the output data structure
    foreach (fastjet::PseudoJet newjets, smallJets){
        //clear inner vector
        vector<const Jet*> InnerVector;
        InnerVector.clear();
        //clear pair

        foreach(PseudoJet constituent, newjets.constituents()){
            if (constituent.user_index()<0){
                //map lookup
                //push back map return into inner vector
                
                MSG_DEBUG("found a ghost at " << constituent.user_index() << " size is " << JetVec.size());
                //InnerVector.push_back(JetMap.find(-constituent.user_index()));
                //the user_index starts at -1; the vector starts at zero. --1-1=0
                InnerVector.push_back(JetVec.at(-constituent.user_index()-1));

            }

        }
        pair<Jet, vector<const Jet*> > JetPair;
        JetPair = make_pair(newjets, InnerVector);
        FinalObject.push_back(JetPair);
        //push back into outer vector


        

    }

    return FinalObject;
}


/// Normalise histograms etc., after the run
void MC_BOOSTEDHBB::finalize() {

    // normalize to 1
    double norm = 1000*crossSection()/sumOfWeights();
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


    cutflow->scaleW(norm);


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
        histos1D[chan][name]["pt"] = bookHisto(chan + "_" + name + "_pt", name, ptlab, 25, 0, 2000*GeV);
        histos1D[chan][name]["eta"] = bookHisto(chan + "_" + name + "_eta", name, "$\\eta$", 25, -5, 5);
        histos1D[chan][name]["m"] = bookHisto(chan + "_" + name + "_m", name, mlab, 25, 0, 1000*GeV);

        histos2D[chan][name]["m_vs_pt"] = bookHisto(chan + "_" + name + "_m_vs_pt", name,
                ptlab, 25, 0, 2000*GeV,
                mlab, 25, 0, 1000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomPair(const string& name) {
    // pairs of particles also are "particles"
    bookFourMom(name);

    foreach (const string& chan, channels) {
        // extra histograms for pairs of particles
        histos1D[chan][name]["dr"] = bookHisto(chan + "_" + name + "_dr", drlab, "", 25, 0, 5);

        histos2D[chan][name]["dr_vs_pt"] = bookHisto(chan + "_" + name + "_dr_vs_pt", name,
                ptlab, 25, 0, 2000*GeV,
                drlab, 25, 0, 5);

        histos2D[chan][name]["pt1_vs_pt2"] = bookHisto(chan + "_" + name + "_pt1_vs_pt2", name,
                ptlab, 25, 0, 2000*GeV,
                ptlab, 25, 0, 2000*GeV);

        // pt balance
        histos1D[chan][name]["pt1_minus_pt2"] = bookHisto(chan + "_" + name + "_pt1_minus_pt2", name,
                ptlab, 25, -1000*GeV, 1000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomComp(const string& name) {

    foreach (const string& chan, channels) {
        histos1D[chan][name]["dr"] = bookHisto(chan + "_" + name + "_dr", drlab, name, 25, 0, 0.5);

        histos1D[chan][name]["pt1_minus_pt2"] = bookHisto(chan + "_" + name + "_pt1_minus_pt2", name,
                ptlab, 25, -100*GeV, 100*GeV);

        histos1D[chan][name]["pt1_by_pt2"] = bookHisto(chan + "_" + name + "_pt1_by_pt2", name,
                "$p_{T,1}/p_{T,2}$" , 25, 0, 3);

        histos2D[chan][name]["dr_vs_dpt"] = bookHisto(chan + "_" + name + "_dr_vs_dpt", name,
                ptlab, 25, -100*GeV, 100*GeV,
                drlab, 25, 0, 0.5);

        histos2D[chan][name]["pt1_vs_pt2"] = bookHisto(chan + "_" + name + "_pt1_vs_pt2", name,
                ptlab, 25, 0, 2000*GeV,
                ptlab, 25, 0, 2000*GeV);
    }

    return;
}


void MC_BOOSTEDHBB::bookFourMomColl(const string& name) {
    bookFourMom(name);

    // bookFourMom(name + "0");
    // bookFourMom(name + "1");

    foreach (const string& chan, channels)
        histos1D[chan][name]["n"] = bookHisto(chan + "_" + name + "_n", "multiplicity", "", 10, 0, 10);

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

    return;
}


Jets MC_BOOSTEDHBB::bTagged(const Jets& js) {
    Jets bjs;
    foreach(const Jet& j, js)
        // TODO
        // update for new rivet version
        if (j.bTags().size()) bjs.push_back(j);

    return bjs;
}

//@}

} // Rivet
