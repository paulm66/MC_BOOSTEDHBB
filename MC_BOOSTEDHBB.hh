// -*- C++ -*-
#ifndef RIVET_MC_BOOSTEDHBB_HH
#define RIVET_MC_BOOSTEDHBB_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    class MC_BOOSTEDHBB : public Analysis {
        public:
            /// Constructor
            MC_BOOSTEDHBB()
                : Analysis("MC_BOOSTEDHBB") {

                    return;
                }

            /// Book histograms and initialise projections before the run
            void init();

            /// Perform the per-event analysis
            void analyze(const Event& event);

            /// Normalise histograms etc., after the run
            void finalize();

        private:

            vector<string> jetColls;
            map<string, double> minJetPtCut;
            map<string, map<string, Histo1DPtr> > histos1D;
            map<string, map<string, Histo2DPtr> > histos2D;

            Histo1DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax);

            Histo2DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel, int nybins, double ymin, double ymax);

            void bookFourMom(const string &name);
            void bookFourMomPair(const string &name);
            void bookFourMomColl(const string &name);

            template <class T>
            void fillFourMom(const string &name, const T &part, double weight);

            template <class T, class S>
            void fillFourMomPair(const string &name,const T &p1, const S &p2, double weight);

            template <class T>
            void fillFourMomColl(const string &name, const vector<T> &parts, double weight);

    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);
}

#endif
