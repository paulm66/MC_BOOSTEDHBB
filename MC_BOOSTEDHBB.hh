// -*- C++ -*-
#ifndef RIVET_MC_BOOSTEDHBB_HH
#define RIVET_MC_BOOSTEDHBB_HH

#include "Rivet/Analysis.hh"

namespace Rivet {

    class MC_BOOSTEDHBB : public Analysis {
        public:
            /// Constructor
            MC_BOOSTEDHBB()
                : Analysis("MC_BOOSTEDHBB"),
                    cutBits(CUTSLEN, false) {

                    return;
                }

            /// Book histograms and initialise projections before the run
            void init();

            /// Perform the per-event analysis
            void analyze(const Event& event);

            /// Normalise histograms etc., after the run
            void finalize();

        private:


            Histo1DPtr cutflow;
            enum cuts {
                NONE,               // 0
                ZLL,                // 1
                WLNU,               // 2
                ZNUNU,              // 3
                ONEAKT10JET,        // 4
                BOOSTEDHB,          // 5
                BOOSTEDHBB,         // 6
                TWOAKT04JETS,       // 7
                RESOLVEDHB,         // 8
                RESOLVEDHBB,        // 9
                ZLLBOOSTEDHB,       // 10
                ZLLBOOSTEDHBB,      // 11
                ZLLRESOLVEDHB,      // 12
                ZLLRESOLVEDHBB,     // 13
                WLNUBOOSTEDHB,      // 14
                WLNUBOOSTEDHBB,     // 15
                WLNURESOLVEDHB,     // 16
                WLNURESOLVEDHBB,    // 17
                ZNUNUBOOSTEDHB,     // 18
                ZNUNUBOOSTEDHBB,    // 19
                ZNUNURESOLVEDHB,    // 20
                ZNUNURESOLVEDHBB,   // 21
                CUTSLEN
            };

            vector<bool> cutBits;


            vector<string> channels;
            map<string, map<string, map<string, Histo1DPtr> > > histos1D;
            map<string, map<string, map<string, Histo2DPtr> > > histos2D;

            void bookChannel(const string& channel);

            Histo1DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax);

            Histo2DPtr bookHisto(const string& name, const string& title,
                    const string& xlabel, int nxbins, double xmin, double xmax,
                    const string& ylabel, int nybins, double ymin, double ymax);

            void bookFourMom(const string& name);
            void bookFourMomPair(const string& name);
            void bookFourMomComp(const string& name);
            void bookFourMomColl(const string& name);

            void fillFourMom(const string& channel,
                    const string& name,
                    const FourMomentum& p,
                    double weight);

            void fillFourMomPair(const string& channel,
                    const string& name,
                    const FourMomentum& p1,
                    const FourMomentum& p2,
                    double weight);

            void fillFourMomComp(const string& channel,
                    const string& name,
                    const FourMomentum& p1,
                    const FourMomentum& p2,
                    double weight);

            template <class T>
            void fillFourMomColl(const string& channel,
                    const string& name,
                    const vector<T>& ps,
                    double weight);


            Jets bTagged(const Jets& js);
    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_BOOSTEDHBB);
}

#endif
