// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include <utility>

std::string SF(const char*str, int i)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i);
    return std::string(buff);
}

namespace Rivet {


std::pair<int,int> getNbNbBarLocal(int pid)
{
    int nB = 0, nBbar = 0;
    if (PID::isHadron(pid) && PID::hasBottom(pid))  {
        if(PID::isMeson(pid)) {
            int q1 = (abs(pid) /10) % 10;
            int q2 = (abs(pid) /100) % 10;
            //cout << "AA "<<pid<<" " << q1  << " "<< q2 << " "<< endl;
            if(q1 == 5 && q2 == 5) { //cout << "di-bottom" << endl; 
                ++nB; ++nBbar;
            }
            else {
                if (pid > 0) ++nBbar; //antiB
                else        ++nB; //B
            }
        }
        else if(PID::isBaryon(pid)) {
            if (pid > 0) ++nBbar; //antiB
            else ++nB; //B
        }

    }
    return std::make_pair(nB, nBbar);
}







std::pair<int,int> getNbNbBar(const Jet &jet)
{
    int i = 0;
    int nB = 0, nBbar = 0;
    if (jet.bTagged() ) {
        for( const Particle &p: jet.particles()) {
            if (PID::isHadron(p.pdgId()) && PID::hasBottom(p.pdgId())) {
                auto var = getNbNbBarLocal(p.pdgId());
                nB    += var.first;
                nBbar += var.second;
                //cout << "RR "<< ++i <<" "<< p.pdgId() << endl;
            }

            HepMC::GenVertex* gv = p.genParticle()->production_vertex();
            if (gv) {
                for (const GenParticle* pi  : Rivet::particles(gv, HepMC::ancestors)) {
                    auto var = getNbNbBarLocal(pi->pdg_id());
                    nB    += var.first;
                    nBbar += var.second;
                }
            }

        }
    }
    return std::make_pair(nB, nBbar);

}



    // This analysis is a derived from the class Analysis:
    class bjetsHL : public Analysis {


        private:
            BinnedHistogram<double> _hist_allInclusive;
            BinnedHistogram<double> _hist_bInclusive, _hist_bbBarInclusive, _hist_nobInclusive;


        public:
            // @name Constructors, init, analyze, finalize
            // @{

            // Constructor
            bjetsHL()
                : Analysis("bjetsHL") {
                }

            // Book histograms and initialize projections:
            void init() {

                const FinalState fs;

                //FinalState calofs(Cuts::abseta < 5);
                //addProjection(MissingMomentum(calofs), "CaloMET");

                //VetoedFinalState vfs;
                //addProjection(VisibleFinalState(Cuts::abseta < 5.0), "vfs");

                // Initialize the projectors:
                //addProjection(FastJets(fs, FastJets::ANTIKT, 0.4),"JetsAK4");
                declare(FastJets(fs, FastJets::ANTIKT, 0.4),"JetsAK4");

                //addProjection(HeavyHadrons(Cuts::abseta < 3.5 && Cuts::pT > 1*GeV), "BHadrons");

                std::vector<double> PtHistobinning = {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};



                // Book histograms:

                vector<double> yBins = {0., 0.5, 1, 1.5, 2, 2.4};

                for(int i = 0; i < yBins.size()-1; ++i) {
                    _hist_allInclusive.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y01",i+1), PtHistobinning));
                    _hist_bInclusive.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y02",i+1), PtHistobinning));
                    _hist_bbBarInclusive.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y03",i+1), PtHistobinning));
                    _hist_nobInclusive.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y04",i+1), PtHistobinning));

                }

            }

            // Analysis
            void analyze(const Event &event) {

                const double weight = event.weight();      
                const FastJets &fjAK4 = applyProjection<FastJets>(event,"JetsAK4");      
                const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(18*GeV, 7000.0*GeV) && Cuts::absrap < 4.7);



                int jId = 0;
                for (const Jet& j : jetsAK4) {

                    _hist_allInclusive.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);

                    if(j.bTagged()) {
                        _hist_bInclusive.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);

                        int nB, nBbar;
                        tie(nB, nBbar) = getNbNbBar(j);

                        if(nB != 0 && nBbar != 0)
                            _hist_bbBarInclusive.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);
                        else if(nB == 0 && nBbar == 0)
                            _hist_nobInclusive.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);

                        //if(nB == 0 && nBbar == 0)
                            //cout << "Jet " << jId++ << " " << j.particles().size() << " "<< j.bTagged()<< " "<< j.containsBottom() <<" : "<<  " = " << nB <<" "<< nBbar << endl;
                    }

                }

            }

            // Finalize
            void finalize() {
                cout<<"cross Section: "<<crossSection()<<endl;

                _hist_allInclusive.scale(crossSection()/sumOfWeights()/2, this);
                _hist_bInclusive.scale(crossSection()/sumOfWeights()/2, this);
                _hist_bbBarInclusive.scale(crossSection()/sumOfWeights()/2, this);
                _hist_nobInclusive.scale(crossSection()/sumOfWeights()/2, this);
            }
            //@}


    };

    // This global object acts as a hook for the plugin system. 
    AnalysisBuilder<bjetsHL> plugin_bjetsHL;

}
