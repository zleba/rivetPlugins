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

    const vector<double> PtBins = {300, 500, 700, 900, 1100, 1300, 1500, 1700};

    // This analysis is a derived from the class Analysis:
    class bjetsDeltaHL : public Analysis {


        private:

            BinnedHistogram<double> _hist_DeltaPhiAA, _hist_DeltaPhiBA, _hist_DeltaPhiBB;
            BinnedHistogram<double> _hist_DeltaYAA, _hist_DeltaYBA, _hist_DeltaYBB;
            BinnedHistogram<double> _hist_PhiStarAA, _hist_PhiStarBA, _hist_PhiStarBB;


        public:
            // @name Constructors, init, analyze, finalize
            // @{

            // Constructor
            bjetsDeltaHL()
                : Analysis("bjetsDeltaHL") {
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



                for(unsigned i = 0; i < PtBins.size()-1; ++i) {
                    _hist_DeltaPhiAA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y01",i+1), 20, 0, M_PI));
                    _hist_DeltaPhiBA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y01",i+1), 20, 0, M_PI));
                    _hist_DeltaPhiBB.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x03-y01",i+1), 20, 0, M_PI));

                    _hist_DeltaYAA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y02",i+1), 20, 0, 5));
                    _hist_DeltaYBA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y02",i+1), 20, 0, 5));
                    _hist_DeltaYBB.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x03-y02",i+1), 20, 0, 5));

                    _hist_PhiStarAA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y03",i+1), 40, 0, 1));
                    _hist_PhiStarBA.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y03",i+1), 40, 0, 1));
                    _hist_PhiStarBB.addHistogram(PtBins[i], PtBins[i+1],
                                   bookHisto1D(SF("d0%d-x03-y03",i+1), 40, 0, 1));

                }

            }


            double phiStar(const Jet &lminusJ, const Jet &lplusJ)
            {
                const FourMomentum &lminus = lminusJ.momentum();
                const FourMomentum &lplus = lplusJ.momentum();
                const double phi_acop = M_PI - deltaPhi(lminus, lplus);
                const double costhetastar = tanh( 0.5 * (lminus.eta() - lplus.eta()) );
                const double sin2thetastar = (costhetastar > 1) ? 0.0 : (1.0 - sqr(costhetastar));
                const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
                return phistar;
            }


            // Analysis
            void analyze(const Event &event) {

                const double weight = event.weight();      
                const FastJets &fjAK4 = applyProjection<FastJets>(event,"JetsAK4");      
                const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(100*GeV, 7000.0*GeV) && Cuts::absrap < 2.4);


                if(jetsAK4.size() < 2) return;
                auto jet1 = jetsAK4[0];
                auto jet2 = jetsAK4[1];

                for(unsigned i = 0; i < PtBins.size()-1; ++i) {
                    
                    if(jet1.pt() < PtBins[i] || jet1.pt() > PtBins[i+1]) continue;
                    //if(jet2.pt() < PtBins[i] || jet2.pt() > PtBins[i+1]) continue;


                    double dPhi    = deltaPhi(jet1, jet2);
                    double dY      = deltaRap(jet1, jet2);
                    double PhiStar = phiStar (jet1, jet2);

                    //cout <<"pars " << dPhi << " " << dY << " "<< PhiStar << endl;

                    //Any
                    _hist_DeltaPhiAA.fill(jet1.pt(), dPhi, weight);
                    _hist_DeltaYAA.  fill(jet1.pt(),   dY, weight);
                    _hist_PhiStarAA. fill(jet1.pt(),  PhiStar, weight);

                    if(jet1.bTagged() && jet2.bTagged()) {
                        //Both b-jets
                        _hist_DeltaPhiBB.fill(jet1.pt(), dPhi, weight);
                        _hist_DeltaYBB.  fill(jet1.pt(),   dY, weight);
                        _hist_PhiStarBB. fill(jet1.pt(),  PhiStar, weight);
                    }
                    else if(jet1.bTagged() || jet2.bTagged()) {
                        //One b-jet
                        _hist_DeltaPhiBA.fill(jet1.pt(), dPhi, weight);
                        _hist_DeltaYBA.  fill(jet1.pt(),   dY, weight);
                        _hist_PhiStarBA. fill(jet1.pt(),  PhiStar, weight);
                    }
                }


            }

            // Finalize
            void finalize() {
                cout<<"cross Section: "<<crossSection()<<endl;

                _hist_DeltaPhiAA.scale(crossSection()/sumOfWeights(), this);
                _hist_DeltaPhiBA.scale(crossSection()/sumOfWeights(), this);
                _hist_DeltaPhiBB.scale(crossSection()/sumOfWeights(), this);
                _hist_DeltaYAA.scale(crossSection()/sumOfWeights(), this);
                _hist_DeltaYBA.scale(crossSection()/sumOfWeights(), this);
                _hist_DeltaYBB.scale(crossSection()/sumOfWeights(), this);
                _hist_PhiStarAA.scale(crossSection()/sumOfWeights(), this);
                _hist_PhiStarBA.scale(crossSection()/sumOfWeights(), this);
                _hist_PhiStarBB.scale(crossSection()/sumOfWeights(), this);




            }
            //@}


    };

    // This global object acts as a hook for the plugin system. 
    AnalysisBuilder<bjetsDeltaHL> plugin_bjetsDeltaHL;

}
