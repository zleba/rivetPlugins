// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

std::string SF(const char*str, int i)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i);
    return std::string(buff);
}


namespace Rivet {

    // This analysis is a derived from the class Analysis:
    class CMS_2011_S9086218_BJetak413TeV : public Analysis {


        private:
            BinnedHistogram<double> _hist_sigmaAK4;
            BinnedHistogram<double> _hist_ifleadingsigmaAK4;
            BinnedHistogram<double> _hist_leadingsigmaAK4;

            BinnedHistogram<double> _hist_sigmaAK7;
            BinnedHistogram<double> _hist_ifleadingsigmaAK7;
            BinnedHistogram<double> _hist_leadingsigmaAK7;   

        public:
            // @name Constructors, init, analyze, finalize
            // @{

            // Constructor
            CMS_2011_S9086218_BJetak413TeV()
                : Analysis("CMS_2011_S9086218_BJetak413TeV") {
                }

            // Book histograms and initialize projections:
            void init() {

                const FinalState fs;

                FinalState calofs(Cuts::abseta < 5);
                addProjection(MissingMomentum(calofs), "CaloMET");

                VetoedFinalState vfs;
                addProjection(VisibleFinalState(Cuts::abseta < 5.0), "vfs");

                // Initialize the projectors:
                addProjection(FastJets(fs, FastJets::ANTIKT, 0.4),"JetsAK4");
                addProjection(FastJets(fs, FastJets::ANTIKT, 0.7),"JetsAK7");

                addProjection(HeavyHadrons(Cuts::abseta < 3.5 && Cuts::pT > 1*GeV), "BHadrons");

                double Ptbinning[] = {18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

                std::vector<double> PtHistobinning;//(Ptbinning, );

                for(int i=0;i< sizeof(Ptbinning)/sizeof(Ptbinning[0]);i++){
                    PtHistobinning.push_back(Ptbinning[i]);
                }

                // Book histograms:

                double yBins[] = {0., 0.5, 1, 1.5, 2, 2.4};

                for(int i = 0; i < 5; ++i) {
                    _hist_sigmaAK4.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y01",i+1), PtHistobinning));
                    _hist_leadingsigmaAK4.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y02",i+1), PtHistobinning));
                    _hist_ifleadingsigmaAK4.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x01-y03",i+1), PtHistobinning));

                    _hist_sigmaAK7.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y01",i+1), PtHistobinning));
                    _hist_leadingsigmaAK7.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y02",i+1), PtHistobinning));
                    _hist_ifleadingsigmaAK7.addHistogram(yBins[i], yBins[i+1],
                                   bookHisto1D(SF("d0%d-x02-y03",i+1),PtHistobinning));
                }

                /*
                _hist_sigmaAK4.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y01",PtHistobinning));
                _hist_sigmaAK4.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y01",PtHistobinning));
                _hist_sigmaAK4.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y01",PtHistobinning));
                _hist_sigmaAK4.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y01",PtHistobinning));
                _hist_sigmaAK4.addHistogram(2.0, 2.4, bookHisto1D("d05-x01-y01",PtHistobinning));

                _hist_leadingsigmaAK4.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y02",PtHistobinning));
                _hist_leadingsigmaAK4.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y02",PtHistobinning));
                _hist_leadingsigmaAK4.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y02",PtHistobinning));
                _hist_leadingsigmaAK4.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y02",PtHistobinning));
                _hist_leadingsigmaAK4.addHistogram(2.0, 2.4, bookHisto1D("d05-x01-y02",PtHistobinning));

                _hist_ifleadingsigmaAK4.addHistogram(0.0, 0.5, bookHisto1D("d01-x01-y03",PtHistobinning));
                _hist_ifleadingsigmaAK4.addHistogram(0.5, 1.0, bookHisto1D("d02-x01-y03",PtHistobinning));
                _hist_ifleadingsigmaAK4.addHistogram(1.0, 1.5, bookHisto1D("d03-x01-y03",PtHistobinning));
                _hist_ifleadingsigmaAK4.addHistogram(1.5, 2.0, bookHisto1D("d04-x01-y03",PtHistobinning));
                _hist_ifleadingsigmaAK4.addHistogram(2.0, 2.4, bookHisto1D("d05-x01-y03",PtHistobinning));

                _hist_sigmaAK7.addHistogram(0.0, 0.5, bookHisto1D("d01-x02-y01",PtHistobinning));
                _hist_sigmaAK7.addHistogram(0.5, 1.0, bookHisto1D("d02-x02-y01",PtHistobinning));
                _hist_sigmaAK7.addHistogram(1.0, 1.5, bookHisto1D("d03-x02-y01",PtHistobinning));
                _hist_sigmaAK7.addHistogram(1.5, 2.0, bookHisto1D("d04-x02-y01",PtHistobinning));
                _hist_sigmaAK7.addHistogram(2.0, 2.4, bookHisto1D("d05-x02-y01",PtHistobinning));

                _hist_leadingsigmaAK7.addHistogram(0.0, 0.5, bookHisto1D("d01-x02-y02",PtHistobinning));
                _hist_leadingsigmaAK7.addHistogram(0.5, 1.0, bookHisto1D("d02-x02-y02",PtHistobinning));
                _hist_leadingsigmaAK7.addHistogram(1.0, 1.5, bookHisto1D("d03-x02-y02",PtHistobinning));
                _hist_leadingsigmaAK7.addHistogram(1.5, 2.0, bookHisto1D("d04-x02-y02",PtHistobinning));
                _hist_leadingsigmaAK7.addHistogram(2.0, 2.4, bookHisto1D("d05-x02-y02",PtHistobinning));

                _hist_ifleadingsigmaAK7.addHistogram(0.0, 0.5, bookHisto1D("d01-x02-y03",PtHistobinning));
                _hist_ifleadingsigmaAK7.addHistogram(0.5, 1.0, bookHisto1D("d02-x02-y03",PtHistobinning));
                _hist_ifleadingsigmaAK7.addHistogram(1.0, 1.5, bookHisto1D("d03-x02-y03",PtHistobinning));
                _hist_ifleadingsigmaAK7.addHistogram(1.5, 2.0, bookHisto1D("d04-x02-y03",PtHistobinning));
                _hist_ifleadingsigmaAK7.addHistogram(2.0, 2.4, bookHisto1D("d05-x02-y03",PtHistobinning));
                */


            }

            // Analysis
            void analyze(const Event &event) {

                const double weight = event.weight();      
                const FastJets &fjAK4 = applyProjection<FastJets>(event,"JetsAK4");      
                const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(18*GeV, 7000.0*GeV) && Cuts::absrap < 4.7);

                const FastJets &fjAK7 = applyProjection<FastJets>(event,"JetsAK7");      
                const Jets& jetsAK7 = fjAK7.jetsByPt(Cuts::ptIn(18*GeV, 7000.0*GeV) && Cuts::absrap < 4.7);

                const Particles& bHadrons = applyProjection<HeavyHadrons>(event, "BHadrons").bHadrons();

                int leadingentry=1;
                int btagleadingentry=1;

                //MET cut                                                                                                                                                                                             
                /*Particles vfs_particles = applyProjection<VisibleFinalState>(event, "vfs").particles();

                  FourMomentum pTmiss;
                  FourMomentum pTvisible;
                  foreach ( const Particle & p, vfs_particles ) {
                  pTmiss -= p.momentum();
                  }
                  double eTmiss = pTmiss.pT();

                  const MissingMomentum& mmcalo = applyProjection<MissingMomentum>(event, "CaloMET");*/

                foreach (const Jet& j, jetsAK4) {

                    bool btag=false;
                    bool btagleading=false;

                    int i = 0;
                    if(j.containsBottom() ) {
                        foreach( const GenParticle *p, j.constituents()) {
                            cout << "RR "<< ++i <<" "<< p->pdg_id() << endl;
                        }
                    }
                    /*foreach (const GenParticle* p, particles(event.genEvent())) {

                      int aid = abs(p->pdg_id());
                      if (aid/100 == 5 || aid/1000==5) {
                    // 2J+1 == 1 (mesons) or 2 (baryons)
                    if (aid%10 == 1 || aid%10 == 2) {
                    // No B decaying to B
                    if (aid != 5222 && aid != 5112 && aid != 5212 && aid != 5322) {
                    double DeltaR=(j.momentum().eta()-p->momentum().eta())*(j.momentum().eta()-p->momentum().eta())+deltaPhi(j.momentum().phi(),p->momentum().phi())*deltaPhi(j.momentum().phi(),p->momentum().phi());
                    if(sqrt(DeltaR)<0.4){
                    btag=true;
                    }
                    }
                    }
                    }
                    }*/

                    /*foreach (const GenParticle* p, particles(event.genEvent())) {

                      const PdgId pid = p->pdg_id();
                      if (abs(pid) == 5) { 
                      double difference=(j.momentum().eta()-p->momentum().eta())*(j.momentum().eta()-p->momentum().eta())+deltaPhi(j.momentum().phi(),p->momentum().phi())*deltaPhi(j.momentum().phi(),p->momentum().phi());
                      if(sqrt(difference)<0.3){
                      btag=true;
                      }
                      }
                      }*/

                    foreach (const Particle& b, bHadrons){
                        if (deltaR(j, b) < 0.4) { btag = true; break; }
                    }

                    if(btagleadingentry){
                        btagleadingentry=0;
                        foreach (const Particle& b, bHadrons){
                            if (deltaR(j, b) < 0.4) { btagleading = true; break; }
                        }	  
                    }

                    if (btag){

                        _hist_sigmaAK4.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);

                        if(leadingentry) {
                            _hist_leadingsigmaAK4.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);
                            leadingentry=0;
                        }
                    }

                    if (btagleading){
                        _hist_ifleadingsigmaAK4.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);
                    }


                }

                //AK7

                foreach (const Jet& j, jetsAK7) {

                    bool btag=false;
                    bool btagleading=false;

                    foreach (const Particle& b, bHadrons){
                        if (deltaR(j, b) < 0.7) { btag = true; break; }
                    }

                    if(btagleadingentry){
                        btagleadingentry=0;
                        foreach (const Particle& b, bHadrons){
                            if (deltaR(j, b) < 0.7) { btagleading = true; break; }
                        }	  
                    }

                    if (btag){

                        _hist_sigmaAK7.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);

                        if(leadingentry) {
                            _hist_leadingsigmaAK7.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);
                            leadingentry=0;
                        }
                    }

                    if (btagleading){
                        _hist_ifleadingsigmaAK7.fill(fabs(j.momentum().rapidity()), j.momentum().pT(), weight);
                    }


                }


            }

            // Finalize
            void finalize() {
                cout<<"cross Section: "<<crossSection()<<endl;
                _hist_sigmaAK7.scale(crossSection()/sumOfWeights()/2, this);
                _hist_leadingsigmaAK7.scale(crossSection()/sumOfWeights()/2, this);
                _hist_ifleadingsigmaAK7.scale(crossSection()/sumOfWeights()/2, this);

                _hist_sigmaAK4.scale(crossSection()/sumOfWeights()/2, this);
                _hist_leadingsigmaAK4.scale(crossSection()/sumOfWeights()/2, this);
                _hist_ifleadingsigmaAK4.scale(crossSection()/sumOfWeights()/2, this);
            }
            //@}


    };

    // This global object acts as a hook for the plugin system. 
    AnalysisBuilder<CMS_2011_S9086218_BJetak413TeV> plugin_CMS_2011_S9086218_BJetak413TeV;

}
