// -*- C++ -*-
#include <cassert>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include <math.h>
using std::to_string;
namespace Rivet {


class b_jets : public Analysis {

    const int nPDFs, nScaleWgts;
    vector<BinnedHistogram<double>> JetInc, JetB, JetN, Jetb, Jetn;

    public:

    b_jets()
        : Analysis("b_jets"), nPDFs(101), nScaleWgts(8)
    { }


    void init() {

        FinalState fs;
        FastJets antikt(fs, FastJets::ANTIKT, 0.4);
        addProjection(antikt, "ANTIKT");
        addProjection(HeavyHadrons(Cuts::abseta < 3.5 && Cuts::pT > 1*GeV), "BHadrons");
        addProjection(FinalPartons(Cuts::abseta < 3.5 && Cuts::pT > 1*GeV), "finalpartons");

        vector<double> Ptbinning {43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116},
            Ybinning {0.,0.5,1.0,1.5,2.0,2.4};

        for (int i = 0; i < nPDFs+nScaleWgts; ++i) {
            BinnedHistogram<double> el_JetInc, el_JetB, el_JetN, el_Jetb, el_Jetn;

            for(size_t y = 0; y < Ybinning.size()-1; ++y) {
                string incName = "ptJetInc_y" + to_string(y) + "_" + to_string(i),
                       BName = "ptJetB_y"   + to_string(y) + "_" + to_string(i),
                       NName = "ptJetN_y"   + to_string(y) + "_" + to_string(i),
                       bName = "ptJetb_y"   + to_string(y) + "_" + to_string(i),
                       nName = "ptJetn_y"   + to_string(y) + "_" + to_string(i);

                el_JetInc.addHistogram(Ybinning[y], Ybinning[y+1], bookHisto1D(incName,Ptbinning));
                el_JetB  .addHistogram(Ybinning[y], Ybinning[y+1], bookHisto1D(  BName,Ptbinning));
                el_JetN  .addHistogram(Ybinning[y], Ybinning[y+1], bookHisto1D(  NName,Ptbinning));
                el_Jetb  .addHistogram(Ybinning[y], Ybinning[y+1], bookHisto1D(  bName,Ptbinning));
                el_Jetn  .addHistogram(Ybinning[y], Ybinning[y+1], bookHisto1D(  nName,Ptbinning));
            }
            JetInc.push_back(el_JetInc);
            JetB  .push_back(el_JetB  );
            JetN  .push_back(el_JetN  );
            Jetb  .push_back(el_Jetb  );
            Jetn  .push_back(el_Jetn  );
        }
    }

    void analyze(const Event& event) {
        //const double weight = event.weight();
        //cout << "----" << endl;

        //cout << this << endl;

        vector<double> weights;
        weights.push_back( event.genEvent()->weights()["nnpdf_1_1"] );

        // PDFs
        for (int i=1;i< nPDFs;i++){
            string Result = to_string(260000+i);
            weights.push_back( event.genEvent()->weights()[Result] );
        }

        // scale variations
        const vector<string> scale_variations { "nnpdf_05_05", "nnpdf_05_1", "nnpdf_1_05", "nnpdf_1_1", "nnpdf_1_2", "nnpdf_2_1", "nnpdf_2_2", "nnpdf_1_1" };
        for (string s: scale_variations)
            weights.push_back(event.genEvent()->weights()[s]);

        static bool check = true;
        if (check) {
            check = false;
            assert((int) weights.size() == nPDFs+nScaleWgts);
            assert(JetInc.size() == weights.size());
            assert(JetB  .size() == weights.size());
            assert(JetN  .size() == weights.size());
            assert(Jetb  .size() == weights.size());
            assert(Jetn  .size() == weights.size());

            for (int i=0; i<nPDFs+nScaleWgts; ++i)
                cout << "weight\t" << 260000+i << '\t' << weights[i]<<endl;
        }

        const Jets& jets = applyProjection<JetAlg>(event, "ANTIKT").jetsByPt();
        const Particles& bHadrons = applyProjection<HeavyHadrons>(event, "BHadrons").bHadrons();
        const Particles& finalpartons = applyProjection<FinalPartons>(event, "finalpartons").particlesByPt();

        foreach (const Jet& jet, jets) {

            bool bflavour = false;
            foreach (const Particle& finalparton, finalpartons){
                if (deltaR(jet, finalparton) < 0.4 && abs(finalparton.pdgId()) == 5) {
                    bflavour = true;
                    //cout << "hurra" << endl;
                    break;
                }
            }
            //foreach (const Particle& particle, jet.particles()){
            //    //cout << "particle ID: " << particle.pdgId() << endl;
            //    if (abs(particle.pdgId()) == 5) {
            //        bflavour = true;
            //        //cout << "hurra" << endl;
            //        break;
            //    }
            //}

            bool Bflavour = false;
            foreach (const Particle& bHadron, bHadrons){
                if (deltaR(jet, bHadron) < 0.4) {
                    Bflavour = true;
                    //cout << "hurra" << endl;
                    break;
                }
            }

            double absrap = jet.momentum().absrapidity(),
                   pt     = jet.momentum().pT();
            //cout << pt << '\t' << absrap << '\t' << bflavour << endl;
            for(int i = 0; i < nPDFs+nScaleWgts; ++i) {
                double w = weights[i];
                //cout << i << '\t' << w << endl;

                // inclusive jet
                JetInc[i].fill(absrap, pt, w);
                // bjet or njet
                (Bflavour ? JetB : JetN)[i].fill(absrap, pt, w);
                (bflavour ? Jetb : Jetn)[i].fill(absrap, pt, w);
            }
        }
    }



    void finalize() {

        //normalize(hist_dn_dphi_with_respectJ1_30_DijetCorr175_180); 
        //normalize(hist_dn_dphi_J1_30_DijetCorr175_180); 

        for(int i = 0; i < nPDFs+nScaleWgts; ++i) {

            JetInc[i].scale(crossSection()/sumOfWeights()/2., this);
            JetB[i].scale(crossSection()/sumOfWeights()/2., this);
            JetN[i].scale(crossSection()/sumOfWeights()/2., this);
            Jetb[i].scale(crossSection()/sumOfWeights()/2., this);
            Jetn[i].scale(crossSection()/sumOfWeights()/2., this);

        }

        /*normalize(hist_DijetCorr175_180);
          normalize(hist_DijetCorr175_180_non_null_30);
          normalize(hist_DijetCorr175_180_null_30);*/
    }



};



    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(b_jets);
    //AnalysisBuilder<Job1> plugin_Job1;
}

