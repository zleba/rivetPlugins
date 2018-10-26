// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include <utility>
#include <cstdlib>

/*
template<typename T>
std::string SF(const char*str, T i)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i);
    return std::string(buff);
}
*/

std::string SF(const char*str, const char *i)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i);
    return std::string(buff);
}

std::string SF(const char*str, const char *i, int j)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i, j);
    return std::string(buff);
}
std::string SF(const char*str, const char *i, int j, int k)
{
    char buff[100];
    snprintf(buff, sizeof(buff), str, i, j, k);
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


double isbTagged(double pT, bool isB, bool isC)
{
    double w = 1;

    if(isB) {
        w = 0.5;
    }
    else if(isC) {
        w = 0.1;
    }
    else {
        w = 0.001;
    }

    double r = rand()/(RAND_MAX+0.0);
    
    return (r < w);

}





double bTaggSF(double pT, bool isB, bool isC, int shH, int shL)
{
    double sf = 1;

    if(isB) {
        const double sfB = 0.05;
        sf = 1 + shH*sfB;

    }
    else if(isC) {
        const double sfC = 0.10;
        sf = 1 + shH*sfC;
    }
    else {
        const double sfL = 0.10;
        sf = 1 + shL*sfL;
    }

    return sf;

}

    int nSys = 7;



    // This analysis is a derived from the class Analysis:
    class bjetsNew : public Analysis {


        private:
            //Histo1DPtr _hist_AK4pT, _hist_AK4bpT, _hist_AK4bbpT;
            //Histo1DPtr _hist_AK8pT, _hist_AK8bpT, _hist_AK8bbpT;

            vector<Histo1DPtr> _hist_AK4pT;
            vector<Histo1DPtr> _hist_AK8pT;
            vector<BinnedHistogram<double>> _hist_AK4DeltaPhi;
            vector<BinnedHistogram<double>> _hist_AK8DeltaPhi;


            vector<vector<Histo1DPtr>> _hist_AK4pTRec;
            vector<vector<Histo1DPtr>> _hist_AK8pTRec;
            vector<vector<BinnedHistogram<double>>> _hist_AK4DeltaPhiRec;
            vector<vector<BinnedHistogram<double>>> _hist_AK8DeltaPhiRec;


            //BinnedHistogram<double> _hist_AK4bDeltaPhi;
            //BinnedHistogram<double> _hist_AK4bbDeltaPhi;

            //BinnedHistogram<double> _hist_AK8DeltaPhi;
            //BinnedHistogram<double> _hist_AK8bDeltaPhi;
            //BinnedHistogram<double> _hist_AK8bbDeltaPhi;


        public:
            // @name Constructors, init, analyze, finalize
            // @{

            // Constructor
            bjetsNew()
                : Analysis("bjetsNew") {
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
                declare(FastJets(fs, FastJets::ANTIKT, 0.8),"JetsAK8");

                //addProjection(HeavyHadrons(Cuts::abseta < 3.5 && Cuts::pT > 1*GeV), "BHadrons");


                std::vector<double> EtaHistobinning;
                const int nEta = 18;
                for(int i = 0; i <= nEta; ++i) {
                    double eta = M_PI/2 + (M_PI/2*i)/nEta;
                    EtaHistobinning.push_back(eta);
                }


                std::vector<double> PtHistobinning = {400., 500., 600., 700., 800., 950., 1100., 1250., 1500., 1750., 2000., 2500., 3000.};
                std::vector<double> PtBinsEta = {400., 800., 1200., 1600., 2000};

                std::vector<string> flavType = {"", "b", "bb"};

                _hist_AK4pT.resize(flavType.size());
                _hist_AK8pT.resize(flavType.size());
                _hist_AK4DeltaPhi.resize(flavType.size());
                _hist_AK8DeltaPhi.resize(flavType.size());

                _hist_AK4pTRec.resize(flavType.size());
                _hist_AK8pTRec.resize(flavType.size());
                _hist_AK4DeltaPhiRec.resize(flavType.size());
                _hist_AK8DeltaPhiRec.resize(flavType.size());

                for(unsigned i = 0; i < flavType.size(); ++i) {
                    _hist_AK4pTRec[i].resize(nSys);
                    _hist_AK8pTRec[i].resize(nSys);
                    _hist_AK4DeltaPhiRec[i].resize(nSys);
                    _hist_AK8DeltaPhiRec[i].resize(nSys);
                }


                for(unsigned i = 0; i < flavType.size(); ++i) {
                    const char * f = flavType[i].data();
                    _hist_AK4pT[i] = bookHisto1D(SF("AK4%spT",f), PtHistobinning);
                    _hist_AK8pT[i] = bookHisto1D(SF("AK8%spT",f), PtHistobinning);

                    for(int s = 0; s < nSys; ++s) {
                        _hist_AK4pTRec[i][s] = bookHisto1D(SF("AK4%spTRec%d",f,s), PtHistobinning);
                        _hist_AK8pTRec[i][s] = bookHisto1D(SF("AK8%spTRec%d",f,s), PtHistobinning);
                    }

                    for(unsigned j = 0; j < PtBinsEta.size()-1; ++j) {

                        _hist_AK4DeltaPhi[i].addHistogram(PtBinsEta[j], PtBinsEta[j+1],
                            bookHisto1D(SF("AK4%sDeltaPhi%d",f, j+1), EtaHistobinning));
                        _hist_AK8DeltaPhi[i].addHistogram(PtBinsEta[j], PtBinsEta[j+1],
                            bookHisto1D(SF("AK8%sDeltaPhi%d",f, j+1), EtaHistobinning));

                        for(int s = 0; s < nSys; ++s) {
                            _hist_AK4DeltaPhiRec[i][s].addHistogram(PtBinsEta[j], PtBinsEta[j+1],
                                bookHisto1D(SF("AK4%sDeltaPhi%dRec%d",f, j+1,s), EtaHistobinning));
                            _hist_AK8DeltaPhiRec[i][s].addHistogram(PtBinsEta[j], PtBinsEta[j+1],
                                bookHisto1D(SF("AK8%sDeltaPhi%dRec%d",f, j+1,s), EtaHistobinning));
                        }
                    }
                }

            }


            void FillJets(const Jets& jets, double weight,  vector<Histo1DPtr> &_hist_pT,  vector<BinnedHistogram<double>>  &_hist_DeltaPhi)
            {
                // jets
                if(jets.size() < 2) return; //Inside of PS

                auto jet1 = jets[0];
                auto jet2 = jets[1];

                double pt1 = jet1.momentum().pT();
                double pt2 = jet2.momentum().pT();
                
                double dPhi    = deltaPhi(jet1, jet2);

                //Inclusive 
                _hist_pT[0]->fill(pt1, weight);
                _hist_DeltaPhi[0].fill(pt1, dPhi, weight);
                //cout << "Filling delta " << pt1 <<" "<< dPhi << endl;

                //Both B 
                if(jet1.bTagged() && jet2.bTagged()) {
                    _hist_pT[2]->fill(pt1, weight);
                    _hist_DeltaPhi[2].fill(pt1, dPhi, weight);
                }
                if(jet1.bTagged() || jet2.bTagged()) {
                    if(jet1.bTagged())
                        _hist_pT[1]->fill(pt1, weight);
                    else
                        _hist_pT[1]->fill(pt2, weight);

                    _hist_DeltaPhi[1].fill(pt1, dPhi, weight);
                }

            }

            void FillJetsRec(const Jets& jets, double weight,  vector<vector<Histo1DPtr>> &_hist_pT,  vector<vector<BinnedHistogram<double>>>  &_hist_DeltaPhi)
            {

                // jets
                if(jets.size() < 2) return; //Inside of PS

                auto jet1 = jets[0];
                auto jet2 = jets[1];
                double pt1 = jet1.momentum().pT();
                double pt2 = jet2.momentum().pT();

                double pt1Org = jet1.momentum().pT();
                double pt2Org = jet2.momentum().pT();
                double dPhi    = deltaPhi(jet1, jet2);
                
                bool isB1 = jet1.bTagged(), isC1 = jet1.cTagged();
                bool isB2 = jet2.bTagged(), isC2 = jet2.cTagged();

                bool isTag1 =  isbTagged(pt1, isB1, isC1);
                bool isTag2 =  isbTagged(pt2, isB2, isC2);


                for(int sys = 0; sys < nSys; ++sys) {
                    int shH = 0, shL = 0;
                    double pt1 = pt1Org;
                    double pt2 = pt2Org;

                    if(sys == 1) {
                        shH = 1;
                        shL = 0;
                    }
                    else if(sys == 2) {
                        shH = -1;
                        shL =  0;
                    }
                    else if(sys == 3) {
                        shH = 0;
                        shL = 1;
                    }
                    else if(sys == 4) {
                        shH =  0;
                        shL = -1;
                    }
                    else if(sys == 5) {
                        pt1 *= 1.01;
                        pt2 *= 1.01;
                    }
                    else if(sys == 6) {
                        pt1 *= 0.99;
                        pt2 *= 0.99;
                    }
                    //cout << "sys " << sys << endl;

                    double sf1 = bTaggSF(pt1, isB1, isC1, shH, shL);
                    double sf2 = bTaggSF(pt2, isB2, isC2, shH, shL);

                    //Inclusive 
                    _hist_pT[0][sys]->fill(pt1, weight);
                    _hist_DeltaPhi[0][sys].fill(pt1, dPhi, weight);

                    //Both B 
                    if(isTag1 && isTag2) {
                        _hist_pT[2][sys]->fill(pt1, weight*sf1*sf2);
                        _hist_DeltaPhi[2][sys].fill(pt1, dPhi, weight*sf1*sf2);
                    }
                    if(isTag1 || isTag2) {
                        double sf = 1;
                        if(isTag1) {
                            sf = sf1;
                            _hist_pT[1][sys]->fill(pt1, weight*sf);
                        }
                        else {
                            sf = sf2;
                            _hist_pT[1][sys]->fill(pt2, weight*sf);
                        }
                        _hist_DeltaPhi[1][sys].fill(pt1, dPhi, weight*sf);
                    }

                }

            }



            // Analysis
            void analyze(const Event &event) {


                auto particles  = event.allParticles();
                //for(int i = 0; i < particles.size(); ++i) {
                cout << particles[2].pid() <<" "<< particles[3].pid() << endl;

                const double weight = event.weight();      
                const FastJets &fjAK4 = applyProjection<FastJets>(event,"JetsAK4");
                const FastJets &fjAK8 = applyProjection<FastJets>(event,"JetsAK8");
                const Jets& jetsAK4 = fjAK4.jetsByPt(Cuts::ptIn(400*GeV, 6500.0*GeV) && Cuts::absrap < 2.5);
                const Jets& jetsAK8 = fjAK8.jetsByPt(Cuts::ptIn(400*GeV, 6500.0*GeV) && Cuts::absrap < 2.5);

                //GenLevel
                FillJets(jetsAK4, weight, _hist_AK4pT, _hist_AK4DeltaPhi);
                FillJets(jetsAK8, weight, _hist_AK8pT, _hist_AK8DeltaPhi);

                //RecLevel
                FillJetsRec(jetsAK4, weight, _hist_AK4pTRec, _hist_AK4DeltaPhiRec);
                FillJetsRec(jetsAK8, weight, _hist_AK8pTRec, _hist_AK8DeltaPhiRec);

            }

            // Finalize
            void finalize() {

                double factor = crossSection()/sumOfWeights();

                for(unsigned i = 0; i < _hist_AK4DeltaPhi.size(); ++i) {
                    _hist_AK4DeltaPhi[i].scale(factor, this);
                    _hist_AK8DeltaPhi[i].scale(factor, this);
                    scale(_hist_AK4pT[i],factor);
                    scale(_hist_AK8pT[i],factor);
                }

            }
            //@}


    };

    // This global object acts as a hook for the plugin system. 
    AnalysisBuilder<bjetsNew> plugin_bjetsNew;

}
