// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include <iostream>

using namespace std;

namespace Rivet {


  /// @brief Add a short analysis description here
  class BDIJETS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BDIJETS);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      //declare(FinalState(Cuts::abseta < 4.7), "FS");
      const FinalState cnfs(-4.7, 4.7);
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      //_h_XXXX = bookProfile1D(1, 1, 1);
      _h_DeltaEtaNN = bookHisto1D("DeltaEtaNN", 10, -4, 4);
      _h_DeltaPhiNN = bookHisto1D("DeltaPhiNN", 10, 0, M_PI);

      _h_DeltaEtaBN = bookHisto1D("DeltaEtaBN", 10, -4, 4);
      _h_DeltaPhiBN = bookHisto1D("DeltaPhiBN", 10, 0, M_PI);

      _h_DeltaEtaBB = bookHisto1D("DeltaEtaBB", 10, -4, 4);
      _h_DeltaPhiBB = bookHisto1D("DeltaPhiBB", 10, 0, M_PI);

      //_h_ZZZZ = bookCounter(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        double weight = event.weight();

        const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(200*GeV);
        
        if(jets.size() < 2) return;

        const bool isB1 = jets[0].bTagged();
        const bool isB2 = jets[1].bTagged();

        if(!isB1 && !isB2) { 
            _h_DeltaEtaNN->fill(jets[0].eta() - jets[1].eta(), weight);
            _h_DeltaPhiNN->fill(deltaPhi(jets[0].phi(),jets[1].phi()), weight);
        }
        else if(isB1 && isB2) {
            _h_DeltaEtaBB->fill(jets[0].eta() - jets[1].eta(), weight);
            _h_DeltaPhiBB->fill(deltaPhi(jets[0].phi(),jets[1].phi()), weight);
            //cout << "Both" << endl;
        }
        else {
            //cout << "One" << endl;
            const Jet &bJet = isB1 ? jets[0] : jets[1];
            const Jet &nJet = isB1 ? jets[1] : jets[0];
            _h_DeltaEtaBN->fill(bJet.eta() - nJet.eta(), weight);
            _h_DeltaPhiBN->fill(deltaPhi(bJet.phi(),nJet.phi()), weight);
        }



      /// @todo Do the event by event analysis here

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_YYYY); // normalize to unity
      //scale(_h_ZZZZ, crossSection()/picobarn/sumOfWeights()); // norm to cross section

        const double xsPerWeight = crossSection()/picobarn/sumOfWeights();
        scale({_h_DeltaEtaNN, _h_DeltaPhiNN, _h_DeltaEtaBN, _h_DeltaPhiBN, _h_DeltaEtaBB, _h_DeltaPhiBB}, xsPerWeight);





    }

    //@}


  private:


    /// @name Histograms
    //@{
    Profile1DPtr _h_XXXX;
    CounterPtr _h_ZZZZ;
    Histo1DPtr _h_DeltaEtaNN, _h_DeltaPhiNN, _h_DeltaEtaBN, _h_DeltaPhiBN, _h_DeltaEtaBB, _h_DeltaPhiBB;
    //@}





  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BDIJETS);


}
