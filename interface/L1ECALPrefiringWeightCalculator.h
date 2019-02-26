#ifndef __L1ECALPREFIRINGWEIGHTCALCULATOR__
#define __L1ECALPREFIRINGWEIGHTCALCULATOR__

#include "bpkFrameWork/bprimeKit/interface/format.h"

#include "TFile.h"
#include "TH2F.h"

#include<string>
#include<map>

class L1ECALPrefiringWeightCalculator{

    public:
        L1ECALPrefiringWeightCalculator(const std::string& L1Maps, const std::string& DataEra, bool UseJetEMPt = false, const double& PrefiringRateSystematicUncty = 0.2);
        ~L1ECALPrefiringWeightCalculator();

        void Calculation (const PhotonInfoBranches& photons, const JetInfoBranches& jets);
        double GetPrefiringWeight(const std::string& type = "Contral") { return _weight[type]; }

    private:

        TFile* _file;
        TH2F* _h_prefmap_photon;
        TH2F* _h_prefmap_jet;
        bool _UseJetEMPt;
        double _PrefiringRateSystematicUncty;
        std::map<std::string, double> _weight;
        
        double GetPrefiringRate(const double& pt, const double& eta, TH2F* h_prefmap, const int& flag);
};

#endif
