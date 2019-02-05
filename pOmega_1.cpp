#include "pOmega_Dimi.h"
#include "pOmega_1.h"

#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"

#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"

double pOmegaSideBand(const double& Momentum, const double* SourcePar, const double* PotPar){
    return SourcePar[0]+Momentum*SourcePar[1]+pow(Momentum,2.)*SourcePar[2]+pow(Momentum,3.)*SourcePar[3];
}

void SimpleFitter1(){

    const double GaussSourceSize = 1.36;
    const unsigned NumMomBins = 25;
    const double kMin = 0;
    const double kMax = 20.*double(NumMomBins);

    const double Mass_p = 938.272;
    const double Mass_Omega = 1672.45;
    const double RedMass =(Mass_p*Mass_Omega)/(Mass_p+Mass_Omega);

    const double lambda0 = 0.85*0.85;
    const double lambdaFeed = 0.85*(0.1+0.04+0.01);
    //side-band
    const double lambdaSB = 0.15;

    printf("The following lambda parameters will be used:\n");
    printf("lambda0 = %.3f\n",lambda0);
    printf("lambdaSB = %.3f\n",lambdaSB);
    printf("lambdaFeed = %.3f\n",lambdaFeed);
    printf("SUM: %.3f\n",lambda0+lambdaSB+lambdaFeed);

    CATS Kitty;
    CATSparameters cPars_pp(CATSparameters::tSource,1,true);
    cPars_pp.SetParameter(0,GaussSourceSize);

    Kitty.SetAnaSource(GaussSource, cPars_pp);
    Kitty.SetUseAnalyticSource(true);

    Kitty.SetExcludeFailedBins(false);
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,1);
    Kitty.SetSpin(1,2);
    Kitty.SetChannelWeight(0, 0);
    Kitty.SetChannelWeight(1, 1);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3334);
    Kitty.SetRedMass( RedMass );

    CATSparameters cPotPars5S2(CATSparameters::tPotential,8,true);
    //cPotPars5S2.SetParameter(0,pOmega_Tetsuo);
    cPotPars5S2.SetParameter(0,pOmega_Lattice);
    cPotPars5S2.SetParameter(1,12);
    //Kitty.SetShortRangePotential(1,0,fDlmPot,cPotPars5S2);

    //AB_pp.SetNotifications(CATS::nSilent);
    Kitty.SetMaxNumThreads(1);
    Kitty.KillTheCat();

    //! DATA FILE
    TFile* DataFile = pOmega_FileName!=""?new TFile(pOmega_FileName, "read"):NULL;
    TH1F* DataHisto = DataFile?(TH1F*)DataFile->Get(pOmega_HistoName):NULL;
    if(!DataHisto){
        printf("\033[1;31mERROR:\033[0m No data was loaded!\n");
        if(DataFile) delete DataFile;
        return;
    }

    TFile* FileReso = new TFile(SigmaMatrixFileName, "read");
    TH2F* hReso_pOmega;
    hReso_pOmega = (TH2F*)FileReso->Get(SigmaMatrixHistoName);
    if(!hReso_pOmega){
        printf("\033[1;31mERROR:\033[0m No resolution matrix was loaded!\n");
        if(DataFile) delete DataFile;
        if(FileReso) delete FileReso;
        return;
    }

    DLM_Ck* Ck_pOmega = new DLM_Ck(1,0,Kitty);
    Ck_pOmega->Update();
    DLM_Ck* Ck_pOmegaSB = new DLM_Ck(4,0,NumMomBins,kMin,kMax,pOmegaSideBand);
    Ck_pOmegaSB->SetSourcePar(0,1.36289);
    Ck_pOmegaSB->SetSourcePar(1,-1.51621*1e-3);
    Ck_pOmegaSB->SetSourcePar(2,2.17304*1e-6);
    Ck_pOmegaSB->SetSourcePar(3,-1.07613*1e-9);
    Ck_pOmegaSB->Update();
//! Correct the number of contributions
    DLM_CkDecomposition CkDec_pOmega("pOmega",2,*Ck_pOmega,hReso_pOmega);
    DLM_CkDecomposition CkDec_pOmegaSB("pOmegaSB",0,*Ck_pOmegaSB,NULL);

    CkDec_pOmega.AddContribution(0,lambdaFeed,DLM_CkDecomposition::cFeedDown);
    CkDec_pOmega.AddContribution(1,lambdaSB,DLM_CkDecomposition::cFake,&CkDec_pOmegaSB);

    CkDec_pOmega.Update();
    CkDec_pOmegaSB.Update();

    DLM_Fitter1* fitter = new DLM_Fitter1(1);

    fitter->SetSystem(0,*DataHisto,1,CkDec_pOmega,kMin,kMax,kMax,kMax);
    fitter->SetSeparateBL(0,false);

    fitter->SetParameter("pOmega",DLM_Fitter1::p_a,1.0,0.7,1.3);
    fitter->FixParameter("pOmega",DLM_Fitter1::p_b,0);
    fitter->FixParameter("pOmega",DLM_Fitter1::p_c,0);
    fitter->FixParameter("pOmega",DLM_Fitter1::p_Cl,-1);

    printf("Fitting in progress...\n");

    fitter->GoBabyGo();

    printf("pp(p_a) = %.3f+/-%.3f\n",fitter->GetParameter("pOmega",DLM_Fitter1::p_a),fitter->GetParError("pOmega",DLM_Fitter1::p_a));
    printf("pp(p_b) = %.3f+/-%.3f\n 1/GeV",fitter->GetParameter("pOmega",DLM_Fitter1::p_b)*1000.,fitter->GetParError("pOmega",DLM_Fitter1::p_b)*1000.);
    printf("pp(p_sor0) = %.3f+/-%.3f\n",fitter->GetParameter("pOmega",DLM_Fitter1::p_sor0),fitter->GetParError("pOmega",DLM_Fitter1::p_sor0));
    printf("chi2/ndf = %.3f/%i = %.3f\n",fitter->GetChi2(),fitter->GetNdf(),fitter->GetChi2Ndf());
    printf("pval = %.3f\n",fitter->GetPval());
    printf("nsigma = %.3f\n",sqrt(2)*TMath::ErfcInverse(fitter->GetPval()));

    TGraph FitResult_pOmega;
    FitResult_pOmega.SetName(TString::Format("FitResult_pOmega"));
    fitter->GetFitGraph(0, FitResult_pOmega);
    TFile* fOut = new TFile(pOmega_OutputFolder+"Output_pOmega.root","recreate");
    FitResult_pOmega.Write();
    DataHisto->Write();

    delete Ck_pOmega;
    delete Ck_pOmegaSB;
    delete fitter;
    if(DataFile) delete DataFile;
    if(FileReso) delete FileReso;
    delete fOut;
}


