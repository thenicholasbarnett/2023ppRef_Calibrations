// imports
#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Timer.h"
#include "JetUncertainty.h"
#include "JetCorrector.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

void save_h1d(TH1D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void save_h1d_1(TH1D *h, TString alg, TString xtitle, TString ytitle, TString hname){
    // making a canvas to make whatever adjustments on
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    // making a copy of the histogram to edit
    TH1D *h_c = (TH1D*)h->Clone(hname);
    h_c->Draw("e1p");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    h_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
    // making a legend for the canvas
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry((TObject*)0, alg, "");
    l->AddEntry(h_c,hname,"pl");
    l->Draw("same");
    // saving the canvas to the file and as a png in the plots folder
    c->Write();
    h->Write();
    c->SaveAs(hname+".png");
    // deleting c and h_c
    delete c;
    delete h_c;
}

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// the script all runs in this function
void pt_balance_generator_ratio_3_26_2024()
{
    // determining the inital time of the script
    Timer timer = Timer();
    timer.Start();
    gStyle->SetOptStat(0);
    
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();

    // creating some binning parameters
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    
    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 7;
    // the low and high pt values for each pt slice
    double ptlow[ptslicenum] = {80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {100,120,140,180,220,300,500};
    double ptedges[ptslicenum+1] = {80,100,120,140,180,220,300,500};
    
    // eta slices
    // number of eta slices
    const Int_t etaslicenum = 8;
    // the low and high eta values for each eta slice
    double etalow[etaslicenum] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9};
    double etahigh[etaslicenum] = {-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};
    double etaedges[etaslicenum+1] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};

    // the only plots we are grabbing from the root files
    TH1D *etaslices_of_ptslicesA[ptslicenum][etaslicenum];
    TH1D *ptslices_of_etaslicesA[etaslicenum][ptslicenum];

    // imported graphs
    // for mc R values tgraphs from files
    TGraph *gMCetaslicesR[etaslicenum];
    TGraph *gMCptslicesR[ptslicenum];
    // for data R values tgraphs from files
    TGraph *gDATAetaslicesR[etaslicenum];
    TGraph *gDATAptslicesR[ptslicenum];

    // R values and their errors
    // MC
    // eta bins of pt slices
    Double_t MCptslicesR[ptslicenum][etaslicenum];
    Double_t MCptslicesRerr[ptslicenum][etaslicenum];
    // pt bins of eta slices
    Double_t MCetaslicesR[etaslicenum][ptslicenum];
    Double_t MCetaslicesRerr[etaslicenum][ptslicenum];
    // DATA
    // eta bins of pt slices
    Double_t DATAptslicesR[ptslicenum][etaslicenum];
    Double_t DATAptslicesRerr[ptslicenum][etaslicenum];
    // pt bins of eta slices
    Double_t DATAetaslicesR[etaslicenum][ptslicenum];
    Double_t DATAetaslicesRerr[etaslicenum][ptslicenum];

    // new hists
    // for mc R values 1d hists
    TH1D *hMCetaslicesR[etaslicenum];
    TH1D *hMCptslicesR[ptslicenum];
    // for data R values 1d hists
    TH1D *hDATAetaslicesR[etaslicenum];
    TH1D *hDATAptslicesR[ptslicenum];
    // for data/mc values 1d hists
    TH1D *hetaslicesR[etaslicenum];
    TH1D *hptslicesR[ptslicenum];

    // setting up the bins for the new hists
    // pt slices have eta bins
    for(unsigned int p=0; p<ptslicenum; p++){
        hMCptslicesR[p] = new TH1D("hMCptslicesR","",etaslicenum,etaedges);
        hDATAptslicesR[p] = new TH1D("hDATAptslicesR","",etaslicenum,etaedges);
    }
    // eta slices have pt bins
    for(unsigned int q=0; q<etaslicenum; q++){
        hMCetaslicesR[q] = new TH1D("hMCetaslicesR","",ptslicenum,ptedges);
        hDATAetaslicesR[q] = new TH1D("hDATAetaslicesR","",ptslicenum,ptedges);
        hetaslicesR[q] = new TH1D("hetaslicesR","",ptslicenum,ptedges);
    }

    // pointing fi to the file holding the jet info of interest
    TFile *fimc = TFile::Open("pt_balance_2023_ppRef_MC_3-26-2024.root","read");
    TFile *fidata = TFile::Open("pt_balance_2023_ppRef_zerobias1_1_3-26-2024.root","read");

    // getting the etaslicesR graphs
    for(unsigned int p=0; p<ptslicenum; p++){
        TString hnamea = Form("pt_%.0f_%.0f",ptlow[p],pthigh[p]);
        gMCptslicesR[p] = (TGraph*)fimc->Get("R_"+hnamea);
        gDATAptslicesR[p] = (TGraph*)fidata->Get("R_"+hnamea);
    }

    // getting the ptslicesR graphs
    for(unsigned int q=0; q<etaslicenum; q++){
        // making the title for the traph
        if(etalow[q]<0){
            TString htitle3 = Form("eta_%.0f_%.0f_n",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
            gMCetaslicesR[q] = (TGraph*)fimc->Get("R_"+htitle3);
            gDATAetaslicesR[q] = (TGraph*)fidata->Get("R_"+htitle3);
        }
        if(etalow[q]>0||etalow[q]==0){
            TString htitle3 = Form("eta_%.0f_%.0f",etalow[q]*10,etahigh[q]*10);
            gMCetaslicesR[q] = (TGraph*)fimc->Get("R_"+htitle3);
            gDATAetaslicesR[q] = (TGraph*)fidata->Get("R_"+htitle3);
        }
    }

    // making the graphs bin values and errors into arrays
    for(unsigned int p=0; p<ptslicenum-1; p++){
        Double_t *ysMCptslicesR = gMCptslicesR[p]->GetY();
        Double_t *ysDATAptslicesR = gDATAptslicesR[p]->GetY();
        for(unsigned int q=0; q<etaslicenum-1; q++){
            Double_t *ysMCetaslicesR = gMCptslicesR[q]->GetY();
            Double_t *ysDATAetaslicesR = gDATAptslicesR[q]->GetY();
            // making arrays of bin values and errors from graphs
            // pt slices
            MCptslicesR[p][q] = ysMCptslicesR[q];
            MCptslicesRerr[p][q] = TMath::Abs(gMCptslicesR[p]->GetErrorY(q));
            DATAptslicesR[p][q] = ysDATAptslicesR[q];
            DATAptslicesRerr[p][q] = TMath::Abs(gDATAptslicesR[p]->GetErrorY(q)); 
            // eta slices
            MCetaslicesR[q][p] = ysMCetaslicesR[p];
            MCetaslicesRerr[q][p] = TMath::Abs(gMCetaslicesR[q]->GetErrorY(p));
            DATAetaslicesR[q][p] = ysDATAetaslicesR[p];
            DATAetaslicesRerr[q][p] = TMath::Abs(gDATAetaslicesR[q]->GetErrorY(p));
        }
    }

    // making hists out of the arrays
    for(unsigned int p=0; p<ptslicenum; p++){
        for(unsigned int q=0; q<etaslicenum; q++){
            // making bins in new hists equal to array values
            // eta slices
            hMCetaslicesR[q]->SetBinContent(p, MCetaslicesR[q][p]);
            hMCetaslicesR[q]->SetBinError(p, MCetaslicesRerr[q][p]);
            hDATAetaslicesR[q]->SetBinContent(p, DATAetaslicesR[q][p]);
            hDATAetaslicesR[q]->SetBinError(p, DATAetaslicesRerr[q][p]);
            // pt slices
            hMCptslicesR[p]->SetBinContent(q, MCetaslicesR[p][q]);
            hMCptslicesR[p]->SetBinError(q, MCetaslicesRerr[p][q]);
            hDATAptslicesR[p]->SetBinContent(q, DATAetaslicesR[p][q]);
            hDATAptslicesR[p]->SetBinError(q, DATAetaslicesRerr[p][q]);
        }
    }

    // making a root file to store stuff
    TFile *f1 = new TFile("pt_balance_ratios_3_26_2024.root", "recreate");
    f1->cd();

    // making the ratios
    // pt slices
    for(unsigned int p=0; p<ptslicenum; p++){
        TString hnamep = Form("pt_slice_%.0f_%.0f",ptlow[p],pthigh[p]);
        hptslicesR[p] = (TH1D*)hMCptslicesR[p]->Clone(hnamep);
        hptslicesR[p]->Divide(hMCptslicesR[p],hDATAptslicesR[p],1,1,"B");
        // save_h1d(hptslicesR[p], hnamep);
        save_h1d_1(hptslicesR[p], "AK4PF", "η", "R^{MC}/R^{DATA}", hnamep);
    }
    // eta slices
    for(unsigned int q=0; q<etaslicenum; q++){
        TString hnameq = Form("eta_slice_%.0f_%.0f",etalow[q],etahigh[q]);
        hetaslicesR[q] = (TH1D*)hMCetaslicesR[q]->Clone(hnameq);
        hetaslicesR[q]->Divide(hMCetaslicesR[q],hDATAetaslicesR[q],1,1,"B");
        // save_h1d(hetaslicesR[q], hnameq);
        save_h1d_1(hetaslicesR[q], "AK4PF", "pT [GeV]", "R^{MC}/R^{DATA}", hnameq);
    }
    f1->Close();
    timer.Stop();
    timer.Report();
}