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
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

void save_g(TGraph *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void save_h1d(TH1D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void save_h2d(TH2D *h, TString hname){
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
    c->SaveAs("plots/"+hname+".png");
    // deleting c and h_c
    delete c;
    delete h_c;
}

// the script all runs in this function
void pt_balance_ratio_generator_2023ppRef()
{
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();

    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 10;
    // the low and high pt values for each pt slice
    double ptlow[ptslicenum] = {15,25,50,80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {25,50,80,100,120,140,180,220,300,500};

    // eta slices
    const Int_t etaslicenum = 30;
    double etahigh[etaslicenum] = {-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5,5.2};
    double etalow[etaslicenum] = {-5.2,-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5};
    double etaedges[etaslicenum+1] = {-5.2,-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5,5.2};
    
    // cout<< "line 50" <<endl;

    double ah1d0[3] = {100,-1,1};

    // pt balance study arrays
    // pt slices
    // MC
    // R for each pt slice
    Double_t MCptslicesR[ptslicenum][etaslicenum];
    // the uncertainty or error in R for each pt slice
    Double_t MCptslicesRerr[ptslicenum][etaslicenum];
    // DATA
    Double_t DATAptslicesR[ptslicenum][etaslicenum];
    Double_t DATAptslicesRerr[ptslicenum][etaslicenum];
    // DATA & MC
    // eta values on x axis
    Double_t ptslicesx[etaslicenum];
    // making tgraphs with information so the error in x is half the bin width
    Double_t ptslicesxerr[etaslicenum];

    // Ratios //
    Double_t ptslicesR[ptslicenum][etaslicenum];
    Double_t ptslicesRerr[ptslicenum][etaslicenum];

    // cout<< "line 65" <<endl;

    // the only plots we are grabbing from the root files
    TH2D *hMCptslicesA[ptslicenum];
    TH2D *hDATAptslicesA[ptslicenum];

    // the slices of those plots
    TH1D *hMCetasbins_of_ptslicesA[ptslicenum][etaslicenum];
    TH1D *hDATAetasbins_of_ptslicesA[ptslicenum][etaslicenum];

    // new hists
    // for mc R values 1d hists
    TH1D *hMCptslicesR[ptslicenum];
    // for data R values 1d hists
    TH1D *hDATAptslicesR[ptslicenum];
    // for data/mc values 1d hists
    TH1D *hptslicesR[ptslicenum];
    // final graphs of interest
    TGraph *gptslicesR[ptslicenum];

    // pointing fi to the file holding the jet info of interest
    TFile *fimc = TFile::Open("pt_balance_Rval_2023ppRef_MC_4_15_2024.root","read");
    TFile *fidata = TFile::Open("pt_balance_Rval_2023ppRef_DATA_4_15_2024.root","read");

    // making an output root file to store stuff
    TFile *f1 = new TFile("pt_balance_ratios_4_16_2024.root", "recreate");
    f1->cd();

    for(unsigned int p=0; p<ptslicenum; p++){
        // name of hist being retrieved
        TString bhtitle0 = Form("ptslicesA_%.0f_%.0f",ptlow[p],pthigh[p]);
        // getting the 2d hists
        hMCptslicesA[p] = (TH2D*) fimc->Get(bhtitle0);
        hDATAptslicesA[p] = (TH2D*) fidata->Get(bhtitle0);
        save_h2d(hMCptslicesA[p],"MC_"+bhtitle0);
        save_h2d(hDATAptslicesA[p],"DATA_"+bhtitle0);
    }
    
    // further initializing the histograms
    // looping over pt slices
    for(unsigned int p=0; p<ptslicenum; p++){
        for(unsigned int q=0; q<etaslicenum; q++){
            if(q==0){
                // initializing pt slice hists of R vs eta
                TString hname01 = Form("pt_slice_%.0f_%.0f",ptlow[p],pthigh[p]);
                hMCptslicesR[p] = new TH1D("MC"+hname01,"",etaslicenum,etaedges);
                hDATAptslicesR[p] = new TH1D("DATA"+hname01,"",etaslicenum,etaedges);
            }
            // intializing hists that are the bins of the slices
            if(etalow[q]<0){
                // MC
                TString dhtitle0 = Form("MC_Aptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                hMCetasbins_of_ptslicesA[p][q] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                hMCetasbins_of_ptslicesA[p][q] = hMCptslicesA[p]->ProjectionY("",q,q,"");
                save_h1d(hMCetasbins_of_ptslicesA[p][q], dhtitle0);
                // DATA
                TString dhtitle1 = Form("DATA_Aptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                hDATAetasbins_of_ptslicesA[p][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
                hDATAetasbins_of_ptslicesA[p][q] = hDATAptslicesA[p]->ProjectionY("",q,q,"");
                save_h1d(hDATAetasbins_of_ptslicesA[p][q], dhtitle1);
            }
            if(etalow[q]>0||etalow[q]==0){
                // MC
                TString dhtitle0 = Form("MC_Aptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10);
                hMCetasbins_of_ptslicesA[p][q] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                hMCetasbins_of_ptslicesA[p][q] = hMCptslicesA[p]->ProjectionY("",q,q,"");
                save_h1d(hMCetasbins_of_ptslicesA[p][q], dhtitle0);
                // DATA
                TString dhtitle1 = Form("DATA_Aptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10);
                hDATAetasbins_of_ptslicesA[p][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
                hDATAetasbins_of_ptslicesA[p][q] = hDATAptslicesA[p]->ProjectionY("",q,q,"");
                save_h1d(hDATAetasbins_of_ptslicesA[p][q], dhtitle1);
            }
            
            // finding stuff for tgraphs
            // MC
            Double_t MCptaavg = (hMCetasbins_of_ptslicesA[p][q])->GetMean();
            Double_t MCptaavgerr = (hMCetasbins_of_ptslicesA[p][q])->GetMeanError();
            MCptslicesR[p][q] = ((1+MCptaavg)/(1-MCptaavg)); 
            MCptslicesRerr[p][q] = (MCptaavgerr*2/((1-MCptaavg)*(1-MCptaavg)));
            // DATA
            Double_t DATAptaavg = (hDATAetasbins_of_ptslicesA[p][q])->GetMean();
            Double_t DATAptaavgerr = (hDATAetasbins_of_ptslicesA[p][q])->GetMeanError();
            DATAptslicesR[p][q] = ((1+DATAptaavg)/(1-DATAptaavg)); 
            DATAptslicesRerr[p][q] = (DATAptaavgerr*2/((1-DATAptaavg)*(1-DATAptaavg)));

            // making bins in new hists equal to array values
            // pt slices
            hMCptslicesR[p]->SetBinContent(q, MCptslicesR[p][q]);
            hMCptslicesR[p]->SetBinError(q, MCptslicesRerr[p][q]);
            hDATAptslicesR[p]->SetBinContent(q, DATAptslicesR[p][q]);
            hDATAptslicesR[p]->SetBinError(q, DATAptslicesRerr[p][q]);
        }
    }


    // making the ratios
    // pt slices
    for(unsigned int p=0; p<ptslicenum; p++){
        TString hnamep = Form("pt_slice_%.0f_%.0f",ptlow[p],pthigh[p]);
        // saving the Rval hists
        save_h1d(hMCptslicesR[p], "R_MC_"+hnamep);
        save_h1d(hDATAptslicesR[p], "R_DATA_"+hnamep);
        // making the ratio of the Rval hists (R_MC / R_DATA)
        hptslicesR[p] = (TH1D*)hMCptslicesR[p]->Clone();
        hptslicesR[p]->Divide(hMCptslicesR[p],hDATAptslicesR[p],1,1,"B");
        // saving Rval ratio hist
        save_h1d(hptslicesR[p], hnamep);
        save_h1d_1(hptslicesR[p], "AK4PF", "#eta^{Probe}", "R_{MC}/R_{DATA}", hnamep);
    }
    
    // cout<< "line 197" <<endl;

    // making the final tgraphs of interest
    for(unsigned int p=0; p<ptslicenum; p++){
        // declaring arrays to make the tgraph for each tgraph being made
        Double_t ys[etaslicenum];
        Double_t yserr[etaslicenum];
        for(unsigned int q=0; q<etaslicenum; q++){
            // assigning the values of interest to arrays
            // bin q of pt slice p
            // ys is for each p, and q is the bin in the R plot for this p or pt slice
            ys[q] = hptslicesR[p]->GetBinContent(q);
            yserr[q] = hptslicesR[p]->GetBinError(q);
        }
        // making the desired final tgraphs out of the arrays made
        gptslicesR[p] = new TGraphErrors(etaslicenum,ptslicesx,ys,ptslicesxerr,yserr);
        // naming them
        TString hnamea = Form("pt_%.0f_%.0f",ptlow[p],pthigh[p]);
        // saving them
        // save_g(gptslicesR[p],  "R_"+hnamea);
    }
    
    f1->Close();
}