// this code takes in HiForestAOD_11.root, HiForestAOD_186.root, and HiForestAOD_187.root
// and prints out the number of total events in them as well as the number of passed events in them
// where passed events are events with |vz|<15 and have passed the following triggers:
// pPAprimaryVertexFilter, HBHENoiseFilterResultRun2Loose, pBeamScrapingFilter, HLT_HIAK4CaloJet80_v1
// this scripts also prints out how long it takes to run

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

#include "JetUncertainty.h"
#include "JetCorrector.h"

// PLOTTING FUNCTIONS

void ploth1d_1(TH1D *h, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TH1D *h_c = (TH1D*)h->Clone(htitle);
    h_c->Draw("e1p");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
}

void ploth1d_1_s(TH1D *h, TString xtitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TH1D *h_c = (TH1D*)h->Clone(htitle);
    h_c->Draw("e1p");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle("Probability");
    h_c->GetYaxis()->CenterTitle(true);
    h_c->GetXaxis()->SetTitle(xtitle);
    c->SaveAs(htitle);
}

void ploth2d_1(TH2D *h, TString xtitle, TString ytitle, TString htitle, TString option){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TH2D *h_c = (TH2D*)h->Clone(htitle);
    if(option == ""){h_c->Draw("e1p");}
    if(option != ""){h_c->Draw(option);}
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
    c->SaveAs(htitle);

}

void ploth1d_2(TH1D *h0, TString h0name, TH1D *h1, TString h1name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    // cloning h0 and h1
    TH1D *h0_c = (TH1D*)h0->Clone(h1name);
    TH1D *h1_c = (TH1D*)h1->Clone(h1name);
    // h0
    h0_c->Draw("e1p");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->SetTitle(htitle);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    // changing the y axis if needed
    // int maxbin0 = h0_c->GetMaximumBin();
    // int maxbin1 = h1_c->GetMaximumBin();
    // int minbin0 = h0_c->GetMinimumBin();
    // int minbin1 = h1_c->GetMinimumBin();
    // Double_t topmaxbin0 = h0_c->GetYaxis()->GetBinUpEdge(maxbin0);
    // Double_t topmaxbin1 = h1_c->GetYaxis()->GetBinUpEdge(maxbin1);
    // Double_t bottomminbin0 = h0_c->GetYaxis()->GetBinLowEdge(minbin0);
    // Double_t bottomminbin1 = h1_c->GetYaxis()->GetBinLowEdge(minbin1);
    // if(topmaxbin0>topmaxbin1){h0_c->SetMaximum(1.2*topmaxbin0);}
    // if(topmaxbin0<topmaxbin1){h0_c->SetMaximum(1.2*topmaxbin1);}
    // if(bottomminbin0>bottomminbin1){h0_c->SetMinimum(1.2*bottomminbin1);}
    // if(bottomminbin0<bottomminbin1){h0_c->SetMinimum(1.2*bottomminbin0);}
    // setting the x axis title
    if(xtitle == "pt"){
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");    
        c->SetLogy();}
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    h1_c->Draw("e1psame");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->Draw("same");
}

void ploth1d_3(TH1D *h0, TString h0name, TH1D *h1, TString h1name, TH1D *h2, TString h2name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    // h0
    TH1D *h0_c = (TH1D*)h0->Clone(h1name);
    h0_c->Draw("e1p");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->SetTitle(htitle);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");    
        c->SetLogy();}
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    TH1D *h1_c = (TH1D*)h1->Clone(h1name);
    h1_c->Draw("e1psame");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // h1
    TH1D *h2_c = (TH1D*)h2->Clone(h2name);
    h2_c->Draw("e1psame");
    h2_c->SetMarkerStyle(21);
    h2_c->SetMarkerColor(kBlue);
    h2_c->SetLineColor(kBlue);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->AddEntry(h2_c,h2name,"pl");
    l->Draw("same");
}

void plotg1d_1_s(TGraph *h, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TGraph *h_c = (TGraph*)h->Clone(htitle);
    h_c->Draw("ALP");
    h_c->SetTitle("");
    h_c->SetMarkerStyle(20);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}
    c->SaveAs(htitle);
}

void plotg1d_2(TGraph *h0, TString h0name, TGraph *h1, TString h1name, TString xtitle, TString ytitle, TString htitle){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle(htitle);
    TGraph *h0_c = (TGraph*)h0->Clone(htitle);
    TGraph *h1_c = (TGraph*)h1->Clone(htitle);
    // h0
    h0_c->Draw("ap");
    h0_c->SetTitle("");
    h0_c->SetMarkerStyle(20);
    h0_c->SetMarkerColor(kBlack);
    h0_c->SetLineColor(kBlack);
    h0_c->GetXaxis()->CenterTitle(true);
    h0_c->GetYaxis()->SetTitle(ytitle);
    h0_c->GetYaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h0_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h0_c->GetXaxis()->SetTitle(xtitle);}
    // h1
    h1_c->Draw("p same");
    h1_c->SetMarkerStyle(21);
    h1_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kRed);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h0_c,h0name,"pl");
    l->AddEntry(h1_c,h1name,"pl");
    l->Draw("same");

}

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// the script all runs in this function
void Apply_JEC()
{
    // determining the inital time of the script
    clock_t ti = clock();

    // creating the histograms of interest
    TH1::SetDefaultSumw2();

    // creating some binning parameters
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    double vzh1d0[3] = {40,-20,20};
    double pth1d0[3] = {100,100,500};
    double phih1d0[3] = {100,4,4};
    double etah1d0[3] = {50,-5.2,5.2};
    double etah1d1[3] = {25,-1.7,1.7};
    double ah1d0[3] = {10,-1,1};
    // making bins for some pT histograms
    const Int_t bin1 = 5;
    double pth1d1[bin1+1] = {60,80,100,120,200,1000};
    // making a binning for the jer vs refpt plot
    const Int_t bin2 = 12;
    double pth1d2[bin2+1] = {80,90,100,110,120,130,140,150,180,220,300,500};

    // parameters for pt balance responses
    // ptslices
    const Int_t ptslicenum = 7;
    double ptlow[ptslicenum] = {80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {100,120,140,180,220,300,500};
    double ptsliceAtots[ptslicenum];
    double ptsliceAnums[ptslicenum];
    double ptsliceAavgs[ptslicenum];
    TH2D *ptslicesA[ptslicenum];
    TH2D *ptslicesR[ptslicenum];
    // eta slices
    const Int_t etaslicenum = 8;
    double etalow[etaslicenum] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9};
    double etahigh[etaslicenum] = {-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};
    double etasliceAtots[etaslicenum];
    double etasliceAnums[etaslicenum];
    double etasliceAavgs[etaslicenum];
    TH2D *etaslicesA[etaslicenum];
    TH2D *etaslicesR[etaslicenum];
    // slices of eta and pt slices
    TH1D *etaslices_of_ptslicesA[ptslicenum][etaslicenum];
    for(unsigned int q=0; q<ptslicenum; q++){
        ptslicesA[q] = new TH2D(Form("Aptslice_%f_%f",ptlow[q],pthigh[q]),"",etah1d0[0],etah1d0[1],etah1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        ptslicesR[q] = new TH2D(Form("Rptslice_%f_%f",ptlow[q],pthigh[q]),"",etah1d0[0],etah1d0[1],etah1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        for(unsigned int r=0; r<etaslicenum; r++){
            etaslices_of_ptslicesA[q][r] = new TH1D(Form("Aptslice%d_%f_%f__etabin%d_%f_%f",q,ptlow[q],pthigh[q],r,etalow[r],etahigh[r]),"",ah1d0[0],ah1d0[1],ah1d0[2]);
        }
    }
    TH1D *ptslices_of_etaslicesA[etaslicenum][ptslicenum];
    for(unsigned int q=0; q<etaslicenum; q++){
        etaslicesA[q] = new TH2D(Form("Aetaslice_%f_%f",etalow[q],etahigh[q]),"",pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        etaslicesR[q] = new TH2D(Form("Retaslice_%f_%f",etalow[q],etahigh[q]),"",pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        for(unsigned int r=0; r<ptslicenum; r++){
            ptslices_of_etaslicesA[q][r] = new TH1D(Form("Aetaslice%d_%f_%f__ptbin%d_%f_%f",q,etalow[q],etahigh[q],r,ptlow[r],pthigh[r]),"",ah1d0[0],ah1d0[1],ah1d0[2]);
        }
    }

    // established parameters
    TH1D *hvz = new TH1D("vz","",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hjtpt = new TH1D("hjtpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjtcorrpt = new TH1D("hjtpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","",etah1d1[0],etah1d1[1],etah1d1[2]);
    TH1D *hjtphi = new TH1D("hjtphi","",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *hrawpt = new TH1D("hrawpt","",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta_uc = new TH1D("hjteta_uc","",etah1d0[0],etah1d0[1],etah1d0[2]);
    
    // a_, b_ and c_ are the number of total events, accepted events before and after eata cut respectively
    Float_t a_=0, b_=0, c_=0;

    // pointing fi0 and fi1 to the files holding the data of interest
    TFile *fi = TFile::Open("HiForestMiniAOD_10k.root","read");
    // TFile *fi = TFile::Open("HP0_1_25_2024.root","read");

    // declaring variables
    int ppVF, HLT_AKCJ60v1, HLT_AKCJ80v1;
    Float_t vz;
    const Int_t nm = 200000;
    Float_t jtpt[nm];
    Float_t jtcorrpt[nm];
    Float_t jtphi[nm];
    Float_t rawpt[nm];
    Float_t jteta[nm];
    Int_t nref;
    Float_t atot, anum, aavg;
    Float_t tageta = 1.3; 
    Float_t probeeta = 2.5;

    // getting the TTrees from the file
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");
    // TTree *t0 = (TTree*)fi->Get("ak3PFJetAnalyzer/t");
    // TTree *t0 = (TTree*)fi->Get("ak2PFJetAnalyzer/t");
    t0->SetBranchAddress("jtpt",jtpt);
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);
    t0->SetBranchAddress("rawpt",rawpt);
    t0->SetBranchAddress("nref",&nref);

    TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");
    t1->SetBranchAddress("vz",&vz);

    TTree *t2 = (TTree*)fi->Get("hltanalysis/HltTree");
    t2->SetBranchAddress("HLT_AK4CaloJet60_v1",&HLT_AKCJ60v1);
    t2->SetBranchAddress("HLT_AK4CaloJet80_v1",&HLT_AKCJ80v1);

    TTree *t3 = (TTree*)fi->Get("skimanalysis/HltTree");
    t3->SetBranchAddress("pprimaryVertexFilter",&ppVF);

    // for loop for events in trees t0, t1, t2, and t3
    for(unsigned int i=0; i<t0->GetEntries(); i++)
    {
        cout << "event " << i << " is being processed" << endl;

        // getting the entries
        t0->GetEntry(i);
        t1->GetEntry(i);
        t2->GetEntry(i);
        t3->GetEntry(i);

        // adding one to total events for every event
        a_+=1;

        // only events with |vz|<15 and all the triggers of interest are passed
        // if((TMath::Abs(vz)<15)&&(pPApVF==1)&&(HBHENFRR2L==1)&&(pBSF==1)&&(HLT_HIAKCJ80v1==1)&&(pthat>15)&&(refpt[0]>80)&&(TMath::Abs(jteta[0])<1.6)){
        // if((TMath::Abs(vz)<15)&&(ppVF==1)&&(HLT_AKCJ60v1==1)){
        if((TMath::Abs(vz)<15)&&(ppVF==1)&&(HLT_AKCJ80v1==1)){
            
            // adding one to passed events iff all the conditionals are true and the eta cut is NOT applied
            b_+=1;
    
            hvz->Fill(vz);
            
            // looping through all jets in each event
            for(unsigned int j=0; j<nref; j++){

                // filling histograms before eta cut
                hjteta_uc->Fill(jteta[j]);
                
                // eta cut
                // if(TMath::Abs(jteta[j])<1.6){
                // if((TMath::Abs(jteta[j])<1.6)&&(jtpt[j]>pth1d0[1])&&(jtpt[j]<pth1d0[2])){
                if((jtpt[j]>pth1d0[1])&&(jtpt[j]<pth1d0[2])){
                    // filling histograms that have variables with more than one value per event
                    hjtpt->Fill(jtpt[j]);
                    hjteta->Fill(jteta[j]);
                    hjtphi->Fill(jtphi[j]);
                    hrawpt->Fill(rawpt[j]);

                    // cout << "line 363 reached for event "<<i << endl;

                    // getting the corrected jtpt
                    vector<string> Files;
                    Files.push_back("ParallelMC_L2Relative_AK4PF_v0_12-21-2023.txt");
                    JetCorrector JEC(Files);

                    JEC.SetJetPT(rawpt[j]);
                    JEC.SetJetEta(jteta[j]);
                    JEC.SetJetPhi(jtphi[j]);  
                    Float_t jet_pt_corr = JEC.GetCorrectedPT();

                    jtcorrpt[j] = jet_pt_corr;

                    // pt balance
                    // making the iterator value for tag and probe jet, then adjusting value depending on eta values
                    int tagiter = 0;
                    int probeiter = 0;
                    // finding the A values first
                    // only considering leading and subleading jets
                    // tag is j = 1, probe is j = 0
                    // if((jteta[1]<tageta)&&(jteta[0]<5.2)){
                    //     tagiter = 1;
                    //     probeiter = 0;
                    // }
                    // tag is j = 0, probe is j = 1
                    if((jteta[0]<tageta)&&(jteta[1]<5.2)){
                        tagiter = 0;
                        probeiter = 1;
                    }
                    // as it is if both leading and subleading jet eta are < 1.3 then the leading jet is taken to be the tag jet

                    // if((jteta[1]<tageta)&&(jteta[0]>tageta)){
                        
                    //     // printing A value and saving it
                    //     cout << "A is " << (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]) << endl;
                    //     double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);

                    //     // pt slices A value filling
                    //     // each k is a different slice of pt
                    //     for(unsigned int k=0; k<ptslicenum; k++){
                    //         double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    //         if((ptavg>ptlow[k])&&(ptavg<pthigh[k])){
                    //             ptslicesA[k]->Fill(jteta[probeiter],Aval);
                    //             // ptsliceAtots[k] += Aval;
                    //             // ptsliceAnums[k] += 1;
                    //         }
                    //     }

                    //     // eta slices A value filling
                    //     // each k is a different slice of eta
                    //     for(unsigned int k=0; k<etaslicenum; k++){
                    //         double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    //         if((jteta[probeiter]>etalow[k])&&(jteta[probeiter]<etahigh[k])){
                    //             etaslicesA[k]->Fill(ptavg,Aval);
                    //             // etasliceAtots[k] += Aval;
                    //             // etasliceAnums[k] += 1;
                    //         }
                    //     }
                    // }

                    if(jteta[0]<tageta){

                        // printing A value and saving it
                        cout << "A is " << (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]) << endl;
                        double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);

                        // pt slices A value filling
                        // each k is a different slice of pt
                        for(unsigned int k=0; k<ptslicenum; k++){
                            double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                            if((ptavg>ptlow[k])&&(ptavg<pthigh[k])){
                                ptslicesA[k]->Fill(jteta[probeiter],Aval);
                            }
                        }

                        // eta slices A value filling
                        // each k is a different slice of eta
                        for(unsigned int k=0; k<etaslicenum; k++){
                            double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                            if((jteta[probeiter]>etalow[k])&&(jteta[probeiter]<etahigh[k])){
                                etaslicesA[k]->Fill(ptavg,Aval);
                            }
                        }
                    }

                    if((jtcorrpt[j]>pth1d0[1])&&(jtcorrpt[j]<pth1d0[2])){
                        hjtcorrpt->Fill(jtcorrpt[j]);}

                    // int t0 = a_;
                    // int t1 = j;
                    // if((t0%250==0)&&(t1%10 == 0)){
                        // cout << "for event " << a_ << endl<< endl;
                        // cout << "jtpt is " << jtpt[j] << endl;
                        // cout << "jteta is " << jteta[j] << endl;
                        // cout << "jtphi is " << jtphi[j] << endl;
                        // cout << "jtcorrpt is " << jtcorrpt[j] << endl;
                        // cout << "rawpt is " << rawpt[j] << endl;
                    // }
                }
            }   
        }
    }

    fi->Close();
    
    // getting y projection or slice of each eta bin for each pt slice
    for(unsigned int k=0; k<ptslicenum; k++){
        // example of ProjectionY() below
        //myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg]");
        ploth2d_1(ptslicesA[k], "eta", "Aval", Form("plots/ptslicesA[k]_for_k_is_%d.png",k), "");
        for(unsigned int l=0; l<etaslicenum; l++){
            etaslices_of_ptslicesA[k][l] = ptslicesA[k]->ProjectionY("",l,l,"");
            // getting y projection or slice of each pt bin for each eta slice
            ptslices_of_etaslicesA[l][k] = etaslicesA[l]->ProjectionY("",k,k,"");
            // saving these plots
            TString htitle1 = Form("plots/Aptslice%d_%f_%f__etabin%d_%f_%f.png",k,ptlow[k],pthigh[k],l,etalow[l],etahigh[l]);
            ploth1d_1_s(etaslices_of_ptslicesA[k][l], "A value",htitle1);
            TString htitle2 = Form("plots/Aetaslice%d_%f_%f__ptbin%d_%f_%f.png",l,etalow[l],etahigh[l],k,ptlow[k],pthigh[k]);
            ploth1d_1_s(ptslices_of_etaslicesA[l][k], "A value", htitle2);
        }
    } 

    Double_t ptslicesAavg[ptslicenum][etaslicenum];
    Double_t ptslicesAavgerr[ptslicenum][etaslicenum];
    Double_t ptslicesAx[etaslicenum];
    Double_t ptslicesAxerr[etaslicenum];

    Double_t etaslicesAavg[etaslicenum][ptslicenum];
    Double_t etaslicesAavgerr[etaslicenum][ptslicenum];
    Double_t etaslicesAx[ptslicenum];
    Double_t etaslicesAxerr[ptslicenum];

    for(unsigned int k=0; k<etaslicenum; k++){
        ptslicesAx[k] = (etahigh[k] + etalow[k])/2;
        ptslicesAxerr[k] = ptslicesAx[k] - etalow[k];
        ploth2d_1(etaslicesA[k], "pt", "Aval", Form("plots/etaslicesA[k]_for_k_is_%d.png",k), "");

        for(unsigned int l=0; l<ptslicenum; l++){
            // pt slices hists of A vs eta
            ptslicesAavg[l][k] = etaslices_of_ptslicesA[l][k]->GetMean();
            cout << "avg of etaslices_of_ptslicesA[l][k] for l = "<<l<<" and k = "<<k << " is "<<etaslices_of_ptslicesA[l][k]->GetMean() << endl;
            ptslicesAavgerr[l][k] = ptslices_of_etaslicesA[l][k]->GetMeanError();
            // cout << "avg error of etaslices_of_ptslicesA[l][k] for l = "<<l<<" and k = "<<k << " is "<<ptslicesAavgerr[l][k] << endl;
            // eta slices hists of A vs pt
            etaslicesAavg[k][l] = ptslices_of_etaslicesA[k][l]->GetMean();
            // cout << "avg of ptslices_of_etaslicesA[k][l] for l = "<<l<<" and k = "<<k << " is "<<etaslicesAavg[k][l] << endl;
            etaslicesAavgerr[k][l] = ptslices_of_etaslicesA[k][l]->GetMeanError();
            // cout << "avg error of ptslices_of_etaslicesA[k][l] for l = "<<l<<" and k = "<<k <<" is "<< etaslicesAavgerr[k][l] << endl;
        }
    }

    for(unsigned int k=0; k<ptslicenum; k++){
        etaslicesAx[k] = (pthigh[k] + ptlow[k])/2;
        etaslicesAxerr[k] = etaslicesAx[k] - ptlow[k];
    }

    cout << "line 527 reached"<< endl;

    TGraph *getaslicesAavg[ptslicenum]; 
    TGraph *gptslicesAavg[etaslicenum];

    cout << "line 532 reached"<< endl;

    for(unsigned int k=0; k<ptslicenum; k++){
        double ys[etaslicenum];
        double yserr[etaslicenum];
        cout << "line 537 reached for pt slice number " << k << endl;
        for(unsigned int l=0; l<etaslicenum; l++){
            cout << "line 540 reached for l = " << l << endl;
            ys[l] = etaslicesAavg[k][l];
            yserr[l] = etaslicesAavgerr[k][l];
        }
        cout << "line 542 reached for pt slice number " << k << endl;
        getaslicesAavg[k] = new TGraphErrors(ptslicenum,etaslicesAx,ys,etaslicesAxerr,yserr);
        cout << "line 544 reached for pt slice number " << k << endl;
        TString htitle = Form("plots/Aavg_ptslice%d_%f_%f.png",k,ptlow[k],pthigh[k]);
        // plotg1d_1_s(gptslicesAavg[k], "p_T [GeV/c]", "<A>", htitle);
        // plotg1d_1_s(gptslicesAavg[k], "p_T [GeV/c]", "<A>", Form("Aavg_ptslice%d_%f_%f",k,ptlow[k],pthigh[k]));
        // plotg1d_1_s(gptslicesAavg[k], "p_T [GeV/c]", "<A>", sprintf("Aavg_ptslice%d_%f_%f",k,ptlow[k],pthigh[k]));
    }

    cout << "line 545 reached"<< endl;

    for(unsigned int k=0; k<etaslicenum; k++){
        Double_t ys[ptslicenum];
        Double_t yserr[ptslicenum];
        for(unsigned int l=0; l<ptslicenum; l++){
            ys[l] = ptslicesAavg[k][l];
            yserr[l] = ptslicesAavgerr[k][l];
        }
        getaslicesAavg[k] = new TGraphErrors(etaslicenum,ptslicesAx,ys,ptslicesAxerr,yserr);
        getaslicesAavg[k]->SetMinimum(-0.1);
        getaslicesAavg[k]->SetMaximum(0.1);
        plotg1d_1_s(getaslicesAavg[k], "p_T [GeV/c]", "<A>", Form("plots/Aavg_etaslice%d_%f_%f.png",k,etalow[k],etahigh[k]));
    }

    cout << "line 563 reached"<< endl;

    // normalizing with the integral function
    normalizeh(hvz);
    normalizeh(hjtpt);
    normalizeh(hjtcorrpt);
    normalizeh(hrawpt);
    normalizeh(hjteta);
    normalizeh(hjtphi);
    normalizeh(hjteta_uc);

    // PLOTTING
    //
    // getting rid of boxes of stats in plots
    gStyle->SetOptStat(0);

    cout << "line 580 reached"<< endl;

    //ploth1d_1_s(TH1D *h, TString xtitle, TString htitle)
    // for(unsigned int k=0; k<ptslicenum; k++){
    //     //myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg]");
    //     for(unsigned int l=0; l<etaslicenum; l++){
    //         TString htitle = Form("plots/Aptslice%d_%f_%f__etabin%d_%f_%f.png",k,ptlow[k],pthigh[k],l,etalow[l],etahigh[l]);
    //         ploth1d_1_s(etaslices_of_ptslicesA[k][l], "A value",htitle);
    //     }
    // }
    // for(unsigned int k=0; k<etaslicenum; k++){
    //     //myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg]");
    //     for(unsigned int l=0; l<ptslicenum; l++){
    //         TString htitle = Form("plots/Aetaslice%d_%f_%f__ptbin%d_%f_%f.png",k,etalow[k],etahigh[k],l,ptlow[l],pthigh[l]);
    //         ploth1d_1_s(ptslices_of_etaslicesA[k][l], "A value", htitle);
    //     }
    // }


    // // 1D histograms
    // ploth1d_1(hrawpt, "rawpt", "Probability", "rawpt");
    // hjtpt->SetMaximum(0.2);
    // ploth1d_1(hjtpt, "jtpt", "Probability", "jtpt");
    // hjteta->SetMaximum(0.06);
    // ploth1d_1(hjteta, "jteta cut", "Probability", "jteta cut");
    // hjtphi->SetMaximum(0.02);
    // ploth1d_1(hjtphi, "jtphi", "Probability", "jtphi");
    // hjtcorrpt->SetMaximum(0.2);
    // ploth1d_1(hjtcorrpt, "jtcorrpt", "Probability", "jtcorrpt");
    // hjteta_uc->SetMaximum(0.06);
    // ploth1d_1(hjteta_uc, "jteta uncut", "Probability", "jteta uncut");

    // // Overlaid 1D histograms
    // ploth1d_2(hjteta_uc, "uncut jteta", hjteta, "cut jteta", "#eta", "Probability", "#eta^{jet} before and after #eta cut, individually normalized");
    // ploth1d_2(hrawpt, "rawpt", hjtcorrpt, "jtcorrpt", "pt", "Probability", "corrected and uncorrected jet pt");
    // ploth1d_2(hjtpt, "jtpt", hjtcorrpt, "jtcorrpt", "pt", "Probability", "corrected and uncorrected jet pt");

    // PRINT STATEMENTS
    // expressing passed/total as a percent
    // float p0 = b_/a_*100;

    // aavg = atot/anum;
    // Float_t rptrel = (1+aavg)/(1-aavg);
    // cout << "number of A values was " << anum << endl;
    // cout << "<A> = " << aavg << endl;
    // cout << "Rptrel = " << rptrel << endl;

    // printing out the info of interest
    // cout<<"number of total events: "<< a_ << endl;
    // cout<<"number of passed events before event cut: "<< b_ << endl;
    // cout<<"percent of passed events before eta cut: "<< p0 <<"%"<< endl;
    //
    // determining the final time of the script
    clock_t tf = clock();
    // printing out the time between ti and tf in seconds
    cout<<"this script took "<<double(tf-ti)/(double)CLOCKS_PER_SEC<<" seconds to run"<<endl<<endl;
}