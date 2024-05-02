// IMPORTS //

// included header files
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

// // separate header files
// #include "Timer.h"
// #include "JetUncertainty.h"
// #include "JetCorrector.h"

// PLOTTING FUNCTIONS //

// saving a TGraph to an output root file
void save_g(TGraph *h, TString hname){
    h->SetName(hname);
    h->Write();
}

// saving a TH1D to an output root file
void save_h1d(TH1D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

// saving a TH2D to an output root file
void save_h2d(TH2D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

// saving a TH1D to an output root file & a png in a folder
void save_h1d_1(TH1D *h, TString xtitle, TString ytitle, TString hname){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
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
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h_c,hname,"pl");
    l->Draw("same");
    c->Write();
    h->Write();
    c->SaveAs("plots_5_2/"+hname+".png");
    delete c;
    delete h_c;
}

// saving a TH2D to an output root file & a png in a folder
void save_h2d_1(TH2D *h, TString xtitle, TString ytitle, TString hname){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    TH2D *h_c = (TH2D*)h->Clone(hname);
    h_c->Draw("COLZ");
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
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h_c,hname,"pl");
    l->Draw("same");
    c->Write();
    h->Write();
    c->SaveAs("plots_5_2/"+hname+".png");
    delete c;
    delete h_c;
}

// OTHER FUNCTIONS //

// normalizing a TH1D by integrating and scaling by the inverse of the integration value
void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// a way to visualize the jet and its subjets in eta phi space
// for _mom arrays [0] is pt, [1] is eta, [2] is phi
void visualize_subjets(Double_t jet_mom[3], Double_t subjet_mom[3][30], Double_t subjet_conesize, Double_t jet_conesize, TH2D *h, TString hname){    
    
    // making the canvas
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    c->Range(-5.2,-4,5.2,4);
    // c->GetXaxis()->SetTitle("#eta");
    // c->GetYaxis()->SetTitle("#phi");
    // c->GetYaxis()->CenterTitle(true);
    // c->GetXaxis()->CenterTitle(true);
    c->SetGridx();
    c->SetGridy();

    Double_t jet_mom_c[3] = {0};
    Double_t subjet_mom_c[3][30] = {0};
    for(unsigned int i=0; i<3; i++){
        jet_mom_c[i] = jet_mom[i];
        for(unsigned int s=0; s<30; s++){
            subjet_mom_c[i][s] = subjet_mom[i][s];
        }
    }
    
    // // drawing the background hist
    // TH1D *h_c = (TH1D*)h->Clone(hname);
    // h_c->Draw("e1p");
    // h_c->SetMarkerStyle(20);
    // h_c->SetMarkerColor(kBlack);
    // h_c->SetLineColor(kBlack);
    // h_c->SetTitle("");
    // h_c->SetName(hname);
    // h_c->GetYaxis()->SetTitle("#phi");
    // h_c->GetXaxis()->SetTitle("#eta");
    // h_c->GetYaxis()->CenterTitle(true);
    // h_c->GetXaxis()->CenterTitle(true);

    // making the circle for the jet
    auto elj1 = new TEllipse(jet_mom[1], jet_mom[2], jet_conesize, jet_conesize);
    
    // drawing the circle for the jet
    elj1->Draw("");

    // making the circle for the subjets
    if(subjet_mom[1][0]!=-999){
        auto els0 = new TEllipse(subjet_mom[1][0], subjet_mom[2][0], subjet_conesize, subjet_conesize);
        els0->Draw("");
        cout<<"1"<<endl;
        cout<<"pt is "<<subjet_mom[0][0] <<endl;
        cout<<"eta is "<<subjet_mom[1][0] <<endl;
        cout<<"phi is "<<subjet_mom[2][0] <<endl;
    }
    if(subjet_mom[1][1]!=-999){
        auto els1 = new TEllipse(subjet_mom[1][1], subjet_mom[2][1], subjet_conesize, subjet_conesize);
        els1->Draw("");
        cout<<"2"<<endl;
        cout<<"pt is "<<subjet_mom[0][1] <<endl;
        cout<<"eta is "<<subjet_mom[1][1] <<endl;
        cout<<"phi is "<<subjet_mom[2][1] <<endl;
    }
    if(subjet_mom[1][2]!=-999){
        auto els2 = new TEllipse(subjet_mom[1][2], subjet_mom[2][2], subjet_conesize, subjet_conesize);
        els2->Draw("");
        cout<<"3"<<endl;
        cout<<"pt is "<<subjet_mom[0][2] <<endl;
        cout<<"eta is "<<subjet_mom[1][2] <<endl;
        cout<<"phi is "<<subjet_mom[2][2] <<endl;
    }
    if(subjet_mom[1][3]!=-999){
        auto els3 = new TEllipse(subjet_mom[1][3], subjet_mom[2][3], subjet_conesize, subjet_conesize);
        els3->Draw("");
        cout<<"4"<<endl;
    }
    if(subjet_mom[1][4]!=-999){
        auto els4 = new TEllipse(subjet_mom[1][4], subjet_mom[2][4], subjet_conesize, subjet_conesize);
        els4->Draw("");
        cout<<"5"<<endl;
    }
    if(subjet_mom[1][5]!=-999){
        auto els5 = new TEllipse(subjet_mom[1][5], subjet_mom[2][5], subjet_conesize, subjet_conesize);
        els5->Draw("");
        cout<<"6"<<endl;
    }
    if(subjet_mom[1][6]!=-999){
        auto els6 = new TEllipse(subjet_mom[1][6], subjet_mom[2][6], subjet_conesize, subjet_conesize);
        els6->Draw("");
        cout<<"7"<<endl;
    }
    if(subjet_mom[1][7]!=-999){
        auto els7 = new TEllipse(subjet_mom[1][7], subjet_mom[2][7], subjet_conesize, subjet_conesize);
        els7->Draw("");
        cout<<"8"<<endl;
    }
    if(subjet_mom[1][8]!=-999){
        auto els8 = new TEllipse(subjet_mom[1][8], subjet_mom[2][8], subjet_conesize, subjet_conesize);
        els8->Draw("");
        cout<<"9"<<endl;
    }
    if(subjet_mom[1][9]!=-999){
        auto els9 = new TEllipse(subjet_mom[1][9], subjet_mom[2][9], subjet_conesize, subjet_conesize);
        els9->Draw("");
        cout<<"10"<<endl;
    }
    if(subjet_mom[1][10]!=-999){
        auto els10 = new TEllipse(subjet_mom[1][10], subjet_mom[2][10], subjet_conesize, subjet_conesize);
        els10->Draw("");
        cout<<"11"<<endl;
    }
    if(subjet_mom[1][11]!=-999){
        auto els11 = new TEllipse(subjet_mom[1][11], subjet_mom[2][11], subjet_conesize, subjet_conesize);
        els11->Draw("");
        cout<<"12"<<endl;
    }
    if(subjet_mom[1][12]!=-999){
        auto els12 = new TEllipse(subjet_mom[1][12], subjet_mom[2][12], subjet_conesize, subjet_conesize);
        els12->Draw("");
        cout<<"13"<<endl;
    }
    if(subjet_mom[1][13]!=-999){
        auto els13 = new TEllipse(subjet_mom[1][13], subjet_mom[2][13], subjet_conesize, subjet_conesize);
        els13->Draw("");
        cout<<"14"<<endl;
        // if(){}
    }
    if(subjet_mom[1][14]!=-999){
        auto els14 = new TEllipse(subjet_mom[1][14], subjet_mom[2][14], subjet_conesize, subjet_conesize);
        els14->Draw("");
        cout<<"15"<<endl;
    }
    if(subjet_mom[1][15]!=-999){
        auto els15 = new TEllipse(subjet_mom[1][15], subjet_mom[2][15], subjet_conesize, subjet_conesize);
        els15->Draw("");
        cout<<"16"<<endl;
    }
    if(subjet_mom[1][16]!=-999){
        auto els16 = new TEllipse(subjet_mom[1][16], subjet_mom[2][16], subjet_conesize, subjet_conesize);
        els16->Draw("");
        cout<<"17"<<endl;
    }
    if(subjet_mom[1][17]!=-999){
        auto els17 = new TEllipse(subjet_mom[1][17], subjet_mom[2][17], subjet_conesize, subjet_conesize);
        els17->Draw("");
        cout<<"18"<<endl;
    }
    if(subjet_mom[1][18]!=-999){
        auto els18 = new TEllipse(subjet_mom[1][18], subjet_mom[2][18], subjet_conesize, subjet_conesize);
        els18->Draw("");
        cout<<"19"<<endl;
    }
    if(subjet_mom[1][19]!=-999){
        auto els19 = new TEllipse(subjet_mom[1][19], subjet_mom[2][19], subjet_conesize, subjet_conesize);
        els19->Draw("");
        cout<<"20"<<endl;
    }
    if(subjet_mom[1][20]!=-999){
        auto els20 = new TEllipse(subjet_mom[1][20], subjet_mom[2][20], subjet_conesize, subjet_conesize);
        els20->Draw("");
        cout<<"21"<<endl;
    }
    if(subjet_mom[1][21]!=-999){
        auto els21 = new TEllipse(subjet_mom[1][21], subjet_mom[2][21], subjet_conesize, subjet_conesize);
        els21->Draw("");
    }
    if(subjet_mom[1][22]!=-999){
        auto els22 = new TEllipse(subjet_mom[1][22], subjet_mom[2][22], subjet_conesize, subjet_conesize);
        els22->Draw("");
    }
    if(subjet_mom[1][23]!=-999){
        auto els23 = new TEllipse(subjet_mom[1][23], subjet_mom[2][23], subjet_conesize, subjet_conesize);
        els23->Draw("");
    }
    if(subjet_mom[1][24]!=-999){
        auto els24 = new TEllipse(subjet_mom[1][24], subjet_mom[2][24], subjet_conesize, subjet_conesize);
        els24->Draw("");
    }
    if(subjet_mom[1][25]!=-999){
        auto els25 = new TEllipse(subjet_mom[1][25], subjet_mom[2][25], subjet_conesize, subjet_conesize);
        els25->Draw("");
    }
    if(subjet_mom[1][26]!=-999){
        auto els26 = new TEllipse(subjet_mom[1][26], subjet_mom[2][26], subjet_conesize, subjet_conesize);
        els26->Draw("");
    }
    if(subjet_mom[1][27]!=-999){
        auto els27 = new TEllipse(subjet_mom[1][27], subjet_mom[2][27], subjet_conesize, subjet_conesize);
        els27->Draw("");
    }
    if(subjet_mom[1][28]!=-999){
        auto els28 = new TEllipse(subjet_mom[1][28], subjet_mom[2][28], subjet_conesize, subjet_conesize);
        els28->Draw("");
    }
    if(subjet_mom[1][29]!=-999){
        auto els29 = new TEllipse(subjet_mom[1][29], subjet_mom[2][29], subjet_conesize, subjet_conesize);
        els29->Draw(""); 
    }


    // SUBJET LOOP
    for(unsigned int s=0; s<30; s++){

        // only going through subjets that aren't flagged as bad
        if(subjet_mom[0][s]!=-999){

            // making the circle for the subjet
            auto el2 = new TEllipse(subjet_mom[1][s], subjet_mom[2][s], subjet_conesize, subjet_conesize);

            // drawing the circle for the subjet
            el2->Draw();
        }
    }

    // writing out the plot
    c->Write();

    // deleting canvas to avoid abusing memory
    delete c;

}

// MAIN CODE FUNCTION //

void subjet_analysis(){
    
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();

    // // getting rid of auto legends
    // gStyle->SetOptStat(0);

    // MAKING HIST BINS //

    // vz
    double vzh1d0[3] = {40,-20,20};

    // numbers
    double numh1d0[3] = {30,0,30};

    // momenta
    double pth1d0[3] = {100,80,500};
    double etah1d0[3] = {50,-5.2,5.2};
    double phih1d0[3] = {100,4,4};

    // MAKING SPECIFIC PT HIST BINNING //

    // number of pt slices
    const Int_t ptslicenum = 10;

    // the low and high pt values for each pt slice
    double ptlow[ptslicenum] = {15,25,50,80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {25,50,80,100,120,140,180,220,300,500};

    // INITIALIZING EVENT HISTS //

    // vz
    TH1D *hvz = new TH1D("vz","vz",vzh1d0[0],vzh1d0[1],vzh1d0[2]);

    // numbers of jets
    TH1D *hngen = new TH1D("hngen","hngen",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hnref = new TH1D("hnref","hnref",numh1d0[0],numh1d0[1],numh1d0[2]);
    
    // INITIALIZING JET HISTS //
    
    // reco jet
    TH1D *hjtpt = new TH1D("hjtpt","hjtpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","hjteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hjtphi = new TH1D("hjtphi","hjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // ref jet
    TH1D *hrefpt = new TH1D("hrefpt","hrefpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hrefeta = new TH1D("hrefeta","hrefeta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hrefphi = new TH1D("hrefphi","hrefphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // gen jet
    TH1D *hgenpt = new TH1D("hgenpt","hgenpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hgeneta = new TH1D("hgeneta","hgeneta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hgenphi = new TH1D("hgenphi","hgenphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    
    // INITIALIZING SUBJET HISTS //

    // 1D HISTOGRAMS    

    // cone size 0.1 
    // reco
    TH1D *hsubjetnum1 = new TH1D("hsubjetnum1","hsubjetnum1",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hsubjetpt1 = new TH1D("hsubjetpt1","hsubjetpt1",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hsubjeteta1 = new TH1D("hsubjeteta1","hsubjeteta1",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hsubjetphi1 = new TH1D("hsubjetphi1","hsubjetphi1",phih1d0[0],phih1d0[1],phih1d0[2]);
    // gen
    TH1D *hsubjetgennum1 = new TH1D("hsubjetgennum1","hsubjetgennum1",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hsubjetgenpt1 = new TH1D("hsubjetgenpt1","hsubjetgenpt1",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hsubjetgeneta1 = new TH1D("hsubjetgeneta1","hsubjetgeneta1",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hsubjetgenphi1 = new TH1D("hsubjetgenphi1","hsubjetgenphi1",phih1d0[0],phih1d0[1],phih1d0[2]);

    // 2D HISTOGRAMS

    // eta - phi plane
    TH2D *hetaphi = new TH2D("hetaphi","hetaphi",etah1d0[0],etah1d0[1],etah1d0[2],phih1d0[0],phih1d0[1],phih1d0[2]);

    // cone size 0.1
    // number vs pt for reco and gen
    TH2D *hsubjetnum1_jtpt = new TH2D("hsubjetnum1_jtpt","hsubjetnum1_jtpt",pth1d0[0],pth1d0[1],pth1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    TH2D *hsubjetgennum1_genpt = new TH2D("hsubjetgennum1_genpt","hsubjetgennum1_genpt",pth1d0[0],pth1d0[1],pth1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    // // number of gen subjets vs reco subjets for various jet pt slices
    // TH2D *hsubjetgennum1_subjetnum1_ptcut5 = new TH2D("hsubjetgennum1_subjetnum1_ptcut5","hsubjetgennum1_subjetnum1_ptcut5",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    // TH2D *hsubjetgennum1_subjetnum1_ptcut5_60_80 = new TH2D("hsubjetgennum1_subjetnum1_ptcut5_60_80","hsubjetgennum1_subjetnum1_ptcut5_60_80",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    // TH2D *hsubjetgennum1_subjetnum1_ptcut5_80_100 = new TH2D("hsubjetgennum1_subjetnum1_ptcut5_80_100","hsubjetgennum1_subjetnum1_ptcut5_80_100",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    // TH2D *hsubjetgennum1_subjetnum1_ptcut5_100_120 = new TH2D("hsubjetgennum1_subjetnum1_ptcut5_100_120","hsubjetgennum1_subjetnum1_ptcut5_100_120",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    // TH2D *hsubjetgennum1_subjetnum1_ptcut5_120_ = new TH2D("hsubjetgennum1_subjetnum1_ptcut5_120_","hsubjetgennum1_subjetnum1_ptcut5_120_",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    
    // DECLARING VARIABLES //
    
    // EVENT VARIABLES

    // vertex position
    Float_t vz;

    // event weight
    Float_t w;

    // event flag of interest
    int ppVF;
    
    // number of jets in event
    Int_t nref;
    Int_t ngen;
    
    // a big number to make my jet and subjet arrays such that they aren't too small
    const Int_t MAXJETS = 500;
    const Int_t SUBJETA = 15000;

    // JET VARIABLES
    
    // reco jet momenta
    Float_t jtpt[MAXJETS];
    Float_t jtphi[MAXJETS];
    Float_t jteta[MAXJETS];
    
    // gen jet momenta
    Float_t genpt[MAXJETS];
    Float_t geneta[MAXJETS];
    Float_t genphi[MAXJETS];
    
    // ref jet momenta
    Float_t refpt[MAXJETS];
    Float_t refeta[MAXJETS];
    Float_t refphi[MAXJETS];

    // SUBJET VARIABLES
    
    // cone size 0.1
    // reco
    Float_t subjetnum1[MAXJETS];
    Float_t subjetpt1[SUBJETA];
    Float_t subjeteta1[SUBJETA];
    Float_t subjetphi1[SUBJETA];
    // gen
    Float_t subjetgennum1[MAXJETS];
    Float_t subjetgenpt1[SUBJETA];
    Float_t subjetgeneta1[SUBJETA];
    Float_t subjetgenphi1[SUBJETA];

    // COUNTING VARIABLES

    // number of subjets
    Double_t subjetcount = 0;

    // number of subjets with pt > 1 TeV
    Double_t subjetcount_hugept = 0;

    // INPUT //

    // INPUT FILE
    TFile *fi = TFile::Open("HiForestAOD_5_2_2024_0.root","read");
    
    // GETTING TTREES
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");
    TTree *t1 = (TTree*)fi->Get("skimanalysis/HltTree");
    TTree *t2 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");

    // SETTING BRANCH STATUSES //

    // EVENT BRANCHES

    // ak4PFJetAnalyzer tree
    t0->SetBranchStatus("*",0);
    t0->SetBranchStatus("vz",1);
    t0->SetBranchStatus("nref",1);
    t0->SetBranchStatus("ngen",1);

    // skimanalysis/HltTree tree
    t1->SetBranchStatus("*",0);
    t1->SetBranchStatus("pPAprimaryVertexFilter",1);

    // hiEvtAnalyzer/HiTree tree
    t2->SetBranchStatus("*",0);
    t2->SetBranchStatus("weight",1);
    
    // JET BRANCHES

    // reco
    t0->SetBranchStatus("jtpt",1);
    t0->SetBranchStatus("jteta",1);
    t0->SetBranchStatus("jtphi",1);

    // ref
    t0->SetBranchStatus("refpt",1);
    t0->SetBranchStatus("refeta",1);
    t0->SetBranchStatus("refphi",1);

    // gen
    t0->SetBranchStatus("genpt",1);
    t0->SetBranchStatus("geneta",1);
    t0->SetBranchStatus("genphi",1);

    // SUBJET BRANCHES

    // cone size 0.1
    // reco
    t0->SetBranchStatus("subjetnum1",1);
    t0->SetBranchStatus("subjetpt1",1);
    t0->SetBranchStatus("subjeteta1",1);
    t0->SetBranchStatus("subjetphi1",1);
    // gen
    t0->SetBranchStatus("subjetgennum1",1);
    t0->SetBranchStatus("subjetgenpt1",1);
    t0->SetBranchStatus("subjetgeneta1",1);
    t0->SetBranchStatus("subjetgenphi1",1);
    
    // SETTING BRANCH ADDRESSES //

    // EVENT BRANCHES
    t0->SetBranchAddress("nref",&nref);
    t0->SetBranchAddress("ngen",&ngen);
    t0->SetBranchAddress("vz",&vz);
    t2->SetBranchAddress("weight",&w);
    t1->SetBranchAddress("pPAprimaryVertexFilter",&ppVF);

    // JET BRANCHES

    // reco
    t0->SetBranchAddress("jtpt",jtpt);
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);

    // ref
    t0->SetBranchAddress("refpt",refpt);
    t0->SetBranchAddress("refeta",refeta);
    t0->SetBranchAddress("refphi",refphi);

    // gen
    t0->SetBranchAddress("genpt",genpt);
    t0->SetBranchAddress("geneta",geneta);
    t0->SetBranchAddress("genphi",genphi);

    // SUBJET BRANCHES

    // cone size 0.1
    // reco
    t0->SetBranchAddress("subjetnum1",subjetnum1);
    t0->SetBranchAddress("subjetpt1",subjetpt1);
    t0->SetBranchAddress("subjeteta1",subjeteta1);
    t0->SetBranchAddress("subjetphi1",subjetphi1);
    // gen
    t0->SetBranchAddress("subjetgennum1",subjetgennum1);
    t0->SetBranchAddress("subjetgenpt1",subjetgenpt1);
    t0->SetBranchAddress("subjetgeneta1",subjetgeneta1);
    t0->SetBranchAddress("subjetgenphi1",subjetgenphi1);

    // OUTPUT FILE //

    // making output file name
    TString output = "subjet_analysis_5_2_2024_0.root";

    // making and pointing to output file
    TFile *f1 = new TFile(output, "recreate");
    
    // EVENT PROCESSING //

    // EVENT LOOP
    for(unsigned int i=0; i<t0->GetEntries(); i++){

        // // printing event number
        // cout<< "event " << i << " is being processed" << endl;

        // only getting the ttree with the event flag needed to be passed first
        t1->GetEntry(i);

        // getting the other ttree info for each entry iff the event flag is passed
        if(ppVF==1){
            t0->GetEntry(i);

            // only events with |vz|<15 are passed
            if(TMath::Abs(vz)<15){

                // only getting the weights for the completely passed events
                t2->GetEntry(i);

                // filling event hists
                hvz->Fill(vz, w);
                hngen->Fill(ngen, w);
                hnref->Fill(nref, w);
                
                // // printing ngen and nref
                // cout << "ngen is "<<ngen<< " and nref is "<< nref << " for event " <<i<<endl;

                // SUBJET VISUALIZATION TEST

                // // testing the subjet visualization function on only the subleading jet of event 75
                // if(i==50){

                // initializing the inputs for the function, setting all entries in them to 0 initially
                Double_t jet_mom[3] = {0};
                Double_t subjet_mom[3][30] = {0};

                // getting the jet momenta of interest
                jet_mom[0] = jtpt[1];
                jet_mom[1] = jteta[1];
                jet_mom[2] = jtphi[1];

                // subjet loop in jet of interest
                for(unsigned int s=0; s<30; s++){

                    // the subjet index for sub jet s for the subleading jet in this event
                    int subjetindex = 30+s;

                    // getting the subjet momenta of interest
                    subjet_mom[0][s] = subjetpt1[subjetindex];
                    subjet_mom[1][s] = subjeteta1[subjetindex];
                    subjet_mom[2][s] = subjetphi1[subjetindex];
                }

                // staging the output file so Write() saves to the output file
                f1->cd();

                // runnning the visualization function over the jet and corresponding subjets of interest 
                visualize_subjets(jet_mom, subjet_mom, 0.1, 0.4, hetaphi, "subjet_visualization");

                // }

                // staging the input file again so I can still Get() everything needed
                fi->cd();

                // RECO JET LOOP
                for(unsigned int j=0; j<nref; j++){

                    // only looking at reco jet info, and associated subjet info, if the reco jet pt > 60 GeV
                    if(jtpt[j]<60){continue;}

                    // FILLING 1D HISTS

                    // number of subjets in reco jets
                    hsubjetnum1->Fill(subjetnum1[j], w);

                    // reco jet momenta
                    hjtpt->Fill(jtpt[j], w);
                    hjteta->Fill(jteta[j], w);
                    hjtphi->Fill(jtphi[j], w);

                    // ref jet momenta
                    hrefpt->Fill(refpt[j], w);
                    hrefeta->Fill(refeta[j], w);
                    hrefphi->Fill(refphi[j], w);

                    // FILLING 2D HISTS

                    // reco jet pt vs number of subjets
                    hsubjetnum1_jtpt->Fill(jtpt[j], subjetnum1[j], w);

                    // RECO SUBJET LOOP
                    for(unsigned int s=0; s<subjetnum1[j]; s++){

                        // the subjet index for sub jet s for jet j for event i
                        int subjetindex = 30*j+s;

                        // filling reco subjet momenta hists
                        hsubjetpt1->Fill(subjetpt1[subjetindex], w);
                        hsubjeteta1->Fill(subjeteta1[subjetindex], w);
                        hsubjetphi1->Fill(subjetphi1[subjetindex], w);
                    }
                }

                // GEN JET LOOP
                for(unsigned int g=0; g<ngen; g++){

                    // only looking at gen jet info, and associated subjet info, if the gen jet pt > 60 GeV
                    if(genpt[g]<60){continue;}

                    // FILLING 1D HISTS 

                    // gen jet momenta
                    hgenpt->Fill(genpt[g], w);
                    hgeneta->Fill(geneta[g], w);
                    hgenphi->Fill(genphi[g], w);

                    // number of subjets in gen jets
                    hsubjetgennum1->Fill(subjetgennum1[g], w);

                    // FILLING 2D HISTS

                    // gen jet pt vs number of subjets
                    hsubjetgennum1_genpt->Fill(genpt[g], subjetgennum1[g], w);

                    // I think this is an oopsy in my logic, I should look at refpt
                    if(g>nref){continue;}

                    // INITIALIZING EVENT SPECIFIC SUBJET VARIABLES

                    // cone size 0.1
                    // number of subjets in gen jet with pt > 5 GeV
                    Double_t subjetgennum1_ptcut5 = 0;
                    // number of subjets in reco jet with pt > 5 GeV
                    Double_t subjetnum1_ptcut5 = 0;

                    // GEN SUBJET LOOP
                    for(unsigned int s=0; s<subjetgennum1[g]; s++){

                        // subjet index of interest
                        int subjetindex = 30*g+s;

                        // just to make sure I don't go over subjet entries that aren't subjets
                        if((subjetgenpt1[subjetindex]==-999)||(subjetpt1[subjetindex]==-999)){continue;}

                        // cone size 0.1
                        //gen 
                        if(subjetgenpt1[subjetindex]*w>5){subjetgennum1_ptcut5+=1;}
                        // reco
                        if(subjetpt1[subjetindex]*w>5){subjetnum1_ptcut5+=1;}

                        // filling gen subjet momenta hists
                        hsubjetgenpt1->Fill(subjetgenpt1[subjetindex], w);
                        hsubjetgeneta1->Fill(subjetgeneta1[subjetindex], w);
                        hsubjetgenphi1->Fill(subjetgenphi1[subjetindex], w);


                        // adding to count of number of subjets
                        subjetcount +=1;

                        // adding to count of number of subjet with pt > 1 TeV
                        if(subjetgenpt1[subjetindex]>1000){subjetcount_hugept +=1;}

                        // printing any subjet pt that is > 1 TeV
                        // if(subjetgenpt1[subjetindex]>1000){cout<<"gen subjet pt for cone size 0.1 is "<<subjetgenpt1[subjetindex]<<" for event "<<i<<" jet number "<<g<<" and subjet number "<<s<<endl;}
                        // if(subjetpt1[subjetindex]>1000){cout << "reco subjet pt for cone size 0.1 is "<< subjetpt1[subjetindex] << " for event "<<i<<" jet number "<<g<<" and subjet number "<<s<<endl;}
                    }
                    
                    // // 2d plots of subjet num vs subjetgen num
                    // // cone size 0.1
                    // if((subjetgennum1_ptcut5!=0)&&(subjetnum1_ptcut5!=0)){
                    //     hsubjetgennum1_subjetnum1_ptcut5->Fill(subjetgennum1_ptcut5, subjetnum1_ptcut5);
                    // }
                    // if((subjetgennum1_ptcut5!=0)&&(subjetnum1_ptcut5!=0)&&(jtpt[g]>60)&&(jtpt[g]<80)){
                    //     hsubjetgennum1_subjetnum1_ptcut5_60_80->Fill(subjetgennum1_ptcut5, subjetnum1_ptcut5);
                    // }
                    // if((subjetgennum1_ptcut5!=0)&&(subjetnum1_ptcut5!=0)&&(jtpt[g]>80)&&(jtpt[g]<100)){
                    //     hsubjetgennum1_subjetnum1_ptcut5_80_100->Fill(subjetgennum1_ptcut5, subjetnum1_ptcut5);
                    // }
                    // if((subjetgennum1_ptcut5!=0)&&(subjetnum1_ptcut5!=0)&&(jtpt[g]>100)&&(jtpt[g]<120)){
                    //     hsubjetgennum1_subjetnum1_ptcut5_100_120->Fill(subjetgennum1_ptcut5, subjetnum1_ptcut5);
                    // }
                    // if((subjetgennum1_ptcut5!=0)&&(subjetnum1_ptcut5!=0)&&(jtpt[g]>120)){
                    //     hsubjetgennum1_subjetnum1_ptcut5_120_->Fill(subjetgennum1_ptcut5, subjetnum1_ptcut5);
                    // }
                }
            }
        }
    }

    // closing the file I'm getting the information from
    fi->Close();

    // opening the output file to save into it
    f1->cd();

    // NORMALIZING HISTOGRAMS //

    // EVENT HISTS
    normalizeh(hvz);

    // RECO JET HISTS
    normalizeh(hjtpt);
    normalizeh(hjteta);
    normalizeh(hjtphi);

    // REF JET HISTS
    normalizeh(hrefpt);
    normalizeh(hrefeta);
    normalizeh(hrefphi);

    // GEN JET HISTS
    normalizeh(hgenpt);
    normalizeh(hgeneta);
    normalizeh(hgenphi);

    // SUBJET HISTS

    // cone size 0.1
    // reco
    normalizeh(hsubjetpt1);
    normalizeh(hsubjeteta1);
    normalizeh(hsubjetphi1);
    // gen
    normalizeh(hsubjetgenpt1);
    normalizeh(hsubjetgeneta1);
    normalizeh(hsubjetgenphi1);

    // WRITING HISTOGRAMS //

    // EVENT HISTS
    hvz->Write();
    hngen->Write();
    hnref->Write();

    // JET HISTS

    // reco
    hjtpt->Write();
    hjteta->Write();
    hjtphi->Write();

    // ref
    hrefpt->Write();
    hrefeta->Write();
    hrefphi->Write();

    // gen
    hgenpt->Write();
    hgeneta->Write();
    hgenphi->Write();

    // SUBJET 1D HISTS
    
    // cone size 0.1 
    // reco
    hsubjetnum1->Write();
    save_h1d_1(hsubjetpt1, "pt", "Probability", "subjetpt1");
    save_h1d_1(hsubjeteta1, "#eta", "Probability", "hsubjeteta1");
    save_h1d_1(hsubjetphi1, "#phi", "Probability", "hsubjetphi1");
    // gen
    hsubjetgennum1->Write();
    save_h1d_1(hsubjetgenpt1, "pt", "Probability", "subjetgenpt1");
    save_h1d_1(hsubjetgeneta1, "eta", "Probability", "hsubjetgeneta1");
    save_h1d_1(hsubjetgenphi1, "phi", "Probability", "hsubjetgenphi1");

    // SUBJET 2D HISTS

    // // number of gen vs reco subjets in various jet pt regions
    // save_h2d_1(hsubjetgennum1_subjetnum1_ptcut5, "reco subjetnum1", "gen subjetnum1", "subjetgennum1_subjetnum1_pt5cut");
    // save_h2d_1(hsubjetgennum1_subjetnum1_ptcut5_60_80, "reco subjetnum1", "gen subjetnum1", "subjetgennum1_subjetnum1_pt5cut_60_80");
    // save_h2d_1(hsubjetgennum1_subjetnum1_ptcut5_80_100, "reco subjetnum1", "gen subjetnum1", "subjetgennum1_subjetnum1_pt5cut_80_100");
    // save_h2d_1(hsubjetgennum1_subjetnum1_ptcut5_100_120, "reco subjetnum1", "gen subjetnum1", "subjetgennum1_subjetnum1_pt5cut_100_120");
    // save_h2d_1(hsubjetgennum1_subjetnum1_ptcut5_120_, "reco subjetnum1", "gen subjetnum1", "subjetgennum1_subjetnum1_pt5cut_120_");
    
    // reco and gen jet pt vs subjet number
    save_h2d_1(hsubjetnum1_jtpt, "reco subjetnum1", "reco jet pt", "subjetnum1_jtpt");
    save_h2d_1(hsubjetgennum1_genpt, "gen subjetnum1", "gen jet pt", "subjetgennum1_genpt");

    // PRINTING COUNTS OF INTEREST
    cout <<endl<< "the number of subjets looked at was " << subjetcount <<endl;
    cout << "the number of subjets with pt > 1 TeV looked at was " << subjetcount_hugept <<endl<<endl;

}
