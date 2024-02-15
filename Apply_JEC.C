// 

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

// PLOTTING FUNCTIONS

void save_h1d(TH1D *h, TString alg, TString xtitle, TString ytitle, TString hname){
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
    c->SaveAs("plots/"+hname+".png");
    // deleting c and h_c
    delete c;
    delete h_c;
}

void save_h2d(TH2D *h, TString alg, TString xtitle, TString ytitle, TString hname){
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
    l->AddEntry((TObject*)0, alg, "");
    l->AddEntry(h_c,hname,"pl");
    l->Draw("same");
    c->Write();
    c->SaveAs("plots/"+hname+".png");
    delete c;
    delete h_c;
}

void save_g(TGraph *h, TString alg, TString xtitle, TString ytitle, TString hname){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    h->Draw("AP");
    h->SetTitle("");
    h->SetName(hname);
    cout << " passed set name in save graph function " << endl;
    h->SetMaximum(0.1);
    h->SetMinimum(-0.2);
    Double_t w1 = 600;
    Double_t h1 = 600;
    c->SetCanvasSize(w1,h1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    h->SetMarkerStyle(20);
    cout << " passed gstyle in save graph function " << endl;
    TLine *line = new TLine(2,0.4,2,0.6);
    line->Draw("same");
    cout << " passed tline stuff in save graph function " << endl;
    // h->GetXaxis()->CenterTitle();
    // h->GetHistogram()->GetYaxis()->SetTitle("new Y axis title");
    // h->GetYaxis()->SetTitle(ytitle);
    // h->GetXaxis()->SetTitle(xtitle);
    // h->GetYaxis()->CenterTitle(true);
    // h->GetXaxis()->CenterTitle(true);
    cout << " passed axis titles in save graph function " << endl;
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetHeader(alg,"C");
    // l->AddEntry((TObject*)0, alg, "");
    l->AddEntry(h,hname,"pl");
    l->Draw("same");
    c->Write();
    c->SaveAs("plots/"+hname+".png");
    delete c;
    cout << " finished save graph function " <<endl<<endl;
}

void save_g_1(TH1D *h_, TString option1, TGraph *h, TString alg, TString xtitle, TString ytitle, TString hname){
    // canvas stuff
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    Double_t w1 = 600;
    c->SetCanvasSize(w1,w1); 
    // hist to make stuff look like what I want
    h_->GetYaxis()->SetTitle(ytitle);
    h_->GetXaxis()->SetTitle(xtitle);
    h_->GetYaxis()->CenterTitle(true);
    h_->GetXaxis()->CenterTitle(true);
    h_->SetMaximum(0.1);
    h_->SetMinimum(-0.2);
    h_->SetTitle("");
    h_->Draw("e1p");
    // plotting the graph of interest 
    h->SetTitle("");
    h->SetMarkerStyle(20);
    h->SetMaximum(0.1);
    h->SetMinimum(-0.2);
    h->Draw("P SAME");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    if(option1=="pt"){
        TLine *line = new TLine(80,0,500,0);
        line->Draw("same");
        line->SetLineStyle(2);
        }
    if(option1=="eta"){
        TLine *line = new TLine(-5.2,0,5.2,0);
        line->Draw("same");
        line->SetLineStyle(2);
        }
    // legend stuff
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetHeader(alg,"C");
    l->AddEntry(h,hname,"pl");
    l->SetTextSize(50);
    l->Draw("same");
    gStyle->SetOptStat(0);
    c->Write();
    c->SaveAs("plots/"+hname+".png");
    // gStyle->SetOptStat(1);
    delete c;
}

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// the script all runs in this function
void Apply_JEC()
{
    // determining the inital time of the script
    Timer timer = Timer();
    timer.Start();

    // INITIALIZING HISTOGRAMS
    timer.StartSplit("Histogram Initialization");

    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();

    // creating some binning parameters
    //////////////////////////////////////////////////////////////////
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    double vzh1d0[3] = {40,-20,20};
    double pth1d0[3] = {100,80,500};
    double phih1d0[3] = {100,4,4};
    double etah1d0[3] = {50,-5.2,5.2};
    double etah1d1[3] = {25,-1.7,1.7};
    double ah1d0[3] = {100,-1,1};
    
    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 7;
    // the low and high pt values for each pt slice
    double ptlow[ptslicenum] = {80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {100,120,140,180,220,300,500};
    
    // eta slices
    // number of eta slices
    const Int_t etaslicenum = 8;
    // the low and high eta values for each eta slice
    double etalow[etaslicenum] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9};
    double etahigh[etaslicenum] = {-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};
    //////////////////////////////////////////////////////////////////

    // Making some hists for pt balance studies
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // the actual A and R values
    TH2D *ptslicesA[ptslicenum];
    TH2D *etaslicesA[etaslicenum];
    
    // eta bin slices for the pt slice hists of A vs eta_probe
    TH1D *etaslices_of_ptslicesA[ptslicenum][etaslicenum];
    // pt bin slices for the eta slice hists of A vs pt_avg
    TH1D *ptslices_of_etaslicesA[etaslicenum][ptslicenum];

    // <A> as well as R vs pt and eta for eta and pt slices respectively
    TGraph *getaslicesAavg[ptslicenum];
    TGraph *getaslicesR[ptslicenum]; 
    TGraph *gptslicesAavg[etaslicenum]; 
    TGraph *gptslicesR[etaslicenum];

    // further initializing the histograms
    // looping over pt slices
    for(unsigned int q=0; q<ptslicenum; q++){
        // A vs eta for each pt slice with title having the high and low pt for the slice
        TString chtitle0 = Form("A_ptslice_%.0f_%.0f",ptlow[q],pthigh[q]);
        ptslicesA[q] = new TH2D(chtitle0,chtitle0,etah1d0[0],etah1d0[1],etah1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        for(unsigned int r=0; r<etaslicenum; r++){
            // avoiding looping through the etaslicenum separately by only doing it on the first q value
            if(q==0){
                // A vs pt for each eta slice with title having the high and low eta for the slice
                TString ahtitle0 = Form("A_etaslice_%.0f_%.0f",etalow[r]*10,etahigh[r]*10);
                etaslicesA[r] = new TH2D(ahtitle0,ahtitle0,pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
            }
            // intializing hists thatre are the bins of the slices
            TString dhtitle0 = Form("Aptslice%d_%.0f_%.0f__etabin%d_%.0f_%.0f",q,ptlow[q],pthigh[q],r,etalow[r]*10,etahigh[r]*10);
            etaslices_of_ptslicesA[q][r] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
            TString dhtitle1 = Form("Aetaslice%d_%.0f_%.0f__ptbin%d_%.0f_%.0f",r,etalow[r]*10,etahigh[r]*10,q,ptlow[q],pthigh[q]);
            ptslices_of_etaslicesA[r][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // initializing histograms for general parameters
    ////////////////////////////////////////////////////////////////////////////
    TH1D *hvz = new TH1D("vz","vz",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hjtpt = new TH1D("hjtpt","hjtpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjtcorrpt = new TH1D("hjtcorrpt","hjtcorrpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","hjteta",etah1d1[0],etah1d1[1],etah1d1[2]);
    TH1D *hjtphi = new TH1D("hjtphi","hjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *hrawpt = new TH1D("hrawpt","hrawpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta_uc = new TH1D("hjteta_uc","hjteta_uc",etah1d0[0],etah1d0[1],etah1d0[2]);
    ////////////////////////////////////////////////////////////////////////////

    // initializing hists for canvases later
    //////////////////////////////////////////////////////////////////////
    TH1D *hpt = new TH1D("hpt","hpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *heta = new TH1D("heta","heta",etah1d0[0],etah1d0[1],etah1d0[2]);
    //////////////////////////////////////////////////////////////////////

    // INITIALIZING PARAMETERS
    timer.StartSplit("Parameter Initialization");

    // number of jets in the barrel with a pt above 80 GeV
    int pt80barreljetnum = 0;

    // a_ and b_ are the number of total events and accepted events respectively
    Float_t a_=0, b_=0, c_=0;

    // pointing fi to the file holding the jet info of interest
    // TFile *fi = TFile::Open("HiForestMiniAOD_10k.root","read");
    // TFile *fi = TFile::Open("HP0_1_25_2024.root","read");
    TFile *fi = TFile::Open("HP_All_2_4_2024.root","read");
    // TFile *fi = TFile::Open("nbarnett@lxplus.cern.ch:eos/user/n/nbarnett/PPRefHardProbes1/crab_foresting_run373710_HP1_Uncorrected_1-25-2024_0/240125_170921/0000/HP_All_2_4_2024.root","read");

    // declaring variables
    ////////////////////////////////////////////////////////////////////////////
    // triggers of interest
    int ppVF, HLT_AKCJ60v1, HLT_AKCJ80v1;
    // vertex position
    Float_t vz;
    // whatever ptcut I decide
    Float_t ptcut = 60;
    // a big number to make my arrays such that they aren't too small
    const Int_t nm = 200000;
    // uncorrected jet pt
    Float_t rawpt[nm];
    // incorrect corrected jet pt
    Float_t jtpt[nm];
    // correct corrected jet pt
    Float_t jtcorrpt[nm];
    // jet phi and pseudorapadity
    Float_t jtphi[nm];
    Float_t jteta[nm];
    // number of jets in event
    Int_t nref;
    // only making tag jets in the barrel, which means a pseudorapadity < tageta
    Float_t tageta = 1.3;
    ////////////////////////////////////////////////////////////////////////////

    // pt balance study arrays
    ////////////////////////////////////////////////////////////////////////////
    // pt slices
    // <A> for each pt slice
    Double_t ptslicesAavg[ptslicenum][etaslicenum];
    // the uncertainty or error in <A> for each pt slice
    Double_t ptslicesAavgerr[ptslicenum][etaslicenum];
    // eta values on x axis
    Double_t ptslicesAx[etaslicenum];
    // making tgraphs with information so the error in x is half the bin width
    Double_t ptslicesAxerr[etaslicenum];
    // cooresponding values for eta slices 
    Double_t etaslicesAavg[etaslicenum][ptslicenum];
    Double_t etaslicesAavgerr[etaslicenum][ptslicenum];
    Double_t etaslicesAx[ptslicenum];
    Double_t etaslicesAxerr[ptslicenum];
    ////////////////////////////////////////////////////////////////////////////

    timer.StartSplit("Setting TTree addresses");

    // getting the TTrees from the file

    // toggling which clustering alg info to look
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");
    // TTree *t0 = (TTree*)fi->Get("ak3PFJetAnalyzer/t");
    // TTree *t0 = (TTree*)fi->Get("ak2PFJetAnalyzer/t");

    // general parameters
    //////////////////////////////////////
    t0->SetBranchAddress("jtpt",jtpt);
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);
    t0->SetBranchAddress("rawpt",rawpt);
    t0->SetBranchAddress("nref",&nref);
    //////////////////////////////////////

    // event cut info in the following trees
    ///////////////////////////////////////////////////////////
    TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");
    t1->SetBranchAddress("vz",&vz);

    TTree *t2 = (TTree*)fi->Get("hltanalysis/HltTree");
    t2->SetBranchAddress("HLT_AK4CaloJet60_v1",&HLT_AKCJ60v1);
    t2->SetBranchAddress("HLT_AK4CaloJet80_v1",&HLT_AKCJ80v1);

    TTree *t3 = (TTree*)fi->Get("skimanalysis/HltTree");
    t3->SetBranchAddress("pprimaryVertexFilter",&ppVF);
    ///////////////////////////////////////////////////////////

    // for loop going over events in the trees
    // for(unsigned int i=0; i<t0->GetEntries(); i++){
    for(unsigned int i=0; i<100; i++){

        // timer 0
        int i_0 = i;
        if(i_0%50==0){timer.StartSplit(Form("event_%d_until line 416",i));}

        cout<< "event " << i << " is being processed" << endl;

        // getting the entries in every event
        ////////////////
        t0->GetEntry(i);
        t1->GetEntry(i);
        t2->GetEntry(i);
        t3->GetEntry(i);
        ////////////////

        // adding one to total events for every event
        a_+=1;

        // EVENT CUT
        // only events with |vz|<15 and all the triggers of interest are passed
        if((TMath::Abs(vz)<15)&&(ppVF==1)&&(HLT_AKCJ60v1==1)){
            
            // adding one to passed events iff all the conditionals are true and the eta cut is NOT applied
            b_+=1;

            // filling the vertex position hist
            hvz->Fill(vz);
            
            // looping through all jets in each event
            for(unsigned int j=0; j<nref; j++){

                // filling histograms before eta and pt cut
                ////////////////////////////////////////////////
                hjteta_uc->Fill(jteta[j]);
                
                if((rawpt[j]>pth1d0[1])&&(rawpt[j]<pth1d0[2])){
                    hrawpt->Fill(rawpt[j]);
                }
                ////////////////////////////////////////////////

                // getting the corrected jtpt
                //////////////////////////////////////////////////////////////////
                vector<string> Files;
                Files.push_back("ParallelMC_L2Relative_AK4PF_v0_12-21-2023.txt");
                JetCorrector JEC(Files);

                JEC.SetJetPT(rawpt[j]);
                JEC.SetJetEta(jteta[j]);
                JEC.SetJetPhi(jtphi[j]);  
                Float_t jet_pt_corr = JEC.GetCorrectedPT();
                //////////////////////////////////////////////////////////////////

                // saving the corrected jet pt
                jtcorrpt[j] = jet_pt_corr;
                
                // timer 1
                int i_1 = i;
                if((i_1%50==0)&&(j==0)){timer.StartSplit(Form("event_%d_jet[%d]_until_line_434",i,j));}

                // Filling some hists
                //////////////////////////////////////////////////////
                if((jtcorrpt[j]>pth1d0[1])&&(jtcorrpt[j]<pth1d0[2])){
                    hjtcorrpt->Fill(jtcorrpt[j]);
                }

                if((jtcorrpt[j]>80)&&(jteta[j]<1.3)){
                    pt80barreljetnum+=1;
                }

                // only look at pt balance studies if jtcorrpt > ptcut
                if(jtcorrpt[j]>ptcut){
                    // filling histograms that have variables with more than one value per event
                    hjtpt->Fill(jtpt[j]);
                    hjteta->Fill(jteta[j]);
                    hjtphi->Fill(jtphi[j]);
                }
                //////////////////////////////////////////////////////

                // Some print statements
                //////////////////////////////////////////////////////
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
                //////////////////////////////////////////////////////
            }

            // timer 2
            int i_2 = i;
            if(i_2%50==0){timer.StartSplit(Form("pt balance stuff for event %d",i));}

            // PT BALANCE
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // making the iterator values for tag, probe, and third leading jet (if there is one)
            // tag jet must be leading, and probe jet must be subleading
            // then adjusting the iterator values depending on pt order of jets in the event
            int tagiter = 0;
            int probeiter = 1;

            // finding the leading jet, which is the possible tag jet
            // looping through the jets in the event
            for(unsigned int j=0; j<nref; j++){
                // only when the jet pt is highest will it replace the current iterator
                // in the case this is never true the original leading jet would still be the leading jet and iter would be 0
                if(jtcorrpt[j]>jtcorrpt[tagiter]){
                    tagiter=j;
                    // if now the leading is the original subleading jet
                    // then the subleading jet is changed to another iter before being found
                    // this will be true iff tagiter = 1 or the leading jet is the original subleading jet index
                    if(probeiter==tagiter){
                        probeiter=0;
                    }
                }
            }

            // finding the subleading jet, which is the possible probe jet
            // looping through the jets in the event
            for(unsigned int j=0; j<nref; j++){
                // only when the jet pt is larger than current subleading jet and smaller than leading jet
                // in the case this is never true the original subleading, or possibly leading, jet would be the subleading jet and iter would be 1, or possibly 0
                if((jtcorrpt[j]<jtcorrpt[tagiter])&&(jtcorrpt[j]>jtcorrpt[probeiter])){
                    probeiter=j;
                }
            }

            // finding the third leading jet, iff there are at least three jets
            if(nref>2){
                int thirditer = 2;
                // start assuming the third leading jet is the third leading jet still, unless either the leading jet or subleading jet is already
                // only looking at the first three jets
                for(unsigned int q=0; q<3; q++){
                    // one of the first three jets isn't the leading or subleading jet and we initialize the third jet to be that one
                    if((q!=probeiter)&&(q!=tagiter)){
                        thirditer = q;
                    }
                }
                // looping through the jets in the event
                for(unsigned int j=2; j<nref; j++){
                    // determining if another jet between index 3 and nref-1 exists with higher pt than the current third jet, but only if it is less pt than and isn't the probe or tag jet
                    if((jtcorrpt[j]>jtcorrpt[thirditer])&&(jtcorrpt[j]<jtcorrpt[probeiter])&&(jtcorrpt[j]<jtcorrpt[tagiter])&&(j!=probeiter)&&(j!=tagiter)){
                        thirditer = j;
                    }
                }              
            }
        
            // finding the A values iff the leading jet has eta < 1.3 and subleading jet passes the pt cut and has eta < 5.2 
            if((TMath::Abs(jteta[tagiter]<tageta))&&(TMath::Abs(jteta[probeiter]<5.2))&&(jtcorrpt[probeiter]>ptcut)){
                
                // printing A value and saving it
                cout << "A is " << (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]) << endl;
                double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);

                // pt slices A value filling
                // each k is a different slice of pt
                for(unsigned int k=0; k<ptslicenum; k++){
                    // average momentum between the probe and tag jet 
                    // these are sliced originally to get A vs eta for different pt slices
                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    // if the pt avg is within a certain slice range then add it to the pt slice hists
                    if((ptavg>ptlow[k])&&(ptavg<pthigh[k])){
                        // ptslicesA are A vs eta_probe hists for different pt ranges or slices
                        ptslicesA[k]->Fill(jteta[probeiter],Aval);
                    }
                }

                // eta slices A value filling
                // each k is a different slice of eta
                for(unsigned int k=0; k<etaslicenum; k++){
                    // ptavg is the x axis in one desired type of plot
                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    // if the eta is within a certain slice range then add it to the eta slice hists
                    if((jteta[probeiter]>etalow[k])&&(jteta[probeiter]<etahigh[k])){
                        // etaslicesA are A vs pt_avg hists for different eta ranges
                        etaslicesA[k]->Fill(ptavg,Aval);
                    }
                }
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }

    timer.StartSplit("after event loop");

    // closing the file I'm getting the information from
    fi->Close();
    cout <<endl<< "the number of jets in the barrel with pT above 80 GeV was " << pt80barreljetnum << endl<<endl;

    // normalizing some hists with the integral function
    ///////////////////////
    normalizeh(hvz);
    normalizeh(hjtpt);
    normalizeh(hjtcorrpt);
    normalizeh(hrawpt);
    normalizeh(hjteta);
    normalizeh(hjtphi); 
    ///////////////////////

    cout<<"line 514"<<endl<<endl;

    // making a new file to store all the histograms of interest in
    // TFile f("pt_balance.root","recreate");
    TFile *f = new TFile("pt_balance.root", "recreate");
    f->cd();

    // writing the base variable hists to this new file
    cout<<"line 521"<<endl<<endl;
    ///////////////////
    TString alghere = "ak4pf";
    cout<<"line 524"<<endl<<endl;
    save_h1d(hvz, alghere, "vz", "Probability", "hvz");
    cout<<"line 526"<<endl<<endl;
    save_h1d(hjtcorrpt, alghere, "p_{T}^{corrected jet}", "Probability", "hjtcorrpt");
    save_h1d(hjteta, alghere, "η", "Probability", "hjeta");
    save_h1d(hjtphi, alghere, "φ", "Probability", "hjtphi");
    save_h1d(hrawpt, alghere, "p_{T}", "Probability", "hrawpt");
    
    ///////////////////

    // looping throught the pt slices
    for(unsigned int k=0; k<ptslicenum; k++){

        // saving the A vs eta plots for each pt slice
        TString bhtitle0 = Form("ptslicesA_ptbin%d_%.0f_%.0f",k,ptlow[k],pthigh[k]);
        save_h2d(ptslicesA[k], alghere, "η", "A", bhtitle0);
        cout<<"line 537"<<endl<<endl;

        // making the x axis points for the eta slices be the center of each pt bin
        etaslicesAx[k] = (pthigh[k] + ptlow[k])/2;
        // making the x axis points error for the eta slices be half the width of each pt bin
        etaslicesAxerr[k] = etaslicesAx[k] - ptlow[k];

        // looping through the eta slices
        for(unsigned int l=0; l<etaslicenum; l++){

            // conditional below acts like an separated l loop 
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(k==0){
                // saving eta slice hists of A vs pt
                TString bhtitle1 = Form("etaslicesA_etabin%d_%.0f_%.0f",l,etalow[l]*10,etahigh[l]*10);
                save_h2d(etaslicesA[l], alghere, "p_{T}^{avg}", "A", bhtitle1);
                // pt axis for tgraphs
                ptslicesAx[l] = (etahigh[l] + etalow[l])/2;
                ptslicesAxerr[l] = ptslicesAx[l] - etalow[l];
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // projecting the slices onto the the A axis for each x axis bin for each slice
            //////////////////////////////////////////////////////////////////////////////////
            // example of ProjectionY(): myhist->ProjectionY(" ",firstxbin,lastxbin,"[cutg]");

            // getting y projection or slice of each eta bin for each pt slice
            etaslices_of_ptslicesA[k][l] = ptslicesA[k]->ProjectionY("",l,l,"");

            // getting y projection or slice of each pt bin for each eta slice
            ptslices_of_etaslicesA[l][k] = etaslicesA[l]->ProjectionY("",k,k,"");
            //////////////////////////////////////////////////////////////////////////////////
            
            // saving the projection hists
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // eta bins of pt slices
            // making the title
            TString htitle1 = Form("Aptslice%d_%.0f_%.0f_etabin%d_%.0f_%.0f",k,ptlow[k],pthigh[k],l,etalow[l]*10,etahigh[l]*10);
            // saving the A value projections for eta bins of the pt slices
            save_h1d(etaslices_of_ptslicesA[k][l], alghere, "A", "probability", htitle1);
            cout<<"line 577"<<endl<<endl;

            // pt bins of eta slices
            TString htitle2 = Form("Aetaslice%d_%.0f_%.0f_ptbin%d_%.0f_%.0f",l,etalow[l]*10,etahigh[l]*10,k,ptlow[k],pthigh[k]);
            // saving the A value projections for eta bins of the pt slices
            save_h1d(ptslices_of_etaslicesA[l][k], alghere, "A", "Probability", htitle2);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // finding stuff for tgraphs
            /////////////////////////////////////////////////////////////////////////
            // pt slices hists of <A> vs eta
            Double_t ptaavg = (etaslices_of_ptslicesA[k][l])->GetMean();
            Double_t ptaavgerr = (etaslices_of_ptslicesA[k][l])->GetMeanError();
            ptslicesAavg[k][l] = ptaavg;
            ptslicesAavgerr[k][l] = ptaavgerr;

            // eta slices hists of <A> vs pt
            Double_t etaaavg = (ptslices_of_etaslicesA[l][k])->GetMean();
            Double_t etaaavgerr = (ptslices_of_etaslicesA[l][k])->GetMeanError();
            etaslicesAavg[l][k] = etaaavg;
            etaslicesAavgerr[l][k] = etaaavgerr;


            /////////////////////////////////////////////////////////////////////////
            
            // print statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cout<<"<A> for η bin "<<etalow[l]<<" to "<<etahigh[l]<<" of pt slice "<<ptlow[k]<<" to "<<pthigh[k]<<" is "<<ptslicesAavg[k][l]<<"±"<<ptslicesAavgerr[k][l]<<endl;
            cout<<"<A> for pt bin "<<ptlow[l]<<" to "<<pthigh[l]<<" of η slice "<<etalow[k]<<" to "<<etahigh[k]<<" is "<<etaslicesAavg[k][l]<<"±"<<etaslicesAavgerr[k][l]<<endl<<endl;
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }

    // actually finding the averages of the A values for the tgraphs
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    // <A> vs pt for each eta slice
    for(unsigned int k=0; k<ptslicenum; k++){

        // making 1D arrays for the TGraph function inputs
        Double_t ys[etaslicenum];
        Double_t yserr[etaslicenum];

        // Assigning all the values to the 1D arrays
        for(unsigned int l=0; l<etaslicenum; l++){
            ys[l] = etaslicesAavg[k][l];
            yserr[l] = etaslicesAavgerr[k][l];
        }

        // making the tgraphs
        getaslicesAavg[k] = new TGraphErrors(ptslicenum,etaslicesAx,ys,etaslicesAxerr,yserr);
        cout<<endl<<"line 623"<<endl<<endl;
        // making the title for the traph
        TString htitle3 = Form("<A>_eta_%.0f_%.0f",etalow[k]*10,etahigh[k]*10);
        // saving the <A> vs pT
        cout<<"line 627"<<endl<<endl;
        // gptslicesAavg[k]->SetTitle(htitle3);
        // cout<<"line 629"<<endl<<endl;
        // gptslicesAavg[k]->SetName(htitle3);
        // cout<<"line 631"<<endl<<endl;
        // gptslicesAavg[k]->Write();
        cout<<"line 633"<<endl<<endl;
        // save_g(getaslicesAavg[k], alghere, "p_T [GeV/c]", "<A>", htitle3);
        save_g_1(hpt, "pt", getaslicesAavg[k], alghere, "p_T [GeV/c]", "<A>", htitle3);
        cout<<"line 635"<<endl<<endl;
    }

    // <A> vs eta for each pt slice
    // same as previous loop but for switched eta and pt
    for(unsigned int k=0; k<etaslicenum; k++){
        Double_t ys[ptslicenum];
        Double_t yserr[ptslicenum];
        for(unsigned int l=0; l<ptslicenum; l++){
            ys[l] = ptslicesAavg[k][l];
            yserr[l] = ptslicesAavgerr[k][l];
        }
        gptslicesAavg[k] = new TGraphErrors(etaslicenum,ptslicesAx,ys,ptslicesAxerr,yserr);
        // getaslicesAavg[k]->SetMinimum(-0.1);
        // getaslicesAavg[k]->SetMaximum(0.1);
        TString hname = Form("<A>_pt_%.0f_%.0f",ptlow[k],pthigh[k]);
        // getaslicesAavg[k]->SetTitle(hname);
        // getaslicesAavg[k]->SetName(hname);
        // getaslicesAavg[k]->Write();
        // save_g(gptslicesAavg[k], alghere, "η", "<A>", hname);
        save_g_1(heta, "eta", gptslicesAavg[k], alghere, "η", "<A>", hname);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////

    // closing the root file storing all of the histograms
    f->Close();

    // FINAL PRINT STATEMENTS

    // expressing passed/total as a percent
    // float p0 = b_/a_*100;

    // printing out the event cut amount
    // cout<<"number of total events: "<< a_ << endl;
    // cout<<"number of passed events before event cut: "<< b_ << endl;
    // cout<<"percent of passed events before eta cut: "<< p0 <<"%"<< endl;
    //
    // determining the final time of the script
    timer.Stop();
    cout<<endl;
    timer.Report();
}