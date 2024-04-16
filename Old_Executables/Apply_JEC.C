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

void save_h2d_1(TH2D *h, TString alg, TString xtitle, TString ytitle, TString hname){
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
    h->Write();
    c->SaveAs("plots/"+hname+".png");
    delete c;
    delete h_c;
}

void save_g_1(TH1D *h_, TString option1,TString option2, TGraph *h, TString alg, TString xtitle, TString ytitle, TString hname){
    // canvas stuff
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    // hist to make stuff look like what I want
    h_->GetYaxis()->SetTitle(ytitle);
    h_->GetXaxis()->SetTitle(xtitle);
    h_->GetYaxis()->CenterTitle(true);
    h_->GetXaxis()->CenterTitle(true);
    h_->SetTitle("");
    h_->Draw("e1p");
    // plotting the graph of interest 
    h->SetMarkerStyle(20);
    h->Draw("P SAME");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    if(option2=="A"){
        h_->SetMaximum(0.1);
        h_->SetMinimum(-0.1);
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
    }
    if(option2=="R"){
        h_->SetMaximum(1.1);
        h_->SetMinimum(0.9);
        if(option1=="pt"){
            TLine *line = new TLine(80,1,500,1);
            line->Draw("same");
            line->SetLineStyle(2);
        }
        if(option1=="eta"){
            TLine *line = new TLine(-5.2,1,5.2,1);
            line->Draw("same");
            line->SetLineStyle(2);
        }
    }

    // legend stuff
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetHeader(alg,"C");
    l->AddEntry(h,hname,"p");
    l->SetTextSize(50);
    l->Draw();
    c->Write();
    c->SaveAs("plots/"+hname+".png");
    delete c;
}

void save_g(TGraph *h, TString hname, TFile *f){
    h->SetName(hname);
    // f->cd();
    h->Write();
    // f->Close();
}

void save_h1d(TH1D *h, TString hname, TFile *f){
    h->SetName(hname);
    // f->cd();
    h->Write();
    // f->Close();
}

void save_h2d(TH2D *h, TString hname, TFile *f){
    h->SetName(hname);
    // f->cd();
    h->Write();
    // f->Close();
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
    gStyle->SetOptStat(0);

    //cout<<"189 line"<<endl;

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

    // for generating random numbers
    TRandom2 *rand = new TRandom2(1);
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

    //cout<<"240 line"<<endl;

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
                if(etalow[r]<0){
                    TString ahtitle0 = Form("A_etaslice_%.0f_%.0f_n",TMath::Abs(etalow[r]*10),TMath::Abs(etahigh[r]*10));
                    etaslicesA[r] = new TH2D(ahtitle0,ahtitle0,pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
                }
                if(etalow[r]>0||etalow[r]==0){
                    TString ahtitle0 = Form("A_etaslice_%.0f_%.0f",etalow[r]*10,etahigh[r]*10);
                    etaslicesA[r] = new TH2D(ahtitle0,ahtitle0,pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
                }
            }
            // intializing hists that are the bins of the slices
            if(etalow[r]<0){
                TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[q],pthigh[q],TMath::Abs(etalow[r]*10),TMath::Abs(etahigh[r]*10));
                TString dhtitle1 = Form("A_etaslice_%.0f_%.0f_n__ptbin_%.0f_%.0f",TMath::Abs(etalow[r]*10),TMath::Abs(etahigh[r]*10),ptlow[q],pthigh[q]);
                etaslices_of_ptslicesA[q][r] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                ptslices_of_etaslicesA[r][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
            }
            if(etalow[r]>0||etalow[r]==0){
                TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[q],pthigh[q],etalow[r]*10,etahigh[r]*10);
                TString dhtitle1 = Form("A_etaslice_%.0f_%.0f__ptbin_%.0f_%.0f",etalow[r]*10,etahigh[r]*10,ptlow[q],pthigh[q]);
                etaslices_of_ptslicesA[q][r] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                ptslices_of_etaslicesA[r][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
            }
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
    TFile *fi = TFile::Open("HP0_1_25_2024.root","read");
    // TFile *fi = TFile::Open("forested_MC_ppRef_2-6-2023_1s.root","read");
    // TFile *fi = TFile::Open("/eos/user/n/nbarnett/PPRefHardProbes1/crab_foresting_run373710_HP1_Uncorrected_1-25-2024_0/240125_170921/0000/HP_All_2_4_2024.root","read");

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
    // R for each pt slice
    Double_t ptslicesR[ptslicenum][etaslicenum];
    // the uncertainty or error in R for each pt slice
    Double_t ptslicesRerr[ptslicenum][etaslicenum];
    // eta values on x axis
    Double_t ptslicesAx[etaslicenum];
    // making tgraphs with information so the error in x is half the bin width
    Double_t ptslicesAxerr[etaslicenum];
    // cooresponding values for eta slices 
    Double_t etaslicesAavg[etaslicenum][ptslicenum];
    Double_t etaslicesAavgerr[etaslicenum][ptslicenum];
    Double_t etaslicesR[etaslicenum][ptslicenum];
    Double_t etaslicesRerr[etaslicenum][ptslicenum];
    Double_t etaslicesAx[ptslicenum];
    Double_t etaslicesAxerr[ptslicenum];
    ////////////////////////////////////////////////////////////////////////////

    timer.StartSplit("Setting TTree addresses");

    //cout<<"362 line"<<endl;
    // getting the TTrees from the file

    // toggling which clustering alg info to look
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");
    //cout<<"367 line"<<endl;
    // TTree *t0 = (TTree*)fi->Get("ak3PFJetAnalyzer/t");
    // TTree *t0 = (TTree*)fi->Get("ak2PFJetAnalyzer/t");

    // general parameters
    //////////////////////////////////////
    t0->SetBranchAddress("jtpt",jtpt);
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);
    t0->SetBranchAddress("rawpt",rawpt);
    t0->SetBranchAddress("nref",&nref);
    //cout<<"378 line"<<endl;
    //////////////////////////////////////

    // event cut info in the following trees
    ///////////////////////////////////////////////////////////
    TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");
    //cout<<"384 line"<<endl;
    t1->SetBranchAddress("vz",&vz);

    //cout<<"387 line"<<endl;
    TTree *t2 = (TTree*)fi->Get("hltanalysis/HltTree");
    t2->SetBranchAddress("HLT_AK4CaloJet60_v1",&HLT_AKCJ60v1);
    t2->SetBranchAddress("HLT_AK4CaloJet80_v1",&HLT_AKCJ80v1);

    // TTree *t3 = (TTree*)fi->Get("skimanalysis/HltTree");
    // t3->SetBranchAddress("pprimaryVertexFilter",&ppVF);
    ///////////////////////////////////////////////////////////


    //cout<<"397 line"<<endl;
    // for loop going over events in the trees
    // for(unsigned int i=0; i<t0->GetEntries(); i++){
    for(unsigned int i=0; i<100; i++){

        // timer 0
        int i_0 = i;
        if(i_0%50==0){timer.StartSplit(Form("event_%d_until line 416",i));}

        cout<< "event " << i << " is being processed" << endl;

        // getting the entries in every event
        ////////////////
        t2->GetEntry(i);
        // t3->GetEntry(i);
        ////////////////

        // adding one to total events for every event
        a_+=1;

        // EVENT CUT
        // only events with |vz|<15 and all the triggers of interest are passed
        // if((ppVF==1)&&(HLT_AKCJ60v1==1)){
        if(HLT_AKCJ60v1==1){
            t1->GetEntry(i);
            if(TMath::Abs(vz)<15){
                t0->GetEntry(i);
                

                
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
                        // //cout << "for event " << a_ << endl<< endl;
                        // //cout << "jtpt is " << jtpt[j] << endl;
                        // //cout << "jteta is " << jteta[j] << endl;
                        // //cout << "jtphi is " << jtphi[j] << endl;
                        // //cout << "jtcorrpt is " << jtcorrpt[j] << endl;
                        // //cout << "rawpt is " << rawpt[j] << endl;
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
                int leaditer = 0;
                int subleaditer = 1;
                int tagiter = 0;
                int probeiter = 1;

                // finding the leading jet, which is the possible tag jet
                // looping through the jets in the event
                for(unsigned int j=0; j<nref; j++){
                    // only when the jet pt is highest will it replace the current iterator
                    // in the case this is never true the original leading jet would still be the leading jet and iter would be 0
                    if(jtcorrpt[j]>jtcorrpt[leaditer]){
                        leaditer=j;
                        // if now the leading is the original subleading jet
                        // then the subleading jet is changed to another iter before being found
                        // this will be true iff tagiter = 1 or the leading jet is the original subleading jet index
                        if(subleaditer==leaditer){
                            subleaditer=0;
                        }
                    }
                }

                // finding the subleading jet, which is the possible probe jet
                // looping through the jets in the event
                for(unsigned int j=0; j<nref; j++){
                    // only when the jet pt is larger than current subleading jet and smaller than leading jet
                    // in the case this is never true the original subleading, or possibly leading, jet would be the subleading jet and iter would be 1, or possibly 0
                    if((jtcorrpt[j]<jtcorrpt[leaditer])&&(jtcorrpt[j]>jtcorrpt[subleaditer])){
                        subleaditer=j;
                    }
                }

                // pt balance study for the case nref < 3
                // defining tag and probe iters based on leading and subleading jet iters
                // tag and probe must be either leading or subleading jet
                // if the leading jet is in the barrel, tag it
                if(TMath::Abs(jteta[leaditer]<1.3)){
                    tagiter = leaditer;
                    probeiter = subleaditer;
                }
                // if the subleading jet is in the barrel, tag it
                if(TMath::Abs(jteta[subleaditer]<1.3)){
                    tagiter = subleaditer;
                    probeiter = leaditer;
                }
                // in the case both jets are in the barrel we make a random number 
                // if the random number is even or odd the tag jet is the leading or subleading jet respectively
                if((TMath::Abs(jteta[leaditer]<1.3))&&(TMath::Abs(jteta[subleaditer]<1.3))){
                    int checkval1 = rand->Integer(101);
                    if((checkval1%2==0)&&(nref<3)){
                        tagiter = leaditer;
                        probeiter = subleaditer;
                        //cout<<"the random number is "<< checkval1<<" and even, leading jet is tagged"<<endl;
                    }
                    if((checkval1%2!=0)&&(nref<3)){
                        tagiter = subleaditer;
                        probeiter = leaditer;
                        //cout<<"the random number is "<< checkval1<<" and odd, subleading jet is tagged"<<endl;
                    }
                }

                // finding the third leading jet, iff there are at least three jets
                if(nref>2){
                    int thirditer = 2;
                    // start assuming the third leading jet is the third leading jet still, unless either the leading jet or subleading jet is already
                    // only looking at the first three jets
                    for(unsigned int q=0; q<3; q++){
                        // one of the first three jets isn't the leading or subleading jet and we initialize the third jet to be that one
                        if((q!=subleaditer)&&(q!=leaditer)){
                            thirditer = q;
                        }
                    }
                    // looping through the jets in the event
                    for(unsigned int j=2; j<nref; j++){
                        // determining if another jet between index 3 and nref-1 exists with higher pt than the current third jet, but only if it is less pt than and isn't the probe or tag jet
                        if((jtcorrpt[j]>jtcorrpt[thirditer])&&(jtcorrpt[j]<jtcorrpt[subleaditer])&&(jtcorrpt[j]<jtcorrpt[leaditer])&&(j!=subleaditer)&&(j!=leaditer)){
                            thirditer = j;
                        }
                    }
                    // still working if there is a third jet
                    // doing the whole pt balance study in the case there is a third jet
                    if((jtcorrpt[subleaditer]+jtcorrpt[leaditer])>0.4*jtcorrpt[thirditer]){
                        if(TMath::Abs(jteta[leaditer]<1.3)){
                            tagiter = leaditer;
                            probeiter = subleaditer;
                        }
                        if(TMath::Abs(jteta[subleaditer]<1.3)){
                            tagiter = subleaditer;
                            probeiter = leaditer;
                        }
                        if((TMath::Abs(jteta[leaditer]<1.3))&&(TMath::Abs(jteta[subleaditer]<1.3))){
                            int checkval = rand->Integer(101);
                            if(checkval%2==0){
                                tagiter = leaditer;
                                probeiter = subleaditer;
                                // cout<<"the random number is "<< checkval<<" and even, leading jet is tagged"<<endl;
                            }
                            if(checkval%2!=0){
                                tagiter = subleaditer;
                                probeiter = leaditer;
                                // cout<<"the random number is "<< checkval<<" and odd, subleading jet is tagged"<<endl;
                            }
                        }
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
                }
            
                // finding the A values iff the leading jet has eta < 1.3 and subleading jet passes the pt cut and has eta < 5.2 
                if(((TMath::Abs(jteta[leaditer]<1.3))||(TMath::Abs(jteta[subleaditer]<1.3)))&&(jtcorrpt[subleaditer]>ptcut)&&(nref<3)){
                    
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
            }
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        }
    }

    timer.StartSplit("after event loop");

    // closing the file I'm getting the information from
    fi->Close();
    //cout <<endl<< "the number of jets in the barrel with pT above 80 GeV was " << pt80barreljetnum << endl<<endl;

    // normalizing some hists with the integral function
    ///////////////////////
    normalizeh(hvz);
    normalizeh(hjtpt);
    normalizeh(hjtcorrpt);
    normalizeh(hrawpt);
    normalizeh(hjteta);
    normalizeh(hjtphi); 
    ///////////////////////

    //cout<<"line 514"<<endl<<endl;

    // making a new file to store all the histograms of interest in
    // TFile f("pt_balance.root","recreate");
    TFile *f1 = new TFile("pt_balance.root", "recreate");
    f1->cd();

    // writing the base variable hists to this new file
    //cout<<"line 521"<<endl<<endl;
    ///////////////////
    TString alghere = "ak4pf";
    //cout<<"line 524"<<endl<<endl;
    save_h1d(hvz, "hvz", f1);
    //cout<<"line 526"<<endl<<endl;
    save_h1d(hjtcorrpt, "hjtcorrpt", f1);
    save_h1d(hjteta, "hjeta", f1);
    save_h1d(hjtphi, "hjtphi", f1);
    save_h1d(hrawpt, "hrawpt", f1);
    
    ///////////////////

    // looping throught the pt slices
    for(unsigned int k=0; k<ptslicenum; k++){

        // saving the A vs eta plots for each pt slice
        TString bhtitle0 = Form("ptslicesA_ptbin%d_%.0f_%.0f",k,ptlow[k],pthigh[k]);
        save_h2d(ptslicesA[k], bhtitle0, f1);
        //cout<<"line 537"<<endl<<endl;

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
                if(etalow[l]<0){
                    TString bhtitle1 = Form("etaslicesA_etabin_%.0f_%.0f_n",TMath::Abs(etalow[l]*10),TMath::Abs(etahigh[l]*10));
                    save_h2d(etaslicesA[l], bhtitle1, f1);
                }
                if(etalow[l]>0||etalow[l]==0){
                    TString bhtitle1 = Form("etaslicesA_etabin_%.0f_%.0f",etalow[l]*10,etahigh[l]*10);
                    save_h2d(etaslicesA[l], bhtitle1, f1);
                }
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
            //cout << "reached line 735" <<endl;
            if(etalow[l]<0){
                TString htitle1 = Form("Aptslice_%.0f_%.0f_etabin_%.0f_%.0f_n",ptlow[k],pthigh[k],TMath::Abs(etalow[l]*10),TMath::Abs(etahigh[l]*10));
                TString htitle2 = Form("Aetaslice_%.0f_%.0f_n_ptbin_%.0f_%.0f",TMath::Abs(etalow[l]*10),TMath::Abs(etahigh[l]*10),ptlow[k],pthigh[k]);
                save_h1d(etaslices_of_ptslicesA[k][l], htitle1, f1);
                save_h1d(ptslices_of_etaslicesA[l][k], htitle2, f1);
            }
            //cout << "reached line 742" <<endl;
            if(etalow[l]>0||etalow[l]==0){
                TString htitle1 = Form("Aptslice_%.0f_%.0f_etabin_%.0f_%.0f",ptlow[k],pthigh[k],etalow[l]*10,etahigh[l]*10);
                TString htitle2 = Form("Aetaslice_%.0f_%.0f_ptbin_%.0f_%.0f",etalow[l]*10,etahigh[l]*10,ptlow[k],pthigh[k]);
                save_h1d(etaslices_of_ptslicesA[k][l], htitle1, f1);
                save_h1d(ptslices_of_etaslicesA[l][k], htitle2, f1);
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // finding stuff for tgraphs
            /////////////////////////////////////////////////////////////////////////
            // pt slices hists of <A> vs eta
            Double_t ptaavg = (etaslices_of_ptslicesA[k][l])->GetMean();
            Double_t ptaavgerr = (etaslices_of_ptslicesA[k][l])->GetMeanError();
            ptslicesAavg[k][l] = ptaavg;
            ptslicesR[k][l] = ((1+ptaavg)/(1-ptaavg));
            ptslicesAavgerr[k][l] = ptaavgerr;
            ptslicesRerr[k][l] = (ptaavgerr*2/((1-ptaavg)*(1-ptaavg)));

            // eta slices hists of <A> vs pt
            Double_t etaaavg = (ptslices_of_etaslicesA[l][k])->GetMean();
            Double_t etaaavgerr = (ptslices_of_etaslicesA[l][k])->GetMeanError();
            etaslicesAavg[l][k] = etaaavg;
            etaslicesR[l][k] = ((1+etaaavg)/(1-etaaavg));
            etaslicesAavgerr[l][k] = etaaavgerr;
            etaslicesRerr[l][k] = (etaaavgerr*2/((1-etaaavg)*(1-etaaavg)));
            /////////////////////////////////////////////////////////////////////////
            
            // print statements
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //cout<<"<A> for η bin "<<etalow[l]<<" to "<<etahigh[l]<<" of pt slice "<<ptlow[k]<<" to "<<pthigh[k]<<" is "<<ptslicesAavg[k][l]<<"±"<<ptslicesAavgerr[k][l]<<endl;
            cout<<"R for η bin "<<etalow[l]<<" to "<<etahigh[l]<<" of pt slice "<<ptlow[k]<<" to "<<pthigh[k]<<" is "<<ptslicesR[k][l]<<"±"<<ptslicesRerr[k][l]<<endl;
            //cout<<"<A> for pt bin "<<ptlow[l]<<" to "<<pthigh[l]<<" of η slice "<<etalow[k]<<" to "<<etahigh[k]<<" is "<<etaslicesAavg[k][l]<<"±"<<etaslicesAavgerr[k][l]<<endl;
            cout<<"R for pt bin "<<ptlow[l]<<" to "<<pthigh[l]<<" of η slice "<<etalow[k]<<" to "<<etahigh[k]<<" is "<<etaslicesR[k][l]<<"±"<<etaslicesRerr[k][l]<<endl<<endl;
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //cout << " l is " << l << ", k is " << k<<endl<<endl;
        }
    }
    //cout << "reached line 779" <<endl;

    // actually finding the averages of the A values for the tgraphs
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    // <A> vs pt for each eta slice
    for(unsigned int k=0; k<etaslicenum; k++){

        // making 1D arrays for the TGraph function inputs
        Double_t ys[etaslicenum];
        Double_t yserr[etaslicenum];
        Double_t ys1[etaslicenum];
        Double_t yserr1[etaslicenum];

        // Assigning all the values to the 1D arrays
        for(unsigned int l=0; l<etaslicenum; l++){
            ys[l] = etaslicesAavg[k][l];
            yserr[l] = etaslicesAavgerr[k][l];
            ys1[l] = etaslicesR[k][l];
            yserr1[l] = etaslicesRerr[k][l];
        }

        // making the tgraphs
        getaslicesAavg[k] = new TGraphErrors(ptslicenum,etaslicesAx,ys,etaslicesAxerr,yserr);
        getaslicesR[k] = new TGraphErrors(ptslicenum,etaslicesAx,ys1,etaslicesAxerr,yserr1);
        // making the title for the traph
        if(etalow[k]<0){
            TString htitle3 = Form("eta_%.0f_%.0f_n",TMath::Abs(etalow[k]*10),TMath::Abs(etahigh[k]*10));
            // save_g_1(hpt, "pt", "A", getaslicesAavg[k], alghere, "p_{T}  [GeV/c]", "<A>", "<A>_"+htitle3);
            // save_g_1(hpt, "pt", "R", getaslicesR[k], alghere, "p_{T}  [GeV/c]", "R", "R_"+htitle3);
            save_g(getaslicesR[k],  "R_"+htitle3, f1);
        }
        if(etalow[k]>0||etalow[k]==0){
            TString htitle3 = Form("eta_%.0f_%.0f",etalow[k]*10,etahigh[k]*10);
            // save_g(getaslicesAavg[k], alghere, "p_T [GeV/c]", "<A>", htitle3);
            // save_g_1(hpt, "pt", "A", getaslicesAavg[k], alghere, "p_{T}  [GeV/c]", "<A>", "<A>_"+htitle3);
            // save_g_1(hpt, "pt", "R", getaslicesR[k], alghere, "p_{T}  [GeV/c]", "R", "R_"+htitle3);
            save_g(getaslicesR[k],  "R_"+htitle3, f1);
        }
    }

    // <A> vs eta for each pt slice
    // same as previous loop but for switched eta and pt
    for(unsigned int k=0; k<ptslicenum; k++){
        Double_t ys[ptslicenum];
        Double_t yserr[ptslicenum];
        Double_t ys1[ptslicenum];
        Double_t yserr1[ptslicenum];
        for(unsigned int l=0; l<ptslicenum; l++){
            ys[l] = ptslicesAavg[k][l];
            yserr[l] = ptslicesAavgerr[k][l];
            ys1[l] = ptslicesR[k][l];
            yserr1[l] = ptslicesRerr[k][l];
        }
        gptslicesAavg[k] = new TGraphErrors(etaslicenum,ptslicesAx,ys,ptslicesAxerr,yserr);
        gptslicesR[k] = new TGraphErrors(etaslicenum,ptslicesAx,ys1,ptslicesAxerr,yserr1);
        // getaslicesAavg[k]->SetMinimum(-0.1);
        // getaslicesAavg[k]->SetMaximum(0.1);
        TString hnamea = Form("pt_%.0f_%.0f",ptlow[k],pthigh[k]);
        // getaslicesAavg[k]->SetTitle(hname);
        // getaslicesAavg[k]->SetName(hname);
        // getaslicesAavg[k]->Write();
        // save_g(gptslicesAavg[k], alghere, "η", "<A>", hname);
        // save_g_1(heta, "eta", "A", gptslicesAavg[k], alghere, "η", "<A>", "<A>_"+hname);
        // save_g_1(heta, "eta", "R", gptslicesR[k], alghere, "η", "R", "R_"+hname);
        save_g(gptslicesR[k],  "R_"+hnamea, f1);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    f1->Close();
    // FINAL PRINT STATEMENTS

    // expressing passed/total as a percent
    // float p0 = b_/a_*100;

    // printing out the event cut amount
    // //cout<<"number of total events: "<< a_ << endl;
    // //cout<<"number of passed events before event cut: "<< b_ << endl;
    // //cout<<"percent of passed events before eta cut: "<< p0 <<"%"<< endl;
    //
    // determining the final time of the script
    timer.Stop();
    //cout<<endl;
    timer.Report();
}
