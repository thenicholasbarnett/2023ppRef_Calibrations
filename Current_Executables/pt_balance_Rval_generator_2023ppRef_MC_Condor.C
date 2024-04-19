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

#include "/afs/cern.ch/user/n/nbarnett/public/header_files/JetUncertainty.h"
#include "/afs/cern.ch/user/n/nbarnett/public/header_files/JetCorrector.h"

// PLOTTING FUNCTIONS

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

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// the script all runs in this function
void pt_balance_Rval_generator_2023ppRef_MC_Condor(TString input, TString output)
{   
    // getting rid of legends in hists in root file
    gStyle->SetOptStat(0);
    
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2(); 
    
    // INITIALIZING HISTOGRAMS

    // creating some binning parameters
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    double vzh1d0[3] = {40,-20,20};
    double pth1d0[3] = {100,15,500};
    double phih1d0[3] = {100,4,4};
    double etah1d0[3] = {50,-5.2,5.2};
    double etah1d1[3] = {25,-1.7,1.7};
    double ah1d0[3] = {100,-1,1};
    
    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 10;
    // the low and high pt values for each pt slice
    double ptlow[ptslicenum] = {15,25,50,80,100,120,140,180,220,300};
    double pthigh[ptslicenum] = {25,50,80,100,120,140,180,220,300,500};
    
    // eta slices
    // number of eta slices
    const Int_t etaslicenum = 30;
    // the low and high eta values for each eta slice
    // double oldetas[etaslicenum] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};
    double etahigh[etaslicenum] = {-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5,5.2};
    double etalow[etaslicenum] = {-5.2,-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5};
    

    // for generating random numbers
    TRandom2 *rand = new TRandom2(1);

    // Making some hists for pt balance studies
    // the actual A and R values
    // <A> as well as R vs pt and eta for eta and pt slices respectively
    // eta and pt bin slices for the pt slice hists of A vs eta_probe and pt_avg
    // pt slices
    TH2D *ptslicesA[ptslicenum];
    TH1D *etaslices_of_ptslicesA[ptslicenum][etaslicenum];
    TGraph *gptslicesAavg[ptslicenum];
    TGraph *gptslicesR[ptslicenum];
    // eta slices
    TH2D *etaslicesA[etaslicenum];
    TH1D *ptslices_of_etaslicesA[etaslicenum][ptslicenum];
    TGraph *getaslicesAavg[etaslicenum];
    TGraph *getaslicesR[etaslicenum];
    
    // further initializing the histograms
    // looping over pt slices
    for(unsigned int p=0; p<ptslicenum; p++){
        // A vs eta for each pt slice with title having the high and low pt for the slice
        TString chtitle0 = Form("A_ptslice_%.0f_%.0f",ptlow[p],pthigh[p]);
        ptslicesA[p] = new TH2D(chtitle0,chtitle0,etah1d0[0],etah1d0[1],etah1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
        for(unsigned int q=0; q<etaslicenum; q++){
            // avoiding looping through the etaslicenum separately by only doing it on the first p value
            if(p==0){
                // A vs pt for each eta slice with title having the high and low eta for the slice
                if(etalow[q]<0){
                    TString ahtitle0 = Form("A_etaslice_%.0f_%.0f_n",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                    etaslicesA[q] = new TH2D(ahtitle0,ahtitle0,pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
                }
                if(etalow[q]>0||etalow[q]==0){
                    TString ahtitle0 = Form("A_etaslice_%.0f_%.0f",etalow[q]*10,etahigh[q]*10);
                    etaslicesA[q] = new TH2D(ahtitle0,ahtitle0,pth1d0[0],pth1d0[1],pth1d0[2],ah1d0[0],ah1d0[1],ah1d0[2]);
                }
            }
            // intializing hists that are the bins of the slices
            if(etalow[q]<0){
                TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                TString dhtitle1 = Form("A_etaslice_%.0f_%.0f_n__ptbin_%.0f_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),ptlow[p],pthigh[p]);
                etaslices_of_ptslicesA[p][q] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                ptslices_of_etaslicesA[q][p] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
            }
            if(etalow[q]>0||etalow[q]==0){
                TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10);
                TString dhtitle1 = Form("A_etaslice_%.0f_%.0f__ptbin_%.0f_%.0f",etalow[q]*10,etahigh[q]*10,ptlow[p],pthigh[p]);
                etaslices_of_ptslicesA[p][q] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                ptslices_of_etaslicesA[q][p] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
            }
        }
    }

    // initializing histograms for general parameters
    TH1D *hvz = new TH1D("vz","vz",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hjtcorrpt = new TH1D("hjtcorrpt","hjtcorrpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","hjteta",etah1d1[0],etah1d1[1],etah1d1[2]);
    TH1D *hjtphi = new TH1D("hjtphi","hjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *hrawpt = new TH1D("hrawpt","hrawpt",pth1d0[0],pth1d0[1],pth1d0[2]);

    // initializing hists for canvases later
    TH1D *hpt = new TH1D("hpt","hpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *heta = new TH1D("heta","heta",etah1d0[0],etah1d0[1],etah1d0[2]);

    // INITIALIZING PARAMETERS

    // a_ and b_ are the number of total events and accepted events respectively
    Float_t a_=0, b_=0;

    // declaring variables
    // vertex position
    Float_t vz;
    // whatever ptcutm is decided
    Float_t ptcut = 1;
    // a big number to make my arrays such that they aren't too small
    // need to have more entries in the arrays than number of jets in the event with most jets
    const Int_t nm = 200000;
    // uncorrected jet pt
    Float_t rawpt[nm];
    // corrected jet pt
    Float_t jtcorrpt[nm];
    // jet phi and pseudorapadity
    Float_t jtphi[nm];
    Float_t jteta[nm];
    // number of jets in event
    Int_t nref;

    // pt balance study arrays
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

    // reading and iterating through a list of files instead of a single file
    TFile *fi = TFile::Open(input,"read");

    // getting the TTrees from the file
    // event info of interest in the following tree
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");

    // event cut info in the following tree
    TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");
    
    // turning off all branches and then turning on only the ones I want
    t0->SetBranchStatus("*",0);
    t0->SetBranchStatus("jteta",1);
    t0->SetBranchStatus("jtphi",1);
    t0->SetBranchStatus("rawpt",1);
    t0->SetBranchStatus("nref",1);
    
    // doing the same for the other tree I want
    t1->SetBranchStatus("*",0);
    t1->SetBranchStatus("vz",1);

    // general parameters
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);
    t0->SetBranchAddress("rawpt",rawpt);
    t0->SetBranchAddress("nref",&nref);
    t1->SetBranchAddress("vz",&vz);

    // for loop going over events in the trees
    for(unsigned int i=0; i<t0->GetEntries(); i++){
    //for(unsigned int i=0; i<100; i++){

        // cout<< "event " << i << " is being processed" << endl;

        // adding one to total events for every event
        a_+=1;

        // EVENT CUT
        // only events with |vz|<15 are passed
        t1->GetEntry(i);
        if(TMath::Abs(vz)<15){
            t0->GetEntry(i);
            
            // adding one to passed events iff all the conditionals are true
            b_+=1;

            // filling the vertex position hist
            hvz->Fill(vz);
            if(nref<2){continue;}
            // looping through all jets in each event
            for(unsigned int j=0; j<nref; j++){
                
                if((rawpt[j]>pth1d0[1])&&(rawpt[j]<pth1d0[2])){
                    hrawpt->Fill(rawpt[j]);
                }

                // getting the corrected jet pt
                vector<string> Files;
                Files.push_back("/afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2023ppRef/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");
		        JetCorrector JEC(Files);

                JEC.SetJetPT(rawpt[j]);
                JEC.SetJetEta(jteta[j]);
                JEC.SetJetPhi(jtphi[j]);  
                Float_t jet_pt_corr = JEC.GetCorrectedPT();
                // saving the corrected jet pt
                jtcorrpt[j] = jet_pt_corr;

                // Filling some hists
                if((jtcorrpt[j]>pth1d0[1])&&(jtcorrpt[j]<pth1d0[2])){
                    hjtcorrpt->Fill(jtcorrpt[j]);
                }

                // only look at pt balance studies if jtcorrpt > ptcut
                if(jtcorrpt[j]>ptcut){
                    // filling histograms that have variables with more than one value per event
                    hjteta->Fill(jteta[j]);
                    hjtphi->Fill(jtphi[j]);
                }
            }
            
            // PT BALANCE
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
            if(TMath::Abs(jteta[leaditer])<1.3){
                tagiter = leaditer;
                probeiter = subleaditer;
            }
            // if the subleading jet is in the barrel, tag it
            if(TMath::Abs(jteta[subleaditer])<1.3){
                tagiter = subleaditer;
                probeiter = leaditer;
            }
            // in the case both jets are in the barrel we make a random number 
            // if the random number is even or odd the tag jet is the leading or subleading jet respectively
            if((TMath::Abs(jteta[leaditer])<1.3)&&(TMath::Abs(jteta[subleaditer])<1.3)){
                int checkval1 = rand->Integer(100);
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
                for(unsigned int j=0; j<nref; j++){
                    // determining if another jet between index 3 and nref-1 exists with higher pt than the current third jet, but only if it is less pt than and isn't the probe or tag jet
                    if((jtcorrpt[j]>jtcorrpt[thirditer])&&(jtcorrpt[j]<jtcorrpt[subleaditer])&&(jtcorrpt[j]<jtcorrpt[leaditer])&&(j!=subleaditer)&&(j!=leaditer)){
                        thirditer = j;
                    }
                }
                // still working if there is a third jet
                // doing the whole pt balance study in the case there is a third jet
                if(((jtcorrpt[subleaditer]+jtcorrpt[leaditer])*0.1>jtcorrpt[thirditer])&&(TMath::Abs(jtphi[leaditer]-jtphi[subleaditer])>2.7)){
                    if(TMath::Abs(jteta[leaditer])<1.3){
                        tagiter = leaditer;
                        probeiter = subleaditer;
                    }
                    if(TMath::Abs(jteta[subleaditer])<1.3){
                        tagiter = subleaditer;
                        probeiter = leaditer;
                    }
                    if((TMath::Abs(jteta[leaditer])<1.3)&&(TMath::Abs(jteta[subleaditer])<1.3)){
                        int checkval = rand->Integer(100);
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
                    double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);
                    cout << "A is " << Aval << " for event " << i << endl;
                    // pt slices A value filling
                    // each k is a different slice of pt
                    for(unsigned int p=0; p<ptslicenum; p++){
                        // average momentum between the probe and tag jet 
                        // these are sliced originally to get A vs eta for different pt slices
                        double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                        // if the pt avg is within a certain slice range then add it to the pt slice hists
                        if((ptavg>ptlow[p])&&(ptavg<pthigh[p])){
                            // ptslicesA are A vs eta_probe hists for different pt ranges or slices
                            ptslicesA[p]->Fill(jteta[probeiter],Aval);
                        }
                    }
                    // eta slices A value filling
                    // each k is a different slice of eta
                    for(unsigned int q=0; q<etaslicenum; q++){
                        // ptavg is the x axis in one desired type of plot
                        double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                        // if the eta is within a certain slice range then add it to the eta slice hists
                        if((jteta[probeiter]>etalow[q])&&(jteta[probeiter]<etahigh[q])){
                            // etaslicesA are A vs pt_avg hists for different eta ranges
                            etaslicesA[q]->Fill(ptavg,Aval);
                        }
                    }
                }              
            }
            // finding the A values iff the leading jet has eta < 1.3 and subleading jet passes the pt cut and has eta < 5.2 
            if(((TMath::Abs(jteta[leaditer])<1.3)||(TMath::Abs(jteta[subleaditer])<1.3))&&(jtcorrpt[subleaditer]>ptcut)&&(nref<3)&&(TMath::Abs(jtphi[leaditer]-jtphi[subleaditer])>2.7)){
                // printing A value and saving it
                double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);
                cout << "A is " << Aval << " for event " << i << endl;
                // pt slices A value filling
                // each k is a different slice of pt
                for(unsigned int p=0; p<ptslicenum; p++){
                    // average momentum between the probe and tag jet 
                    // these are sliced originally to get A vs eta for different pt slices
                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    // if the pt avg is within a certain slice range then add it to the pt slice hists
                    if((ptavg>ptlow[p])&&(ptavg<pthigh[p])){
                        // ptslicesA are A vs eta_probe hists for different pt ranges or slices
                        ptslicesA[p]->Fill(jteta[probeiter],Aval);
                    }
                }
                // eta slices A value filling
                // each k is a different slice of eta
                for(unsigned int q=0; q<etaslicenum; q++){
                    // ptavg is the x axis in one desired type of plot
                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                    // if the eta is within a certain slice range then add it to the eta slice hists
                    if((jteta[probeiter]>etalow[q])&&(jteta[probeiter]<etahigh[q])){
                        // etaslicesA are A vs pt_avg hists for different eta ranges
                        etaslicesA[q]->Fill(ptavg,Aval);
                    }
                }
            }
        }
    }
    // closing the file I'm getting the information from
    fi->Close();
    
    // normalizing some hists with the integral function
    normalizeh(hvz);
    normalizeh(hjtcorrpt);
    normalizeh(hrawpt);
    normalizeh(hjteta);
    normalizeh(hjtphi); 

    // making a new file to store all the histograms of interest in
    TFile *f1 = new TFile(output,"recreate");
    f1->cd();

    // writing the base variable hists to this new file
    save_h1d(hvz, "hvz");
    save_h1d(hjtcorrpt, "hjtcorrpt");
    save_h1d(hjteta, "hjeta");
    save_h1d(hjtphi, "hjtphi");
    save_h1d(hrawpt, "hrawpt");
    
    // looping throught the pt slices
    for(unsigned int p=0; p<ptslicenum; p++){

        // saving the A vs eta plots for each pt slice
        TString bhtitle0 = Form("ptslicesA_%.0f_%.0f",ptlow[p],pthigh[p]);
        save_h2d(ptslicesA[p], bhtitle0);

        // making the x axis points for the eta slices be the center of each pt bin
        etaslicesAx[p] = (pthigh[p] + ptlow[p])/2;
        // making the x axis points error for the eta slices be half the width of each pt bin
        etaslicesAxerr[p] = etaslicesAx[p] - ptlow[p];

        // looping through the eta slices
        for(unsigned int q=0; q<etaslicenum; q++){

            // conditional below acts like an separated etanum loop 
            if(p==0){
                // saving eta slice hists of A vs pt
                if(etalow[q]<0){
                    TString bhtitle1 = Form("etaslicesA_%.0f_%.0f_n",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                    save_h2d(etaslicesA[q], bhtitle1);
                }
                if(etalow[q]>0||etalow[q]==0){
                    TString bhtitle1 = Form("etaslicesA_%.0f_%.0f",etalow[q]*10,etahigh[q]*10);
                    save_h2d(etaslicesA[q], bhtitle1);
                }
                // pt axis for tgraphs
                ptslicesAx[q] = (etahigh[q] + etalow[q])/2;
                ptslicesAxerr[q] = ptslicesAx[q] - etalow[q];
            }

            // getting y projection or slice of each eta bin for each pt slice
            etaslices_of_ptslicesA[p][q] = ptslicesA[p]->ProjectionY("",q,q,"");

            // getting y projection or slice of each pt bin for each eta slice
            ptslices_of_etaslicesA[q][p] = etaslicesA[q]->ProjectionY("",p,p,"");
            
            // saving the projection hists
            // eta bins of pt slices
            if(etalow[q]<0){
                TString htitle1 = Form("Aptslice_%.0f_%.0f_etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
                TString htitle2 = Form("Aetaslice_%.0f_%.0f_n_ptbin_%.0f_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),ptlow[p],pthigh[p]);
                save_h1d(etaslices_of_ptslicesA[p][q], htitle1);
                save_h1d(ptslices_of_etaslicesA[q][p], htitle2);
            }
            if(etalow[q]>0||etalow[q]==0){
                TString htitle1 = Form("Aptslice_%.0f_%.0f_etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10);
                TString htitle2 = Form("Aetaslice_%.0f_%.0f_ptbin_%.0f_%.0f",etalow[q]*10,etahigh[q]*10,ptlow[p],pthigh[p]);
                save_h1d(etaslices_of_ptslicesA[p][q], htitle1);
                save_h1d(ptslices_of_etaslicesA[q][p], htitle2);
            }

            // finding stuff for tgraphs
            // pt slices hists of <A> vs eta
            Double_t ptaavg = (etaslices_of_ptslicesA[p][q])->GetMean();
            Double_t ptaavgerr = (etaslices_of_ptslicesA[p][q])->GetMeanError();
            ptslicesAavg[p][q] = ptaavg;
            ptslicesR[p][q] = ((1+ptaavg)/(1-ptaavg));
            ptslicesAavgerr[p][q] = ptaavgerr;
            ptslicesRerr[p][q] = (ptaavgerr*2/((1-ptaavg)*(1-ptaavg)));

            // eta slices hists of <A> vs pt
            Double_t etaaavg = (ptslices_of_etaslicesA[q][p])->GetMean();
            Double_t etaaavgerr = (ptslices_of_etaslicesA[q][p])->GetMeanError();
            etaslicesAavg[q][p] = etaaavg;
            etaslicesR[q][p] = ((1+etaaavg)/(1-etaaavg));
            etaslicesAavgerr[q][p] = etaaavgerr;
            etaslicesRerr[q][p] = (etaaavgerr*2/((1-etaaavg)*(1-etaaavg)));
            
            // print statements
            //cout<<"<A> for η bin "<<etalow[q]<<" to "<<etahigh[q]<<" of pt slice "<<ptlow[p]<<" to "<<pthigh[p]<<" is "<<ptslicesAavg[p][q]<<"±"<<ptslicesAavgerr[p][q]<<endl;
            // cout<<"R for η bin "<<etalow[q]<<" to "<<etahigh[q]<<" of pt slice "<<ptlow[p]<<" to "<<pthigh[p]<<" is "<<ptslicesR[p][q]<<"±"<<ptslicesRerr[p][q]<<endl;
            //cout<<"<A> for pt bin "<<ptlow[q]<<" to "<<pthigh[q]<<" of η slice "<<etalow[p]<<" to "<<etahigh[p]<<" is "<<etaslicesAavg[p][q]<<"±"<<etaslicesAavgerr[p][q]<<endl;
            // cout<<"R for pt bin "<<ptlow[p]<<" to "<<pthigh[p]<<" of η slice "<<etalow[q]<<" to "<<etahigh[q]<<" is "<<etaslicesR[q][p]<<"±"<<etaslicesRerr[q][p]<<endl<<endl;
        }
    }
    // actually finding the averages of the A values for the tgraphs
    // <A> vs pt for each eta slice
    for(unsigned int q=0; q<etaslicenum; q++){
        // making 1D arrays for the TGraph function inputs
        // each eta slice q gets its own y axis
        Double_t ys[ptslicenum];
        Double_t yserr[ptslicenum];
        Double_t ys1[ptslicenum];
        Double_t yserr1[ptslicenum];
        // Assigning all the values to the 1D arrays
        // each pt bin in the y axis is the <A> or R for pt bin p of eta slice q
        for(unsigned int p=0; p<ptslicenum; p++){
            ys[p] = etaslicesAavg[q][p];
            yserr[p] = etaslicesAavgerr[q][p];
            ys1[p] = etaslicesR[q][p];
            yserr1[p] = etaslicesRerr[q][p];
        }
        // making the tgraphs
        getaslicesAavg[q] = new TGraphErrors(ptslicenum,etaslicesAx,ys,etaslicesAxerr,yserr);
        getaslicesR[q] = new TGraphErrors(ptslicenum,etaslicesAx,ys1,etaslicesAxerr,yserr1);
        // making the title for the traph
        if(etalow[q]<0){
            TString htitle3 = Form("eta_%.0f_%.0f_n",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10));
            save_g(getaslicesR[q],  "R_"+htitle3);
        }
        if(etalow[q]>0||etalow[q]==0){
            TString htitle3 = Form("eta_%.0f_%.0f",etalow[q]*10,etahigh[q]*10);
            save_g(getaslicesR[q],  "R_"+htitle3);
        }
    }
    // <A> vs eta for each pt slice
    // same as previous loop but for switched eta and pt
    for(unsigned int p=0; p<ptslicenum; p++){
        Double_t ys[etaslicenum];
        Double_t yserr[etaslicenum];
        Double_t ys1[etaslicenum];
        Double_t yserr1[etaslicenum];
        for(unsigned int q=0; q<etaslicenum; q++){
            ys[q] = ptslicesAavg[p][q];
            yserr[q] = ptslicesAavgerr[p][q];
            ys1[q] = ptslicesR[p][q];
            yserr1[q] = ptslicesRerr[p][q];
        }
        gptslicesAavg[p] = new TGraphErrors(etaslicenum,ptslicesAx,ys,ptslicesAxerr,yserr);
        gptslicesR[p] = new TGraphErrors(etaslicenum,ptslicesAx,ys1,ptslicesAxerr,yserr1);
        TString hnamea = Form("pt_%.0f_%.0f",ptlow[p],pthigh[p]);
        save_g(gptslicesR[p],  "R_"+hnamea);
    }

    f1->Close();
}
