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
    // l->Draw("same");
    c->Write();
    h->Write();
    c->SaveAs("plots_5_19/"+hname+".png");
    delete c;
    delete h_c;
}

void save_h1d_3(TH1D *h1, TH1D *h2, TH1D *h3, TString xtitle, TString ytitle, TString hname, TString hname1, TString hname2, TString hname3){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    TH1D *h1_c = (TH1D*)h1->Clone(hname1);
    TH1D *h2_c = (TH1D*)h2->Clone(hname2);
    TH1D *h3_c = (TH1D*)h3->Clone(hname3);
    h1_c->Draw("e1p");
    h2_c->Draw("same");
    h3_c->Draw("same");
    h1_c->SetMarkerStyle(20);
    h1_c->SetMarkerColor(kBlack);
    h1_c->SetLineColor(kBlack);
    h1_c->SetTitle("");
    h1_c->SetName(hname);
    h1_c->GetYaxis()->SetTitle(ytitle);
    h1_c->GetYaxis()->CenterTitle(true);
    h1_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h1_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h1_c->GetXaxis()->SetTitle(xtitle);}
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h1_c,hname1,"pl");
    l->AddEntry(h2_c,hname2,"pl");
    l->AddEntry(h3_c,hname3,"pl");
    l->Draw("same");
    c->Write();
    // h1->Write();
    c->SaveAs("plots_5_19/"+hname+".png");
    delete c;
}

void save_h1d_2(TH1D *h1, TH1D *h2, TString xtitle, TString ytitle, TString hname, TString hname1, TString hname2){
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);
    TH1D *h1_c = (TH1D*)h1->Clone(hname1);
    TH1D *h2_c = (TH1D*)h2->Clone(hname2);
    h1_c->Draw("e1p");
    h2_c->Draw("same");
    h1_c->SetMarkerStyle(20);
    h1_c->SetMarkerColor(kBlack);
    h1_c->SetLineColor(kBlack);
    h2_c->SetLineColor(kRed);
    h1_c->SetTitle("");
    h1_c->SetName(hname);
    h1_c->GetYaxis()->SetTitle(ytitle);
    h1_c->GetYaxis()->CenterTitle(true);
    h1_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h1_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h1_c->GetXaxis()->SetTitle(xtitle);}
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h1_c,hname1,"pl");
    l->AddEntry(h2_c,hname2,"pl");
    l->Draw("same");
    c->Write();
    // h1->Write();
    c->SaveAs("plots_5_19/"+hname+".png");
    delete c;
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
    // TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    // l->SetBorderSize(0);
    // l->SetFillStyle(0);
    // l->AddEntry(h_c,hname,"pl");
    // l->Draw("same");
    c->Write();
    h->Write();
    c->SaveAs("plots_5_19/"+hname+".png");
    delete c;
    delete h_c;
}

// OTHER FUNCTIONS //

// normalizing a TH1D by integrating and scaling by the inverse of the integration value
void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1/a);
}

// VISUALIZATION FUNCTIONS //

// another try at the visualization tool
void visualize_subjets(Double_t jet_mom[3], Double_t subjet_mom[3][30], TH2D *h, TString hname){
    
    // // getting rid of automatic legend
    // gStyle->SetOptStat(0);

    // radius of subjet
    Double_t r_sj = 0.1;
    
    // radius of jet
    Double_t r_j = 0.4;

    // making the canvas
    TCanvas *c = new TCanvas();
    c->cd();
    // c->Modified();
    // c->SetGrid();
    c->SetTitle("");
    c->SetName(hname);
    // c->Range(-3.2,-3.2,3.2,3.2);
    c->Range(-5.2,-3.2,5.2,3.2);

    // // drawing the background histogram for the correct axis and titles
    // TH2D *h_c = (TH2D*)h->Clone("");
    // h_c->SetTitle("");
    // h_c->GetXaxis()->SetTitle("#eta");
    // h_c->GetYaxis()->SetTitle("#phi");
    // h_c->GetYaxis()->CenterTitle(true);
    // h_c->GetXaxis()->CenterTitle(true);
    // h_c->Draw("COLZ ");

    // try with 1D hist

    // making the circle for the jet
    auto elj1 = new TEllipse(jet_mom[1], jet_mom[2], r_j, r_j);

    // setting the outline color for the jet circle
    elj1->SetLineColor(1);
    
    // // drawing the circle for the jet
    // if(jet_mom[1]>-2.4&&jet_mom[1]<2.4){
    //     elj1->Draw("same");
    // }
    elj1->Draw("same");

    // // printing out jet momenta
    // if((jet_mom[0]!=0)||(jet_mom[0]!=-999)){
    //     cout<<"jet pt is "<<jet_mom[0]<<endl;
    //     cout<<"jet eta is "<<jet_mom[1]<<endl;
    //     cout<<"jet phi is "<<jet_mom[2]<<endl;
    // }

    // pointing to TObjects from TEllipse
    TEllipse *els[30];
    
    for(unsigned int i=0; i<30; i++){

        // only going through subjets that aren't flagged as bad
        if((subjet_mom[0][i]==-999)||(subjet_mom[0][i]==0)||(jet_mom[0]==0)||(jet_mom[0]==-999)){continue;}

        // // printing out jet momenta
        // cout<<"subjet pt is "<<subjet_mom[0][i]<<endl;
        // cout<<"subjet eta is "<<subjet_mom[1][i]<<endl;
        // cout<<"subjet phi is "<<subjet_mom[2][i]<<endl;

        // make and drawing the subjet ellipse
        els[i] = new TEllipse(subjet_mom[1][i], subjet_mom[2][i], r_sj, r_sj);
        els[i]->SetLineColor(2);
        els[i]->Draw("same");
        // if(jet_mom[1]>-2.4&&jet_mom[1]<2.4&&subjet_mom[1][i]>-2.4&&subjet_mom[2][i]<2.4){
        //     els[i]->Draw("same");
        // }
    }

    // // writing out the plot
    // if(jet_mom[1]>-2.4&&jet_mom[1]<2.4){
    //     c->Write();
    // }
    c->Write();

    // deleting canvas and TEllipses to avoid abusing memory
    delete c;
    delete elj1;

}

// another try at the visualization tool
void visualize_subjets1(Double_t jet_mom[3][10], Double_t subjet_mom[3][10][30], int nevent){
    
    // // getting rid of automatic legend
    // gStyle->SetOptStat(0);

    // making sure there is at least one jet and one subjet in the canvas before saving it
    double jnum = 0;
    double sjnum = 0;

    // CONSTANTS // 

    // radius of subjet
    const Double_t r_sj = 0.1;
    
    // radius of jet
    const Double_t r_j = 0.4;
    
    // max number of jets plotted here
    const int nj = 10;
    
    // max number of subjets plotted here
    const int nsj = 30;

    // CANVAS //

    // making the name of the visualization
    TString hname = Form("subjet_visualization_event%d",nevent);

    // making the canvas
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);

    // // Range of interest
    // c->Range(-2.7,-3.2,2.7,3.2);
    
    // entire possible range of values
    c->Range(-5.2,-1*TMath::Pi(),5.2,TMath::Pi());

    // // other possible options that may help with adding eta phi map
    // c->Modified();
    // c->SetGrid();

    // JET ELLIPSES //
    
    // pointing to TObjects for jet TEllipses
    TEllipse *elj[nj];

    // pointing to TObjects for subjet TEllipses
    TEllipse *els[nj][nsj];

    for(unsigned int j=0; j<nj; j++){

        // not drawing the jet if it's flagged as bad
        if((jet_mom[0][j]==0)||(jet_mom[0][j]==-999)){continue;}

        // adding to jet number if there is one
        jnum +=1;
        
        // making the circle for the jet
        elj[j] = new TEllipse(jet_mom[1][j], jet_mom[2][j], r_j, r_j);

        // setting the outline color for the jet circle
        elj[j]->SetLineColor(1);

        // drawing the jet
        elj[j]->Draw("same");
        
        // // printing out jet momenta
        // cout<<"pt is "<<jet_mom[0][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        // cout<<"eta is "<<jet_mom[1][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        // cout<<"phi is "<<jet_mom[2][j]<<" for jet "<<j<<" in event "<<nevent<<endl;

        // printing out jet momenta if it's eta or phi value doesn't make sense
        if((TMath::Abs(jet_mom[1][j])>5.2)||(TMath::Abs(jet_mom[2][j])>TMath::Pi())){
            cout<<"pt is "<<jet_mom[0][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
            cout<<"eta is "<<jet_mom[1][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
            cout<<"phi is "<<jet_mom[2][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        }
    
        // SUBJET ELLIPSES //

        for(unsigned int s=0; s<nsj; s++){

            // only going through subjets that aren't flagged as bad
            if((subjet_mom[0][j][s]==-999)||(subjet_mom[0][j][s]==0)){continue;}

            // adding to jet number if there is one
            sjnum +=1;

            // make and drawing the subjet ellipse
            els[j][s] = new TEllipse(subjet_mom[1][j][s], subjet_mom[2][j][s], r_sj, r_sj);
            els[j][s]->SetLineColor(2);
            els[j][s]->Draw("same");

            // // printing out subjet momenta if it's eta or phi value doesn't make sense
            if((TMath::Abs(subjet_mom[1][j][s])>5.2)||(TMath::Abs(subjet_mom[2][j][s])>TMath::Pi())){
                cout<<"pt is "<<subjet_mom[0][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
                cout<<"eta is "<<subjet_mom[1][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
                cout<<"phi is "<<subjet_mom[2][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
            }
        }
    }

    // saving the canvas to the file of interest, specifically the file that is staged
    if((jnum>0)&&(sjnum>0)){c->Write();}

    // deleting canvas to avoid abusing memory
    delete c;

}

void visualize_subjets2(Double_t jet_mom[3][10], Double_t subjet_mom[3][10][30], Double_t constituent_mom[3][10][200], int nevent){
    
    // // getting rid of automatic legend
    // gStyle->SetOptStat(0);

    // making sure there is at least one jet and one subjet in the canvas before saving it
    double jnum = 0;
    double sjnum = 0;

    // CONSTANTS // 

    // radius of subjet
    const Double_t r_sj = 0.1;
    
    // radius of jet
    const Double_t r_j = 0.4;
    
    // max number of jets plotted here
    const int nj = 10;
    
    // max number of subjets plotted here
    const int nsj = 30;
    
    // max number of jet constituents plotted here
    const int njc = 200;

    // CANVAS //

    // making the name of the visualization
    TString hname = Form("subjet_visualization_event%d",nevent);

    // making the canvas
    TCanvas *c = new TCanvas();
    c->cd();
    c->SetTitle("");
    c->SetName(hname);

    // // Range of interest
    // c->Range(-2.7,-3.2,2.7,3.2);
    
    // entire possible range of values
    c->Range(-5.2,-1*TMath::Pi(),5.2,TMath::Pi());

    // // other possible options that may help with adding eta phi map
    // c->Modified();
    // c->SetGrid();

    // JET ELLIPSES //
    
    // pointing to TObjects for jet TEllipses
    TEllipse *elj[nj];

    // pointing to TObjects for subjet TEllipses
    TEllipse *els[nj][nsj];

    // pointing to TMarkers for jet constituents
    TMarker *mjc[nj][njc];

    for(unsigned int j=0; j<nj; j++){

        // not drawing the jet if it's flagged as bad
        if((jet_mom[0][j]==0)||(jet_mom[0][j]==-999)){continue;}

        // adding to jet number if there is one
        jnum +=1;
        
        // making the circle for the jet
        elj[j] = new TEllipse(jet_mom[1][j], jet_mom[2][j], r_j, r_j);

        // setting the outline color for the jet circle
        elj[j]->SetLineColor(1);

        // drawing the jet
        elj[j]->Draw("same");
        
        // // printing out jet momenta
        // cout<<"pt is "<<jet_mom[0][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        // cout<<"eta is "<<jet_mom[1][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        // cout<<"phi is "<<jet_mom[2][j]<<" for jet "<<j<<" in event "<<nevent<<endl;

        // printing out jet momenta if it's eta or phi value doesn't make sense
        if((TMath::Abs(jet_mom[1][j])>5.2)||(TMath::Abs(jet_mom[2][j])>TMath::Pi())){
            cout<<"pt is "<<jet_mom[0][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
            cout<<"eta is "<<jet_mom[1][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
            cout<<"phi is "<<jet_mom[2][j]<<" for jet "<<j<<" in event "<<nevent<<endl;
        }
    
        // SUBJET ELLIPSES //

        for(unsigned int s=0; s<nsj; s++){

            // only going through subjets that aren't flagged as bad
            if((subjet_mom[0][j][s]==-999)||(subjet_mom[0][j][s]==0)){continue;}

            // adding to jet number if there is one
            sjnum +=1;

            // make and drawing the subjet ellipse
            els[j][s] = new TEllipse(subjet_mom[1][j][s], subjet_mom[2][j][s], r_sj, r_sj);
            els[j][s]->SetLineColor(2);
            els[j][s]->Draw("same");

            // // printing out subjet momenta if it's eta or phi value doesn't make sense
            if((TMath::Abs(subjet_mom[1][j][s])>5.2)||(TMath::Abs(subjet_mom[2][j][s])>TMath::Pi())){
                cout<<"pt is "<<subjet_mom[0][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
                cout<<"eta is "<<subjet_mom[1][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
                cout<<"phi is "<<subjet_mom[2][j][s]<<" for jet "<<j<<" and subjet "<<s<<" in event "<<nevent<<endl;
            }
        }

        // CONSTITUENTS //

        for(unsigned int jc=0; jc < njc; jc++){
            // only going through constituents that aren't flagged as bad
            if((constituent_mom[0][j][jc]==-999)||(constituent_mom[0][j][jc]==0)){continue;}
            mjc[j][jc] = new TMarker(constituent_mom[1][j][jc],constituent_mom[2][j][jc],5);
            mjc[j][jc]->Draw("same");
        }
    }

    // saving the canvas to the file of interest, specifically the file that is staged
    if((jnum>0)&&(sjnum>0)){c->Write();}

    // deleting canvas to avoid abusing memory
    delete c;

}

void visualize_subjets3(Double_t jet_mom[2][3][10], Double_t subjet_mom[2][3][10][30], Double_t constituent_mom[2][3][10][200], int nevent){
    
    // // getting rid of automatic legend
    // gStyle->SetOptStat(0);

    // making sure there is at least one jet and one subjet in the canvas before saving it
    double jnum = 0;
    double sjnum = 0;

    // CONSTANTS // 

    // radius of subjet
    const Double_t r_sj = 0.1;
    
    // radius of jet
    const Double_t r_j = 0.4;
    
    // max number of jets plotted here
    const int nj = 10;
    
    // max number of subjets plotted here
    const int nsj = 30;
    
    // max number of jet constituents plotted here
    const int njc = 200;

    // CANVAS //

    // making the name of the visualizations
    TString hname0 = Form("subjet_visualization_event_%d_reco",nevent);
    TString hname1 = Form("subjet_visualization_event_%d_ref",nevent);
    TString hname2 = Form("subjet_visualization_event_%d",nevent);

    // making the canvas
    TCanvas *c0 = new TCanvas();
    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    c0->SetTitle("");
    c1->SetTitle("");
    c2->SetTitle("");
    c0->SetName(hname0);
    c1->SetName(hname1);
    c2->SetName(hname2);

    // // Range of interest
    // c0->Range(-2.7,-1*TMath::Pi(),2.7,TMath::Pi());
    // c1->Range(-2.7,-1*TMath::Pi(),2.7,TMath::Pi());
    // c2->Range(-2.7,-1*TMath::Pi(),2.7,TMath::Pi());
    
    // entire possible range of values
    c0->Range(-5.2,-1*TMath::Pi(),5.2,TMath::Pi());
    c1->Range(-5.2,-1*TMath::Pi(),5.2,TMath::Pi());
    c2->Range(-5.2,-1*TMath::Pi(),5.2,TMath::Pi());

    // JET ELLIPSES //
    
    // pointing to TObjects for jet TEllipses
    TEllipse *elj[2][nj];

    // pointing to TObjects for subjet TEllipses
    TEllipse *els[2][nj][nsj];

    // pointing to TMarkers for jet constituents
    TMarker *mjc[2][nj][njc];

    // for(unsigned int j=nj; j>-1; j--){
    for(unsigned int j=0; j<nj; j++){

        // not drawing the jet if it's flagged as bad
        if((jet_mom[0][0][j]==0)||(jet_mom[0][0][j]==-999)||(jet_mom[1][0][j]==0)||(jet_mom[1][0][j]==-999)){continue;}

        // adding to jet number
        jnum +=2;
        
        // making the circle for the jet
        // reco
        elj[0][j] = new TEllipse(jet_mom[0][1][j], jet_mom[0][2][j], r_j, r_j);
        // ref
        elj[1][j] = new TEllipse(jet_mom[1][1][j], jet_mom[1][2][j], r_j, r_j);

        // setting the outline color for the jet circles
        elj[0][j]->SetLineColor(2);
        elj[1][j]->SetLineColor(4);

        // drawing the jet
        // reco jet
        c0->cd();
        elj[0][j]->Draw("same");
        // ref jet
        c1->cd();
        elj[1][j]->Draw("same");
        // both jets
        c2->cd();
        elj[0][j]->Draw("same");
        elj[1][j]->Draw("same");

        // printing out jet momenta if it's eta or phi value doesn't make sense
        if((TMath::Abs(jet_mom[0][1][j])>5.2)||(TMath::Abs(jet_mom[0][2][j])>TMath::Pi())||(TMath::Abs(jet_mom[1][1][j])>5.2)||(TMath::Abs(jet_mom[1][2][j])>TMath::Pi())){
            cout<<"check event "<<nevent<<endl;
        }
    
        // SUBJET ELLIPSES //

        // for(unsigned int s=nsj; s>-1; s--){
        for(unsigned int s=0; s<nsj; s++){

            // only going through subjets that aren't flagged as bad
            if((subjet_mom[0][0][j][s]==-999)||(subjet_mom[0][0][j][s]==0)||(subjet_mom[1][0][j][s]==-999)||(subjet_mom[1][0][j][s]==0)){continue;}

            // adding to jet number if there is one
            sjnum +=2;

            // make and drawing the subjet ellipse
            els[0][j][s] = new TEllipse(subjet_mom[0][1][j][s], subjet_mom[0][2][j][s], r_sj, r_sj);
            els[1][j][s] = new TEllipse(subjet_mom[1][1][j][s], subjet_mom[1][2][j][s], r_sj, r_sj);
            els[0][j][s]->SetLineColor(6);
            els[1][j][s]->SetLineColor(7);
            c0->cd();
            els[0][j][s]->Draw("same");
            c1->cd();
            els[1][j][s]->Draw("same");
            c2->cd();
            els[0][j][s]->Draw("same");
            els[1][j][s]->Draw("same");

            // // printing out subjet momenta if it's eta or phi value doesn't make sense
            if((TMath::Abs(subjet_mom[0][1][j][s])>5.2)||(TMath::Abs(subjet_mom[0][2][j][s])>TMath::Pi())||(TMath::Abs(subjet_mom[1][1][j][s])>5.2)||(TMath::Abs(subjet_mom[1][2][j][s])>TMath::Pi())){
                cout<<"check event "<<nevent<<endl;
            }
        }

        // CONSTITUENTS //

        for(unsigned int jc=0; jc < njc; jc++){

            // only going through constituents that aren't flagged as bad
            if((constituent_mom[0][0][j][jc]==-999)||(constituent_mom[0][0][j][jc]==0)||(constituent_mom[1][0][j][jc]==-999)||(constituent_mom[1][0][j][jc]==0)){continue;}
            mjc[0][j][jc] = new TMarker(constituent_mom[0][1][j][jc],constituent_mom[0][2][j][jc],5);
            mjc[1][j][jc] = new TMarker(constituent_mom[1][1][j][jc],constituent_mom[1][2][j][jc],5);
            mjc[0][j][jc]->SetMarkerColor(3);
            if(constituent_mom[0][0][j][jc]<5){mjc[0][j][jc]->SetMarkerColor(8);}
            if(constituent_mom[0][0][j][jc]<2){mjc[0][j][jc]->SetMarkerColor(28);}
            mjc[1][j][jc]->SetMarkerColor(1);
            if(constituent_mom[1][0][j][jc]<5){mjc[1][j][jc]->SetMarkerColor(12);}
            if(constituent_mom[1][0][j][jc]<2){mjc[1][j][jc]->SetMarkerColor(28);}
            c0->cd();
            mjc[0][j][jc]->Draw("same");
            c1->cd();
            mjc[1][j][jc]->Draw("same");
            c2->cd();
            mjc[0][j][jc]->Draw("same");
            mjc[1][j][jc]->Draw("same");
        }
    }

    // saving the canvas to the file of interest, specifically the file that is staged
    if((jnum>0)&&(sjnum>0)){
        c0->Write();
        c1->Write();
        c2->Write();
    }

    // deleting canvas to avoid abusing memory
    delete c0;
    delete c1;
    delete c2;

}

// void minivisualization(){
//     //
// }

// MAIN CODE FUNCTION //

void subjet_analysis(){
    
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();

    // getting rid of auto legends
    gStyle->SetOptStat(0);

    // MAKING HIST BINS //

    // vz
    const double vzh1d0[3] = {40,-20,20};

    // numbers
    const double numh1d0[3] = {30,0,30};
    const double numh1d1[3] = {50,0,50};

    // momenta
    const double pth1d0[3] = {100,80,500};
    const double etah1d0[3] = {50,-5.2,5.2};
    const double phih1d0[3] = {100,-1*TMath::Pi(),TMath::Pi()};

    // MAKING SPECIFIC PT HIST BINNING //

    // number of pt slices
    const Int_t ptslicenum = 10;

    // the low and high pt values for each pt slice
    const double ptlow[ptslicenum] = {15,25,50,80,100,120,140,180,220,300};
    const double pthigh[ptslicenum] = {25,50,80,100,120,140,180,220,300,500};

    // INITIALIZING EVENT HISTS //

    // vz
    TH1D *hvz = new TH1D("vz","vz",vzh1d0[0],vzh1d0[1],vzh1d0[2]);

    // numbers of jets
    TH1D *hngen = new TH1D("hngen","hngen",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hnref = new TH1D("hnref","hnref",numh1d0[0],numh1d0[1],numh1d0[2]);
    
    // INITIALIZING JET HISTS //
    
    // Reco
    TH1D *hjtpt = new TH1D("hjtpt","hjtpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","hjteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hjtphi = new TH1D("hjtphi","hjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // Ref
    TH1D *hrefpt = new TH1D("hrefpt","hrefpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hrefeta = new TH1D("hrefeta","hrefeta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hrefphi = new TH1D("hrefphi","hrefphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // Gen
    TH1D *hgenpt = new TH1D("hgenpt","hgenpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hgeneta = new TH1D("hgeneta","hgeneta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hgenphi = new TH1D("hgenphi","hgenphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    
    // INITIALIZING SUBJET HISTS //

    // 1D HISTOGRAMS 
    
    // reco
    // number of subjets
    TH1D *hsubrecojetnum = new TH1D("hsubrecojetnum","hsubrecojetnum",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hsubrecojetnum_ptcut5 = new TH1D("hsubrecojetnum_ptcut5","hsubrecojetnum_ptcut5",numh1d0[0],numh1d0[1],numh1d0[2]);
    // subjets momenta
    TH1D *hsubrecojetpt = new TH1D("hsubrecojetpt","hsubrecojetpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hsubrecojeteta = new TH1D("hsubrecojeteta","hsubrecojeteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hsubrecojetphi = new TH1D("hsubrecojetphi","hsubrecojetphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    
    // ref
    // number of subjets
    TH1D *hsubrefjetnum = new TH1D("hsubrefjetnum","hsubrefjetnum",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hsubrefjetnum_ptcut5 = new TH1D("hsubrefjetnum_ptcut5","hsubrefjetnum_ptcut5",numh1d0[0],numh1d0[1],numh1d0[2]);
    // subjets momenta
    TH1D *hsubrefjetpt = new TH1D("hsubrefjetpt","hsubrefjetpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hsubrefjeteta = new TH1D("hsubrefjeteta","hsubrefjeteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hsubrefjetphi = new TH1D("hsubrefjetphi","hsubrefjetphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    
    // gen
    // number of subjets
    TH1D *hsubgenjetnum = new TH1D("hsubgenjetnum","hsubgenjetnum",numh1d0[0],numh1d0[1],numh1d0[2]);
    TH1D *hsubgenjetnum_ptcut5 = new TH1D("hsubgenjetnum_ptcut5","hsubgenjetnum_ptcut5",numh1d0[0],numh1d0[1],numh1d0[2]);
    // subjets momenta
    TH1D *hsubgenjetpt = new TH1D("hsubgenjetpt","hsubgenjetpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hsubgenjeteta = new TH1D("hsubgenjeteta","hsubgenjeteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hsubgenjetphi = new TH1D("hsubgenjetphi","hsubgenjetphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // 2D HISTOGRAMS

    // eta - phi plane
    TH2D *heta_phi = new TH2D("heta_phi","heta_phi",etah1d0[0],etah1d0[1],etah1d0[2],phih1d0[0],phih1d0[1],phih1d0[2]);

    // number vs pt for reco
    TH2D *hsubrecojetnum_jtpt = new TH2D("hsubrecojetnum_jtpt","hsubrecojetnum_jtpt",pth1d0[0],pth1d0[1],pth1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);

    // number vs pt cut for reco and ref
    TH2D *hsubrecojetnum_ptcut = new TH2D("hsubrecojetnum_ptcut","hsubrecojetnum_ptcut",numh1d1[0],numh1d1[1],numh1d1[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    TH2D *hsubrefjetnum_ptcut = new TH2D("hsubrefjetnum_ptcut","hsubrefjetnum_ptcut",numh1d1[0],numh1d1[1],numh1d1[2],numh1d0[0],numh1d0[1],numh1d0[2]);

    // number of reco subjets vs number of ref subjets
    TH2D *hsubrecojetnum_subrefjetnum = new TH2D("hsubrecojetnum_subrefjetnum","hsubrecojetnum_subrefjetnum",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    TH2D *hsubrecojetnum_subrefjetnum_ptcut5 = new TH2D("hsubrecojetnum_subrefjetnum_ptcut5","hsubrecojetnum_subrefjetnum_ptcut5",numh1d0[0],numh1d0[1],numh1d0[2],numh1d0[0],numh1d0[1],numh1d0[2]);
    
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
    const Int_t MAXSUBJETS = 30;
    const Int_t MAXCONSTITUENTS = 200;

    // JET VARIABLES
    
    // reco jet momenta
    Float_t jtpt[MAXJETS];
    Float_t jtphi[MAXJETS];
    Float_t jteta[MAXJETS];
    
    // ref jet momenta
    Float_t refpt[MAXJETS];
    Float_t refeta[MAXJETS];
    Float_t refphi[MAXJETS];
    
    // gen jet momenta
    Float_t genpt[MAXJETS];
    Float_t geneta[MAXJETS];
    Float_t genphi[MAXJETS];

    // SUBJET VARIABLES

    // reco
    // number of subjets
    Float_t subrecojetnum[MAXJETS];
    // subjets momenta
    Float_t subrecojetpt[MAXJETS][MAXSUBJETS];
    Float_t subrecojeteta[MAXJETS][MAXSUBJETS];
    Float_t subrecojetphi[MAXJETS][MAXSUBJETS];
    // constituents momenta
    Float_t recojetcpt[MAXJETS][MAXCONSTITUENTS];
    Float_t recojetceta[MAXJETS][MAXCONSTITUENTS];
    Float_t recojetcphi[MAXJETS][MAXCONSTITUENTS];

    // ref
    // number of subjets
    Float_t subrefjetnum[MAXJETS];
    // subjets momenta
    Float_t subrefjetpt[MAXJETS][MAXSUBJETS];
    Float_t subrefjeteta[MAXJETS][MAXSUBJETS];
    Float_t subrefjetphi[MAXJETS][MAXSUBJETS];
    // constituents momenta
    Float_t refjetcpt[MAXJETS][MAXCONSTITUENTS];
    Float_t refjetceta[MAXJETS][MAXCONSTITUENTS];
    Float_t refjetcphi[MAXJETS][MAXCONSTITUENTS];

    // gen
    // number of subjets
    Float_t subgenjetnum[MAXJETS];
    // subjets momenta
    Float_t subgenjetpt[MAXJETS][MAXSUBJETS];
    Float_t subgenjeteta[MAXJETS][MAXSUBJETS];
    Float_t subgenjetphi[MAXJETS][MAXSUBJETS];
    // constituents momenta
    Float_t genjetcpt[MAXJETS][MAXCONSTITUENTS];
    Float_t genjetceta[MAXJETS][MAXCONSTITUENTS];
    Float_t genjetcphi[MAXJETS][MAXCONSTITUENTS];

    // COUNTING VARIABLES

    // number of subjets
    Double_t subjetcount = 0;

    // number of subjets with pt > 1 TeV
    Double_t subjetcount_hugept = 0;

    // INPUT //

    // INPUT FILE
    TFile *fi = TFile::Open("HiForestAOD_5_19_2024_1.root","read");
    
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

    // Reco
    t0->SetBranchStatus("jtpt",1);
    t0->SetBranchStatus("jteta",1);
    t0->SetBranchStatus("jtphi",1);

    // Ref
    t0->SetBranchStatus("refpt",1);
    t0->SetBranchStatus("refeta",1);
    t0->SetBranchStatus("refphi",1);

    // Gen
    t0->SetBranchStatus("genpt",1);
    t0->SetBranchStatus("geneta",1);
    t0->SetBranchStatus("genphi",1);

    // SUBJET BRANCHES

    // Reco
    // number of subjets
    t0->SetBranchStatus("subrecojetnum",1);
    // subjets momenta
    t0->SetBranchStatus("subrecojetpt",1);
    t0->SetBranchStatus("subrecojeteta",1);
    t0->SetBranchStatus("subrecojetphi",1);
    // constituent momenta
    t0->SetBranchStatus("recojetcpt",1);
    t0->SetBranchStatus("recojetceta",1);
    t0->SetBranchStatus("recojetcphi",1);

    // Ref
    // number of subjets
    t0->SetBranchStatus("subrefjetnum",1);
    // subjets momenta
    t0->SetBranchStatus("subrefjetpt",1);
    t0->SetBranchStatus("subrefjeteta",1);
    t0->SetBranchStatus("subrefjetphi",1);
    // constituent momenta
    t0->SetBranchStatus("refjetcpt",1);
    t0->SetBranchStatus("refjetceta",1);
    t0->SetBranchStatus("refjetcphi",1);

    // Gen
    // number of subjets
    t0->SetBranchStatus("subgenjetnum",1);
    // subjets momenta
    t0->SetBranchStatus("subgenjetpt",1);
    t0->SetBranchStatus("subgenjeteta",1);
    t0->SetBranchStatus("subgenjetphi",1);
    // constituent momenta
    t0->SetBranchStatus("genjetcpt",1);
    t0->SetBranchStatus("genjetceta",1);
    t0->SetBranchStatus("genjetcphi",1);
    
    // SETTING BRANCH ADDRESSES //

    // EVENT BRANCHES
    t0->SetBranchAddress("nref",&nref);
    t0->SetBranchAddress("ngen",&ngen);
    t0->SetBranchAddress("vz",&vz);
    t2->SetBranchAddress("weight",&w);
    t1->SetBranchAddress("pPAprimaryVertexFilter",&ppVF);

    // JET BRANCHES

    // Reco
    t0->SetBranchAddress("jtpt",jtpt);
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);

    // Ref
    t0->SetBranchAddress("refpt",refpt);
    t0->SetBranchAddress("refeta",refeta);
    t0->SetBranchAddress("refphi",refphi);

    // Gen
    t0->SetBranchAddress("genpt",genpt);
    t0->SetBranchAddress("geneta",geneta);
    t0->SetBranchAddress("genphi",genphi);

    // SUBJET BRANCHES

    // Reco
    // number of subjets
    t0->SetBranchAddress("subrecojetnum",subrecojetnum);
    // subjets momenta
    t0->SetBranchAddress("subrecojetpt",subrecojetpt);
    t0->SetBranchAddress("subrecojeteta",subrecojeteta);
    t0->SetBranchAddress("subrecojetphi",subrecojetphi);
    // constituent momenta
    t0->SetBranchAddress("recojetcpt",recojetcpt);
    t0->SetBranchAddress("recojetceta",recojetceta);
    t0->SetBranchAddress("recojetcphi",recojetcphi);

    // Ref
    // number of subjets
    t0->SetBranchAddress("subrefjetnum",subrefjetnum);
    // subjets momenta
    t0->SetBranchAddress("subrefjetpt",subrefjetpt);
    t0->SetBranchAddress("subrefjeteta",subrefjeteta);
    t0->SetBranchAddress("subrefjetphi",subrefjetphi);
    // constituent momenta
    t0->SetBranchAddress("refjetcpt",refjetcpt);
    t0->SetBranchAddress("refjetceta",refjetceta);
    t0->SetBranchAddress("refjetcphi",refjetcphi);

    // Gen
    // number of subjets
    t0->SetBranchAddress("subgenjetnum",subgenjetnum);
    // subjets momenta
    t0->SetBranchAddress("subgenjetpt",subgenjetpt);
    t0->SetBranchAddress("subgenjeteta",subgenjeteta);
    t0->SetBranchAddress("subgenjetphi",subgenjetphi);
    // constituent momenta
    t0->SetBranchAddress("genjetcpt",genjetcpt);
    t0->SetBranchAddress("genjetceta",genjetceta);
    t0->SetBranchAddress("genjetcphi",genjetcphi);

    // OUTPUT FILE //

    // making and pointing to output file
    TFile *f1 = new TFile("subjet_analysis_5_21_1.root", "recreate");
    
    // EVENT PROCESSING //
    
    // EVENT LOOP
    // for(unsigned int i=0; i<t0->GetEntries(); i++){
    for(unsigned int i=0; i<100; i++){

        // // printing event number
        cout<< "event " << i << " is being processed" << endl;

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

                // staging the input file again so I can still Get() everything needed
                fi->cd();

                // SUBJET VISUALIZATION TEST

                double visualize = 0;
                // visualize = 1;
                if(visualize!=1){

                    // making max number of visualized jets, subjets, and jet constituents
                    const int MAX_VIZ_JETS = 10;
                    const int MAX_VIZ_SUBJETS = 30;
                    const int MAX_VIZ_CONSTITUENTS = 200;

                    // making arrays for the visualization tool
                    Double_t jet_mom_[2][3][MAX_VIZ_JETS] = {0};
                    Double_t subjet_mom_[2][3][MAX_VIZ_JETS][MAX_VIZ_SUBJETS] = {0};
                    Double_t constituents_mom_[2][3][MAX_VIZ_JETS][MAX_VIZ_CONSTITUENTS] = {0};

                    double counted = 0;

                    // subjet visualization loop
                    for(unsigned int j=0; j< MAX_VIZ_JETS; j++){

                        // jet and ref pt cut
                        if((jtpt[j]<60)||(refpt[j]<60)){continue;}
                        
                        // skipping index value if jet is not in event
                        if((j>nref)||(j==nref)){continue;}

                        counted += 1;

                        // setting values in jet array momentum to the values in the event
                        
                        // Reco
                        jet_mom_[0][0][j] = jtpt[j];
                        jet_mom_[0][1][j] = jteta[j];
                        jet_mom_[0][2][j] = jtphi[j];
                        
                        // Ref
                        jet_mom_[1][0][j] = refpt[j];
                        jet_mom_[1][1][j] = refeta[j];
                        jet_mom_[1][2][j] = refphi[j];

                        for(unsigned int s=0; s< MAX_VIZ_SUBJETS; s++){
                        
                            // subrecojet and subrefjet pt cut
                            if((subrecojetpt[j][s]<5)||(subrefjetpt[j][s]<5)){continue;}

                            counted += 1;

                            // skipping index value if subjet is not in jet, flagged as bad, or under ptcut

                            // Reco
                            if((s>subrecojetnum[j])||(s==subrecojetnum[j])||(subrecojetpt[j][s]==-999)||(subrecojetpt[j][s]==0)){continue;}

                            // Ref
                            if((s>subrefjetnum[j])||(s==subrefjetnum[j])||(subrefjetpt[j][s]==-999)||(subrefjetpt[j][s]==0)){continue;}

                            // setting values in subjet arrays momentum to the values in the jet

                            // Reco
                            subjet_mom_[0][0][j][s] = subrecojetpt[j][s];
                            subjet_mom_[0][1][j][s] = subrecojeteta[j][s];
                            subjet_mom_[0][2][j][s] = subrecojetphi[j][s];

                            // Ref
                            subjet_mom_[1][0][j][s] = subrefjetpt[j][s];
                            subjet_mom_[1][1][j][s] = subrefjeteta[j][s];
                            subjet_mom_[1][2][j][s] = subrefjetphi[j][s];

                        }

                        for(unsigned int c=0; c< MAX_VIZ_CONSTITUENTS; c++){

                            // adding the momenta information to the vizualization arrays

                            // reco
                            constituents_mom_[0][0][j][c] = recojetcpt[j][c];
                            constituents_mom_[0][1][j][c] = recojetceta[j][c];
                            constituents_mom_[0][2][j][c] = recojetcphi[j][c];

                            // ref
                            constituents_mom_[1][0][j][c] = refjetcpt[j][c];
                            constituents_mom_[1][1][j][c] = refjetceta[j][c];
                            constituents_mom_[1][2][j][c] = refjetcphi[j][c];
                        }
                    }
                    if(counted > 2){
                        f1->cd();
                        // // visualize_subjets1(jet_mom_, subjet_mom_, i);
                        // // visualize_subjets2(jet_mom_, subjet_mom_, constituents_mom_, i);
                        visualize_subjets3(jet_mom_,subjet_mom_,constituents_mom_, i);
                        fi->cd();
                    }
                }
                

                // RECO JET LOOP
                for(unsigned int j=0; j<nref; j++){

                    // jet and ref pt cut
                    if((jtpt[j]<60)||(refpt[j]<60)){continue;}

                    // FILLING 1D HISTS

                    // number of subjets

                    // Reco
                    hsubrecojetnum->Fill(subrecojetnum[j], w);

                    // Ref
                    hsubrefjetnum->Fill(subrefjetnum[j], w);

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
                    hsubrecojetnum_jtpt->Fill(jtpt[j], subrecojetnum[j], w);

                    // number of subjets reco vs ref
                    hsubrecojetnum_subrefjetnum->Fill(subrefjetnum[j], subrecojetnum[j], w);

                    // making variable for number of jets with pt cut

                    // number of ptcuts
                    const Int_t numptcut = 50;
                    
                    // Reco
                    double subrecojetnum_ptcut5 = 0;
                    double subrecojetnums[numptcut] = {0};
                    
                    // Ref
                    double subrefjetnum_ptcut5 = 0;
                    double subrefjetnums[numptcut] = {0};

                    // RECO SUBJET LOOP
                    for(unsigned int s=0; s<subrecojetnum[j]; s++){
                        
                        // getting number of subjets under various pt cuts
                        for(unsigned int ct = numh1d1[1]; ct < numh1d1[2]; ct++){
                            if(subrecojetpt[j][s] > ct){subrecojetnums[ct] += 1;}
                        }

                        // Reco subjet pt cut
                        if(subrecojetpt[j][s]>5){
                            
                            // adding one to number of subjets with pt > 5 GeV
                            subrecojetnum_ptcut5 +=1;

                            // filling subrecojet momenta hists
                            hsubrecojetpt->Fill(subrecojetpt[j][s], w);
                            hsubrecojeteta->Fill(subrecojeteta[j][s], w);
                            hsubrecojetphi->Fill(subrecojetphi[j][s], w);

                        }
                    }

                    // REF SUBJET LOOP
                    for(unsigned int s=0; s<subrefjetnum[j]; s++){
                        
                        // getting number of subjets under various pt cuts
                        for(unsigned int ct = numh1d1[1]; ct < numh1d1[2]; ct++){
                            if(subrefjetpt[j][s] > ct){subrefjetnums[ct] += 1;}
                        }

                        // Ref subjet pt cut
                        if(subrefjetpt[j][s]>5){
                            
                            // adding one to number of subjets with pt > 5 GeV
                            subrefjetnum_ptcut5 +=1;

                            // filling subrecojet momenta hists
                            hsubrefjetpt->Fill(subrefjetpt[j][s], w);
                            hsubrefjeteta->Fill(subrefjeteta[j][s], w);
                            hsubrefjetphi->Fill(subrefjetphi[j][s], w);
                        }
                    }

                    // adding the number of subjets with pt cuts to the hist

                    // Reco
                    hsubrecojetnum_ptcut5->Fill(subrecojetnum_ptcut5, w);
                    for(unsigned int ct = numh1d1[1]; ct < numh1d1[2]; ct++){
                        hsubrecojetnum_ptcut->Fill(ct, subrecojetnums[ct], w);
                    }

                    // Ref
                    hsubrefjetnum_ptcut5->Fill(subrefjetnum_ptcut5, w);
                    for(unsigned int ct = numh1d1[1]; ct < numh1d1[2]; ct++){
                        hsubrefjetnum_ptcut->Fill(ct, subrefjetnums[ct], w);
                    }

                    // Reco vs Ref
                    hsubrecojetnum_subrefjetnum_ptcut5->Fill(subrefjetnum_ptcut5, subrecojetnum_ptcut5, w);

                }

                // visualize_subjets1(jet_mom_, subjet_mom_, heta_phi, blah, i);

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
                    hsubgenjetnum->Fill(subgenjetnum[g], w);

                    // number of subjets in gen jet with pt > 5 GeV
                    Double_t subgenjetnum_ptcut5 = 0;

                    // GEN SUBJET LOOP
                    for(unsigned int s=0; s< subgenjetnum[g]; s++){

                        // Gen 
                        if(subgenjetpt[g][s]>5){

                            // adding one to number of subgenjets with pt > 5 GeV
                            subgenjetnum_ptcut5+=1;

                            // filling gen subjet momenta hists
                            hsubgenjetpt->Fill(subgenjetpt[g][s], w);
                            hsubgenjeteta->Fill(subgenjeteta[g][s], w);
                            hsubgenjetphi->Fill(subgenjetphi[g][s], w);
                        }
                    }

                    // filling gen subjet number with pt > 5 GeV
                    hsubgenjetnum_ptcut5->Fill(subgenjetnum_ptcut5, w);
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
    normalizeh(hngen);
    normalizeh(hnref);

    // JET HISTS

    // Reco
    normalizeh(hjtpt);
    normalizeh(hjteta);
    normalizeh(hjtphi);

    // Ref
    normalizeh(hrefpt);
    normalizeh(hrefeta);
    normalizeh(hrefphi);

    // Gen
    normalizeh(hgenpt);
    normalizeh(hgeneta);
    normalizeh(hgenphi);

    // SUBJET HISTS

    // Reco
    // number of subjets
    normalizeh(hsubrecojetnum);
    normalizeh(hsubrecojetnum_ptcut5);
    // subjet momenta
    normalizeh(hsubrecojetpt);
    normalizeh(hsubrecojeteta);
    normalizeh(hsubrecojetphi);

    // Ref
    // number of subjets
    normalizeh(hsubrefjetnum);
    normalizeh(hsubrefjetnum_ptcut5);
    // subjet momenta
    normalizeh(hsubrefjetpt);
    normalizeh(hsubrefjeteta);
    normalizeh(hsubrefjetphi);

    // Gen
    // number of subjets
    normalizeh(hsubgenjetnum);
    normalizeh(hsubgenjetnum_ptcut5);
    // subjet momenta
    normalizeh(hsubgenjetpt);
    normalizeh(hsubgenjeteta);
    normalizeh(hsubgenjetphi);

    // WRITING HISTOGRAMS //

    // EVENT HISTS
    hvz->Write();
    hngen->Write();
    hnref->Write();

    // JET HISTS

    // Reco
    hjtpt->Write();
    hjteta->Write();
    hjtphi->Write();

    // Ref
    hrefpt->Write();
    hrefeta->Write();
    hrefphi->Write();

    // Gen
    hgenpt->Write();
    hgeneta->Write();
    hgenphi->Write();

    double savepngs = 0;
    savepngs = 1;
    if(savepngs != 1){

        // SUBJET 1D HISTS
        
        // Reco
        // number of subjets
        save_h1d_1(hsubrecojetnum, "Number of Subjets in Reco Jets", "Probability", "hsubrecojetnum");
        save_h1d_1(hsubrecojetnum_ptcut5, "Number of Subjets in Reco Jets", "Probability", "hsubrecojetnum_ptcut5");
        // subjet momenta
        save_h1d_1(hsubrecojetpt, "pt", "Probability", "subrecojetpt");
        save_h1d_1(hsubrecojeteta, "#eta", "Probability", "subrecojeteta");
        save_h1d_1(hsubrecojetphi, "#phi", "Probability", "subrecojetphi");
        
        // Ref
        // number of subjets
        save_h1d_1(hsubrefjetnum, "Number of Subjets in Reco Jets", "Probability", "hsubrefjetnum");
        save_h1d_1(hsubrefjetnum_ptcut5, "Number of Subjets in Reco Jets", "Probability", "hsubrefjetnum_ptcut5");
        // subjet momenta
        save_h1d_1(hsubrefjetpt, "pt", "Probability", "subrefjetpt");
        save_h1d_1(hsubrefjeteta, "#eta", "Probability", "subrefjeteta");
        save_h1d_1(hsubrefjetphi, "#phi", "Probability", "subrefjetphi");
        
        // Gen
        // number of subjets
        save_h1d_1(hsubgenjetnum, "Number of Subjets in Reco Jets", "Probability", "hsubgenjetnum");
        save_h1d_1(hsubgenjetnum_ptcut5, "Number of Subjets in Reco Jets", "Probability", "hsubgenjetnum_ptcut5");
        // subjet momenta
        save_h1d_1(hsubgenjetpt, "pt", "Probability", "subgenjetpt");
        save_h1d_1(hsubgenjeteta, "#eta", "Probability", "subgenjeteta");
        save_h1d_1(hsubgenjetphi, "#phi", "Probability", "subgenjetphi");

        // OVERLAY 1D HISTS

        // number of subjets for ref and reco, with and without pt cuts
        save_h1d_2(hsubrecojetnum, hsubrefjetnum, "Number of Subjets", "Probability", "subjetnums", "reco", "gen");
        save_h1d_2(hsubrecojetnum_ptcut5, hsubrefjetnum_ptcut5, "Number of Subjets with pt > 5 GeV", "Probability", "subjetnums_ptcut5", "reco", "gen");

        // subjet momenta for gen and reco
        save_h1d_2(hsubrecojetpt, hsubgenjetpt, "pt", "", "reco_gen_subjet_pt", "p_{T}^{subjet} Reco", "p_{T}^{subjet} Gen");
        save_h1d_2(hsubrecojeteta, hsubgenjeteta, "#eta", "", "reco_gen_subjet_eta", "#eta^{subjet} Reco", "#eta^{subjet} Gen");
        save_h1d_2(hsubrecojetphi, hsubgenjetphi, "#phi", "", "reco_gen_subjet_phi", "#phi^{subjet} Reco", "#phi^{subjet} Gen");

        // SUBJET 2D HISTS

        // number of gen vs reco subjets in various jet pt regions
        save_h2d_1(hsubrecojetnum_subrefjetnum, "subrefjetnum", "subrecojetnum", "hsubrecojetnum_subrefjetnum");
        save_h2d_1(hsubrecojetnum_subrefjetnum_ptcut5, "subrefjetnum", "subrecojetnum", "hsubrecojetnum_subrefjetnum_ptcut5");
        
        // number of subjets vs subjet pt cut
        save_h2d_1(hsubrecojetnum_ptcut, "p_{T}^{subjet} cut in GeV", "number of subjets", "hsubrecojetnum_ptcut");
        save_h2d_1(hsubrefjetnum_ptcut, "p_{T}^{subjet} cut in GeV", "number of subjets", "hsubrefjetnum_ptcut");

        // // reco and gen jet pt vs subjet number
        // save_h2d_1(hsubjetnum1_jtpt, "reco subjetnum1", "reco jet pt", "subjetnum1_jtpt");
        // save_h2d_1(hsubjetgennum1_genpt, "gen subjetnum1", "gen jet pt", "subjetgennum1_genpt");
    }
    // // PRINTING COUNTS OF INTEREST
    // cout <<endl<< "the number of subjets looked at was " << subjetcount <<endl;
    // cout << "the number of subjets with pt > 1 TeV looked at was " << subjetcount_hugept <<endl<<endl;

}
