#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include <string>

#include <TROOT.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TStyle.h>
#include <TPaletteAxis.h>
#include <TTree.h>


/**********************************************************************/

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::stringstream;
using std::vector;
using std::pair;
using std::string;
using std::make_pair;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::setw;

/**********************************************************************/

int main(int, char**);
void draw1();
void draw2();
void draw3();
void draw4();
void set_plot_style();

/**********************************************************************/
int main(int argc, char *argv[]){
	TH1F * tt = new TH1F("t","",10,0,10);
	TTree *ttt = new TTree();
	TLorentzVector* tttt = new TLorentzVector();
	delete tt;
	delete ttt;
	delete tttt;
	
	// draw1();
	// draw2();
	// draw3();
	draw4();
	return 0;
}

/**********************************************************************/
void draw1(){
	TFile * f = new TFile("eff_2007_matrix_1.00.root","read");

	stringstream nn;
	
// 	hpos0 = eff (pi+->pi+)
// 	hpos1 =      K+->pi+
// 	hpos2 =      pi+->K+
// 	hpos3 =      K+->K+
	
	string name[2] = {"hpos","hneg"};

	TH2D* h2;
	vector<double> mom[2][4][7];
	vector<double> eff[2][4][7];
	vector<double> eff_r[2][4][7];
	double th[8];
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<4; i++){  
			nn.str("");
			nn.clear();
			nn << name[c] << i;
			
			h2= (TH2D*)f->Get(nn.str().c_str());
			for(int bx = 1; bx <= h2->GetXaxis()->GetNbins(); bx++){
				for(int by = 1; by <= h2->GetYaxis()->GetNbins(); by++){
					
					if( h2->GetXaxis()->GetBinCenter(bx) < 10. || h2->GetXaxis()->GetBinCenter(bx) > 50.0 ) continue;
					th[by-1] = h2->GetYaxis()->GetBinLowEdge(by) ;
					th[by] = h2->GetYaxis()->GetBinUpEdge(by) ;
					mom[c][i][by-1].push_back(h2->GetXaxis()->GetBinCenter(bx));
					eff[c][i][by-1].push_back(h2->GetBinContent(bx,by));
					eff_r[c][i][by-1].push_back(h2->GetBinError(bx,by));
					
				}
			}
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","",800,600);
	set_plot_style();
	int col[14];
	TH2F* palpal = new TH2F("palpal","",2,0,1,1,0,1);
	
	palpal->Fill(0.005,0.1);
	
	TLegend* leg= new TLegend(0.15,0.15,0.35,0.49);
	leg->SetLineColor(0);
	leg->SetNColumns(1);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.039);
	
	
	
	palpal->Draw("colz");
	gPad->Update();
	
	gPad->Update();
	TPaletteAxis *palette = (TPaletteAxis*)palpal->GetListOfFunctions()->FindObject("palette");
	for(int i = 0; i< 7;i++) col[i] = palette->GetValueColor(i/6.3);
	
	TGraphErrors* gr[2][4][7];
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<4; i++){  // pi, K , p
			for(int k = 0; k<7;k++){
				// gr[c][i][k] = new TGraphErrors(mom[c][i][k].size(),&mom[c][i][k][0],&eff[c][i][k][0],0,&eff_r[c][i][k][0]);
				gr[c][i][k] = new TGraphErrors(mom[c][i][k].size(),&mom[c][i][k][0],&eff[c][i][k][0],0,0);
				gr[c][i][k]->SetMarkerColor(col[k]);
				gr[c][i][k]->SetLineColor(col[k]);
				gr[c][i][k]->SetMarkerStyle(33);
				gr[c][i][k]->SetMarkerSize(1.5);
				nn.str("");
				nn.clear();
				nn<< th[k] << " < #theta < " << th[k+1];
				if(c == 0 && i ==0)leg->AddEntry(gr[c][i][k], nn.str().c_str(),	"p");
				
			}
		}
	}
	
	TH1F* hr;
	
	string charge[2] = {"p","m"};
	string charge2[2] = {"+","-"};
	string is[4] = {"pi","k","pi","k"};
	string id[4] = {"pi","pi","k","k"};
	string is2[4] = {"#pi","k","#pi","k"};
	string id2[4] = {"#pi","#pi","k","k"};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<4; i++){  // pi, K , p
			
			c1->Clear();
			hr = c1->DrawFrame(10,0,50,1);
			hr->GetXaxis()->SetTitleOffset(0.9);
			hr->GetXaxis()->SetTitleSize(0.05);
			hr->GetXaxis()->SetLabelSize(0.05);
			hr->GetXaxis()->SetTitle("p (GeV/#it{c})");
			hr->GetXaxis()->SetNdivisions(504);
			
			hr->GetYaxis()->SetTitleOffset(0.9);
			hr->GetYaxis()->SetTitleSize(0.05);
			hr->GetYaxis()->SetLabelSize(0.05);
			hr->GetYaxis()->SetNdivisions(504);
						
			nn.str("");
			nn << "Eff("<<is2[i] << "^{" <<charge2[c] <<"} #rightarrow " <<id2[i] <<")";
			hr->GetYaxis()->SetTitle(nn.str().c_str());
			leg->Draw();
			for(int k=0; k<7;k++){
				gr[c][i][k]->Draw("lp");
			}
			nn.str("");
			nn << "rich_fig/2007/plot_07_" << is[i] << charge[c] << "_" << id[i] << ".pdf";
			c1->Print(nn.str().c_str());
			delete hr;
		}
	}
	
	
	delete c1;
}

/**********************************************************************/
void draw2(){
	TFile * f = new TFile("rich_2011.root","read");

	stringstream nn;
	
	// string name[2] = {"hpos","hneg"};

	string charge[2] = {"p","m"};
	string charge2[2] = {"+","-"};
	string is[3] = {"pi","k","p"};
	string id[4] = {"pi","k","p","u"};
	string is2[3] = {"#pi","k","p"};
	string id2[4] = {"#pi","k","p","unk"};
	
	TGraphErrors* gr[2][3][4][4];
	double th[5] = {0,0.01,0.04,0.12,0.3};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 1;t<3;t++){
					nn.str("");
					nn.clear();
					nn << is[i] << "_" << charge[c] << "_" << id[j]  << "_" << t;
					cout << nn.str().c_str() << endl;
					gr[c][i][j][t]= (TGraphErrors*)f->Get(nn.str().c_str());
					
				}
			}
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","",800,600);
	set_plot_style();
	int col[4];
	TH2F* palpal = new TH2F("palpal","",2,0,1,1,0,1);
	
	palpal->Fill(0.005,0.1);
	
	TLegend* leg= new TLegend(0.15,0.15,0.35,0.49);
	leg->SetLineColor(0);
	leg->SetNColumns(1);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.039);
	
	
	
	palpal->Draw("colz");
	gPad->Update();
	
	gPad->Update();
	TPaletteAxis *palette = (TPaletteAxis*)palpal->GetListOfFunctions()->FindObject("palette");
	for(int i = 0; i< 4;i++) col[i] = palette->GetValueColor(i/4.3);
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 1;t<3;t++){
					gr[c][i][j][t]->SetMarkerColor(col[t]);
					gr[c][i][j][t]->SetLineColor(col[t]);
					gr[c][i][j][t]->SetMarkerStyle(33);
					gr[c][i][j][t]->SetMarkerSize(1.5);
					nn.str("");
					nn.clear();
					nn<< th[t] << " < #theta < " << th[t+1];
					if(c == 0 && i ==0 && j == 0 )leg->AddEntry(gr[c][i][j][t], nn.str().c_str(),	"p");
				}
			}
		}
	}
	
	TH1F* hr;
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				c1->Clear();
				hr = c1->DrawFrame(10,0,50,1);
				hr->GetXaxis()->SetTitleOffset(0.9);
				hr->GetXaxis()->SetTitleSize(0.05);
				hr->GetXaxis()->SetLabelSize(0.05);
				hr->GetXaxis()->SetTitle("p (GeV/#it{c})");
				hr->GetXaxis()->SetNdivisions(504);
				
				hr->GetYaxis()->SetTitleOffset(0.9);
				hr->GetYaxis()->SetTitleSize(0.05);
				hr->GetYaxis()->SetLabelSize(0.05);
				hr->GetYaxis()->SetNdivisions(504);
							
				nn.str("");
				nn << "Eff("<<is2[i] << "^{" <<charge2[c] <<"} #rightarrow " <<id2[j] <<")";
				hr->GetYaxis()->SetTitle(nn.str().c_str());
				leg->Draw();
				for(int t=1; t<3;t++){
					gr[c][i][j][t]->Draw("lp");
				}
				nn.str("");
				nn << "rich_fig/2011/plot_11_" << is[i] << charge[c] << "_" << id[j] << ".pdf";
				c1->Print(nn.str().c_str());
				delete hr;
			}
		}
	}
	
	
	delete c1;
}

/**********************************************************************/
void draw3(){
	TFile * f = new TFile("rich_2011_more_theta_2.root","read");

	stringstream nn;
	
	// string name[2] = {"hpos","hneg"};

	string charge[2] = {"p","m"};
	string charge2[2] = {"+","-"};
	string is[3] = {"pi","k","p"};
	string id[4] = {"pi","k","p","u"};
	string is2[3] = {"#pi","k","p"};
	string id2[4] = {"#pi","k","p","noID"};
	
	TGraphErrors* gr[2][3][4][4];
	double th[5] = {0,0.01,0.04,0.12,0.3};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 0;t<4;t++){
					nn.str("");
					nn.clear();
					nn << is[i] << "_" << charge[c] << "_" << id[j]  << "_" << t;
					cout << nn.str().c_str() << endl;
					gr[c][i][j][t]= (TGraphErrors*)f->Get(nn.str().c_str());
					
				}
			}
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","",800,600);
	set_plot_style();
	int col[4];
	TH2F* palpal = new TH2F("palpal","",2,0,1,1,0,1);
	
	palpal->Fill(0.005,0.1);
	
	TLegend* leg= new TLegend(0.15,0.15,0.35,0.49);
	leg->SetLineColor(0);
	leg->SetNColumns(1);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.039);
	
	
	
	palpal->Draw("colz");
	gPad->Update();
	
	gPad->Update();
	TPaletteAxis *palette = (TPaletteAxis*)palpal->GetListOfFunctions()->FindObject("palette");
	// for(int i = 0; i< 4;i++) col[i] = palette->GetValueColor(i/4.3);
	for(int i = 0; i< 4;i++) col[i] = palette->GetValueColor(i/4.1);
	
	double min[3][4] = {
		{0.5,0.,0.,0.},
		{0.,0.,0.,0.},
		{0.,0.,0.,0.}
	};
	
	double max[3][4] = {
		{1,0.3,0.3,0.3},
		{0.5,1.,0.5,0.5},
		{0.65,0.3,1.,0.3}
	};
	
	double leg_x[3][4] = {
		{0.15,0.15,0.15,0.15},
		{0.40,0.40,0.40,0.40},
		{0.60,0.60,0.60,0.60}
	};
	
	double leg_y[3][4] = {
		{0.15,0.5,0.5,0.5},
		{0.55,0.15,0.55,0.55},
		{0.55,0.55,0.15,0.55},
	};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 0;t<4;t++){
					for(int p=0;p<gr[c][i][j][t]->GetN();p++){
						double px,py;
						gr[c][i][j][t]->GetPoint(p,px,py);
						px += 0.2*t;
						gr[c][i][j][t]->SetPoint(p,px,py);
					}
					gr[c][i][j][t]->SetMarkerColor(col[t]);
					gr[c][i][j][t]->SetLineColor(col[t]);
					gr[c][i][j][t]->SetMarkerStyle(33);
					gr[c][i][j][t]->SetMarkerSize(1.5);
					nn.str("");
					nn.clear();
					nn<< th[t] << " < #theta < " << th[t+1];
					if(c == 0 && i ==0 && j == 0 )leg->AddEntry(gr[c][i][j][t], nn.str().c_str(),	"p");
				}
			}
		}
	}
	
	TH1F* hr;
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.05);
	lat->SetTextColor(kGray+2);
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				c1->Clear();
				hr = c1->DrawFrame(10,min[i][j],50,max[i][j]);
				
				leg->SetX1NDC(leg_x[i][j]);
				leg->SetX2NDC(leg_x[i][j]+0.2);
				leg->SetY1NDC(leg_y[i][j]);
				leg->SetY2NDC(leg_y[i][j]+0.3);
				hr->GetXaxis()->SetTitleOffset(0.9);
				hr->GetXaxis()->SetTitleSize(0.05);
				hr->GetXaxis()->SetLabelSize(0.05);
				hr->GetXaxis()->SetTitle("p (GeV/#it{c})");
				hr->GetXaxis()->SetNdivisions(504);
				
				hr->GetYaxis()->SetTitleOffset(0.9);
				hr->GetYaxis()->SetTitleSize(0.05);
				hr->GetYaxis()->SetLabelSize(0.05);
				hr->GetYaxis()->SetNdivisions(504);
							
				nn.str("");
				if(is[i] != "p")nn << "#epsilon("<<is2[i] << "^{" <<charge2[c] <<"} #rightarrow " <<id2[j] <<")";
				else{
					if(charge2[c] == "+")nn << "#epsilon("<<is2[i] << " #rightarrow " <<id2[j] <<")";
					if(charge2[c] == "-")nn << "#epsilon(#bar{"<<is2[i] << "} #rightarrow " <<id2[j] <<")";
				}
				hr->GetYaxis()->SetTitle(nn.str().c_str());
				leg->Draw();
				// lat->DrawLatex(leg_x[i][j]+0.02,leg_y[i][j]+0.33,"Preliminary");
				for(int t=0; t<4;t++){
					gr[c][i][j][t]->Draw("lp");
				}
				nn.str("");
				nn << "rich_fig/2011_2/plot_11_2_" << is[i] << charge[c] << "_" << id[j] << ".pdf";
				c1->Print(nn.str().c_str());
				delete hr;
			}
		}
	}
	
	
	delete c1;
}

/**********************************************************************/
void draw4(){
	TFile * f = new TFile("rich_2011_lh07.root","read");

	stringstream nn;
	
	// string name[2] = {"hpos","hneg"};

	string charge[2] = {"p","m"};
	string charge2[2] = {"+","-"};
	string is[3] = {"pi","k","p"};
	string id[4] = {"pi","k","p","u"};
	string is2[3] = {"#pi","k","p"};
	string id2[4] = {"#pi","k","p","noID"};
	
	TGraphErrors* gr[2][3][4][4];
	double th[5] = {0,0.01,0.04,0.12,0.3};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 0;t<4;t++){
					nn.str("");
					nn.clear();
					nn << is[i] << "_" << charge[c] << "_" << id[j]  << "_" << t;
					cout << nn.str().c_str() << endl;
					gr[c][i][j][t]= (TGraphErrors*)f->Get(nn.str().c_str());
					
				}
			}
		}
	}
	
	TCanvas* c1 = new TCanvas("c1","",800,600);
	set_plot_style();
	int col[4];
	TH2F* palpal = new TH2F("palpal","",2,0,1,1,0,1);
	
	palpal->Fill(0.005,0.1);
	
	TLegend* leg= new TLegend(0.15,0.15,0.35,0.49);
	leg->SetLineColor(0);
	leg->SetNColumns(1);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.039);
	
	
	
	palpal->Draw("colz");
	gPad->Update();
	
	gPad->Update();
	TPaletteAxis *palette = (TPaletteAxis*)palpal->GetListOfFunctions()->FindObject("palette");
	// for(int i = 0; i< 4;i++) col[i] = palette->GetValueColor(i/4.3);
	for(int i = 0; i< 4;i++) col[i] = palette->GetValueColor(i/4.1);
	
	double min[3][4] = {
		{0.5,0.,0.,0.},
		{0.,0.,0.,0.},
		{0.,0.,0.,0.}
	};
	
	double max[3][4] = {
		{1,0.4,0.3,0.1},
		{0.5,1.,0.5,0.15},
		{0.65,0.3,1.,0.1}
	};
	
	double leg_x[3][4] = {
		{0.15,0.15,0.15,0.15},
		{0.40,0.40,0.40,0.40},
		{0.60,0.60,0.60,0.60}
	};
	
	double leg_y[3][4] = {
		{0.15,0.5,0.5,0.5},
		{0.55,0.15,0.55,0.55},
		{0.55,0.55,0.15,0.55},
	};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 0;t<4;t++){
					for(int p=0;p<gr[c][i][j][t]->GetN();p++){
						double px,py;
						gr[c][i][j][t]->GetPoint(p,px,py);
						px += 0.2*t;
						gr[c][i][j][t]->SetPoint(p,px,py);
					}
					gr[c][i][j][t]->SetMarkerColor(col[t]);
					gr[c][i][j][t]->SetLineColor(col[t]);
					gr[c][i][j][t]->SetMarkerStyle(33);
					gr[c][i][j][t]->SetMarkerSize(1.5);
					nn.str("");
					nn.clear();
					nn<< th[t] << " < #theta < " << th[t+1];
					if(c == 0 && i ==0 && j == 0 )leg->AddEntry(gr[c][i][j][t], nn.str().c_str(),	"p");
				}
			}
		}
	}
	
	TH1F* hr;
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.05);
	lat->SetTextColor(kGray+2);
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				c1->Clear();
				hr = c1->DrawFrame(10,min[i][j],50,max[i][j]);
				
				leg->SetX1NDC(leg_x[i][j]);
				leg->SetX2NDC(leg_x[i][j]+0.2);
				leg->SetY1NDC(leg_y[i][j]);
				leg->SetY2NDC(leg_y[i][j]+0.3);
				hr->GetXaxis()->SetTitleOffset(0.9);
				hr->GetXaxis()->SetTitleSize(0.05);
				hr->GetXaxis()->SetLabelSize(0.05);
				hr->GetXaxis()->SetTitle("p (GeV/#it{c})");
				hr->GetXaxis()->SetNdivisions(504);
				
				hr->GetYaxis()->SetTitleOffset(0.9);
				hr->GetYaxis()->SetTitleSize(0.05);
				hr->GetYaxis()->SetLabelSize(0.05);
				hr->GetYaxis()->SetNdivisions(504);
							
				nn.str("");
				if(is[i] != "p")nn << "#epsilon("<<is2[i] << "^{" <<charge2[c] <<"} #rightarrow " <<id2[j] <<")";
				else{
					if(charge2[c] == "+")nn << "#epsilon("<<is2[i] << " #rightarrow " <<id2[j] <<")";
					if(charge2[c] == "-")nn << "#epsilon(#bar{"<<is2[i] << "} #rightarrow " <<id2[j] <<")";
				}
				hr->GetYaxis()->SetTitle(nn.str().c_str());
				leg->Draw();
				// lat->DrawLatex(leg_x[i][j]+0.02,leg_y[i][j]+0.33,"Preliminary");
				for(int t=0; t<4;t++){
					gr[c][i][j][t]->Draw("lp");
				}
				nn.str("");
				nn << "rich_fig/2011_lh07/plot_11_lh07_" << is[i] << charge[c] << "_" << id[j] << ".pdf";
				c1->Print(nn.str().c_str());
				delete hr;
			}
		}
	}
	
	
	delete c1;
}

/**********************************************************************/
void set_plot_style(){
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}
/**********************************************************************/
