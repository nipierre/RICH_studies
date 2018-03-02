#include <fstream>
#include <iomanip> 
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

#include <TROOT.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <RooGlobalFunc.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooAddition.h>
#include <RooCategory.h>
#include <RooCBShape.h>
#include <RooChebychev.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooMinuit.h>
#include <RooMy.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooRelBreitWigner.h>
#include <RooSimultaneous.h>
#include <RooVoigtian.h>


/*
 * 
 *	const TMatrixDSym& cor = r->correlationMatrix() ;
 *	const TMatrixDSym& cov = r->covarianceMatrix() ;
 *	// Print correlation, covariance matrix
 *	cor.Print() ;
 *	cov.Print() ;
 * 
 * 
 *	// Make list of model parameters
 *	RooArgSet* params = model.getParameters(x) ;
 *	// Write LaTex table to file
 *	params->printLatex(Sibling(*initParams),OutputFile("rf407_latextables.tex")) ;
 * 
 * 
 * 
 */



// ******************************************************************************************

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
using std::setw;
using std::setprecision;
using std::ifstream;


using namespace RooFit ;

// ******************************************************************************************
// const int Np = 6;
// const int Nt = 7;
// const double p_bins[Np+1]={10. ,15. ,20. ,25. ,30. ,40. ,50.};
// const double t_bins[Nt+1]={0.00,0.01,0.02,0.03,0.04,0.06,0.09,0.12};

const int Np = 13;
const int Nt = 4;
const double p_bins[Np+1]={10.,11.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.,50.};
const double t_bins[Nt+1]={0.00,0.01,0.04,0.12,0.3};
// const double t_bins[Nt+1]={0.00,0.01,0.12};

// const int Np = 1;
// const int Nt = 1;
// const double p_bins[Np+1]={10.,50.};
// const double t_bins[Nt+1]={0.01,0.12};

double N_id[8][5][Np][Nt];
//double R_id[6][4][Np][Nt];

double last_para[30];
double last_para2[Np][Nt][30];

string data_file;
string hist_file = "hist.root";
string hist_file_k0 = "hist.root";
string hist_file_lam = "hist.root";
string hist_file_phi = "hist.root";
string hist_file_rho = "hist.root";
string out_file = "rich.root";
double lh_cut[4][6];
TH1D* h[8][7][Np][Nt];
TH2D* h2[8][5][Np][Nt];
RooFitResult *r[8][Np][Nt];
int lw = 1;
double thr_diff = 0.;
int retry = 20;
bool rpipe = false;

double rho_z1 = 0.;
double rho_z2 = 1.;

bool use_improve = false;
bool use_hesse = true;
bool use_minos = false;
bool use_sidebins = true;
bool first_phi = true;
TFile* input;
TFile* input_k0;
TFile* input_rho;
TFile* input_phi;
TFile* input_lam;

// ******************************************************************************************

int main(int, char**);
void read_options();
void get_input_data_t1();
void get_input_data_t2();
void get_input_data_rho();
void get_input_data2();

RooDataHist* gen_k0(int, int , int, int);
void write_hist();
void create_hist();
void get_plots();
void rm_plots();
void fit_table_k0(int);
void fit_table_phi(int);
void fit_table_lambda(int);
void fit_table_rho(int);
void print_table();
void set_plot_style();
