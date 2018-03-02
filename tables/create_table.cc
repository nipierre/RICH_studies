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
#include <TF1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::stringstream;
using std::vector;
using std::pair;
using std::string;
using std::make_pair;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::fixed;


/***************************************************************/
const int Np = 13;
const int Nt = 4;
const double p_bins[Np+1]={10.,11.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.,50.};
const double t_bins[Nt+1]={0.00,0.01,0.04,0.12,0.3};

vector<double> eff[2][3][4][4]; //charge, Is, id_as, theta
vector<double> eff_err[2][3][4][4]; //charge, Is, id_as, theta

/***************************************************************/
int main(int,char**);
void read_file(string);
void write_output(string);
void clear_vec();

/***************************************************************/
int main(int argc, char *argv[]){
	
	read_file("rich_2011_more_theta_2.root");
	write_output("rich_table_strict.tex");
	
	clear_vec();
	
	read_file("rich_2011_lh07.root");
	write_output("rich_table_lh07.tex");
	
	clear_vec();
	return 0;
}

/***************************************************************/
void read_file(string file){
	TFile * f = new TFile(file.c_str(),"read");
	
	if(!f->IsOpen()) std::exit(1);

	stringstream nn;
	
	string charge[2] = {"p","m"};
	// string charge2[2] = {"+","-"};
	string is[3] = {"pi","k","p"};
	string id[4] = {"pi","k","p","u"};
	// string is2[3] = {"#pi","k","p"};
	// string id2[4] = {"#pi","k","p","noID"};
	
	TGraphErrors* gr;
	
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			for(int j = 0; j<4; j++){  //pi, k ,p, u
				for(int t = 0;t<4;t++){
					nn.str("");
					nn.clear();
					nn << is[i] << "_" << charge[c] << "_" << id[j]  << "_" << t;
					gr= (TGraphErrors*)f->Get(nn.str().c_str());
					
					for(int p=0; p<gr->GetN();p++){
						double t1,t2;
						gr->GetPoint(p, t1, t2);
						eff[c][i][j][t].push_back(t2);
						t1=gr->GetErrorY(p);
						eff_err[c][i][j][t].push_back(t1);
					}
					
					delete gr;
				}
			}
		}
	}
}

/***************************************************************/
void write_output(string file){
	ofstream out;
	out.open(file.c_str());
	out << std::left;
	
	string charge[2] = {"p","m"};
	string charge2[2] = {"+","-"};
	string is[3] = {"pi","k","p"};
	// string id[4] = {"pi","k","p","u"};
	string is2[3] = {"\\pi","k","p"};
	// string id2[4] = {"#pi","k","p","noID"};
	
	for(int c = 0; c<2; c++){ // + , -
		for(int i = 0; i<3; i++){  //pi, k ,p
			
			out << "\\begin{table}[htbp]" << endl
				<< "\\centering" << endl
				<< "\\caption{RICH efficiency tables for $" << is2[i] <<"^" << charge2[c] << "$.}" << endl
				<< "\\small" << endl
				<< "\\begin{tabular}{lcccc}" << endl
				<< "\\toprule" << endl
				<< "\\tableheadline{Range} & "
				<< "\\tableheadline{$\\epsilon(" << is2[i] << "^" <<charge2[c] << "\\rightarrow \\pi^" << charge2[c] << ")$} &"
				<< "\\tableheadline{$\\epsilon(" << is2[i] << "^" <<charge2[c] << "\\rightarrow k^" << charge2[c] << ")$} &"; 
			if(c==0){
				out	<< "\\tableheadline{$\\epsilon(" << is2[i] << "^" <<charge2[c] << "\\rightarrow p)$} &";
			} else {
				out	<< "\\tableheadline{$\\epsilon(" << is2[i] << "^" <<charge2[c] << "\\rightarrow \\bar{p})$} &";
			}
			out	<< "\\tableheadline{$\\epsilon(" << is2[i] << "^" <<charge2[c] << "\\rightarrow \\text{no ID})$} \\\\"  << endl;
			
			for(int t = 0;t<4;t++){
				out << "\\midrule" << endl;
				out << fixed << setprecision(2)<< "\\multicolumn{5}{l}{$" << t_bins[t] << " \\leq \\theta < " << t_bins[t+1] << "$}\\\\"<<endl;
				out << "\\midrule" << endl;
				for(unsigned int p=0; p<eff_err[c][i][0][t].size(); p++){
					
					out << fixed << setprecision(0) << "$" << p_bins[p] << "\\GeV/c^2 \\leq p < " << p_bins[p+1] << "\\GeV/c^2$ & ";
					for(int j = 0; j<4; j++){  //pi, k ,p, u
						int prec = 0;
						for(int l=5; l>=0; l--){
							if(eff_err[c][i][j][t][p] *TMath::Power(10,l) > 1) prec = l;
							
						}
						prec++;
						
						out << "$" 
							<< setw(6) << fixed << setprecision(prec) << eff[c][i][j][t][p]
							<< " \\pm " 
							<< setw(6) << fixed << setprecision(prec) << eff_err[c][i][j][t][p]
							<< "$ ";
						if(j!=3) out << "&"; 
					}
					out << "\\\\" << endl;
				}
			}
			out << "\\bottomrule" << endl
				<< "\\end{tabular}" << endl
				<< "\\label{tab:prob_" << is[i]<<charge[c] << "}" << endl
				<< "\\end{table}" << endl << endl << endl;
		}
	}
	out.close();
}
 
/***************************************************************/
void clear_vec(){
	for(int c=0; c<2; c++){
		for(int i=0; i<3; i++){
			for(int j=0; j<4; j++){
				for(int t=0; t<4; t++){
					eff[c][i][j][t].clear();
					vector<double>().swap(eff[c][i][j][t]);
					
					eff_err[c][i][j][t].clear();
					vector<double>().swap(eff_err[c][i][j][t]);
				}
			}
		}
	}
}

/***************************************************************/
