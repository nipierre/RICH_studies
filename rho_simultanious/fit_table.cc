#include "fit_table.h"

// ******************************************************************************************
//
// Using the seelcted data for Lambda K0 and Phi in a simultanious fit to
// produce the RICH Table with pi,k,p
//
// ******************************************************************************************
const int start = 0;
const int stop = 8;

/**********************************************************************/
int main(int argc, char *argv[]){
	TH1F * tt = new TH1F("t","",10,0,10);
	TTree *ttt = new TTree();
	TLorentzVector* tttt = new TLorentzVector();
	delete tt;
	delete ttt;
	delete tttt;

	if(argc!=2){
		 cerr <<std::left<< "wrong options, possible options are: " << endl
			  << setw(10) << "fit"	 << "do the fit and produce the table" << endl
			  << setw(10) << "plots"  << "create the histograms for fitting"<< endl;
		return 1;
	}

	string mode = argv[1];

	if( mode == "test" ){ read_options(); return 0;}

	if(mode != "fit" && mode != "plots") {
		 cerr <<std::left << "wrong options, possible options are: " << endl
			  << setw(10) << "fit"	 << "do the fit and produce the table" << endl
			  << setw(10) << "plots"  << "create the histograms for fitting"<< endl;
		return 2;
	}

	// cerr << "check charge of particle and numbers (cc)" << endl;

	//RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
	gErrorIgnoreLevel = kWarning;

	if(mode == "plots"){
		read_options();
		create_hist();
		for(int i=0; i<data_nb; i++)
		{
			if(get_inputFile(i+data_ff_nb))
			{
				// t1 : K0 Lambdas
				if(analysis=="K0L") get_input_data_t1();
				// t2 : inc. Phis
				if(analysis=="iPhi") get_input_data_t2();
				// t3 : Rho
				if(analysis=="Rho") get_input_data_rho();
			}
		}
		write_hist();

		return 0;
	}
	if( mode == "fit"){
		read_options();
		get_plots();

		if(fit_type=="all")
		{
			fit_table_k0(0); // pi-
			fit_table_k0(1); // pi+
			fit_table_phi(0); // k-
			fit_table_phi(1); // k+
			fit_table_lambda(0); // p-
			fit_table_lambda(1); // p+
			fit_table_rho(0); // pi-
			fit_table_rho(1); // pi+
		}
		else if(fit_type=="pi")
		{
			fit_table_k0(0); // pi-
			fit_table_k0(1); // pi+
			fit_table_rho(0); // pi-
			fit_table_rho(1); // pi+
		}
		else if(fit_type=="k")
		{
			fit_table_phi(0); // k-
			fit_table_phi(1); // k+
		}
		else if(fit_type=="p")
		{
			fit_table_lambda(0); // p-
			fit_table_lambda(1); // p+
		}

	 	print_table();

		return 0;
	}

	return 0;
}

/**********************************************************************/
bool get_inputFile(int pi)
{
	delete input;
	// input = new TFile(Form("%s/%s_%d.root",data_file.c_str(),data_template.c_str(),pi),"READ");
	input = TFile::Open(Form("%s/%s_%d.root",data_file.c_str(),data_template.c_str(),pi));
	if(input)
		return 1;
	else
		return 0;
}

/**********************************************************************/
void read_options(){
	ifstream stream;
	string var1, var2;

	for(int i = 0; i<6; i++) for(int j = 0; j<4; j++)lh_cut[j][i] = -1;

	for(int i = 0; i<6; i++)
		for(int j = 0; j<5; j++)
			for(int t = 0; t< Nt; t++)
				for(int p = 0; p<Np; p++) N_id[i][j][p][t] = -1;

	stringstream nn;

	map< string ,pair<int,int> > opt;
// 	opt["LH_pi_pi:"] = make_pair(0,0);
	opt["LH_pi_K:"]  = make_pair(0,1);
	opt["LH_pi_p:"]  = make_pair(0,2);
	opt["LH_pi_e:"]  = make_pair(0,3);
	opt["LH_pi_mu:"] = make_pair(0,4);
	opt["LH_pi_bg:"] = make_pair(0,5);

	opt["LH_K_pi:"] = make_pair(1,0);
// 	opt["LH_K_K:"]  = make_pair(1,1);
	opt["LH_K_p:"]  = make_pair(1,2);
	opt["LH_K_e:"]  = make_pair(1,3);
	opt["LH_K_mu:"] = make_pair(1,4);
	opt["LH_K_bg:"] = make_pair(1,5);

	opt["LH_bthr_p_pi:"] = make_pair(2,0);
	opt["LH_bthr_p_K:"]  = make_pair(2,1);
	opt["LH_bthr_p_p:"]  = make_pair(2,2);
	opt["LH_bthr_p_e:"]  = make_pair(2,3);
	opt["LH_bthr_p_mu:"] = make_pair(2,4);
// 	opt["LH_bthr_p_bg:"] = make_pair(2,5);

	opt["LH_athr_p_pi:"] = make_pair(3,0);
	opt["LH_athr_p_K:"]  = make_pair(3,1);
// 	opt["LH_athr_p_p:"]  = make_pair(3,2);
	opt["LH_athr_p_e:"]  = make_pair(3,3);
	opt["LH_athr_p_mu:"] = make_pair(3,4);
	opt["LH_athr_p_bg:"] = make_pair(3,5);

	opt["LH_bthr_p_pi_bg:"] = make_pair(4,0);
	opt["LH_bthr_p_K_bg:"]  = make_pair(4,1);
	opt["LH_bthr_m_pi_bg:"] = make_pair(4,2);
	opt["LH_bthr_m_K_bg:"]  = make_pair(4,3);

	stream.open("options_fit.dat", ifstream::in);

	map<string, pair<int,int> >::iterator it;

	string line;

	if(!stream.good()) { cerr << "No option file" << endl; std::exit(0);}
	while(true){
		if (stream.eof()) break;
		getline(stream,line);
		if( line[0] != '#' && !line.empty() ) {

			nn.str("");
			nn.clear();
			nn.str(line);
			nn >> var1 >> var2;
			cout << var1 << " " << var2 << endl;
			if (var1 == "analysis:")		analysis = var2;
			if (var1 == "data_file:")		data_file = var2;
			if (var1 == "data_template:")		data_template = var2;
			if (var1 == "data_firstfile_nb:")		data_ff_nb = stoi(var2.c_str());
			if (var1 == "data_lastfile_nb:")		data_lf_nb = stoi(var2.c_str());
			if (var1 == "fit_type:")		fit_type = var2;
			if (var1 == "hist_file_k0:")	hist_file_k0 = var2;
			if (var1 == "hist_file_phi:")	hist_file_phi = var2;
			if (var1 == "hist_file_lam:")	hist_file_lam = var2;
			if (var1 == "hist_file_rho:")	hist_file_rho = var2;
			if (var1 == "out_file:")		out_file = var2;
			if (var1 == "line_width:")		stringstream ( var2 ) >> lw;
			if (var1 == "remove_richpipe:"){if(var2=="true") rpipe = true; else rpipe = false;}
			if (var1 == "max_retry:")		stringstream ( var2 ) >> retry;
			if (var1 == "thr_diff:")		stringstream ( var2 ) >> thr_diff;

			if (var1 == "rho_zlow:")		stringstream ( var2 ) >> rho_z1;
			if (var1 == "rho_zhigh:")		stringstream ( var2 ) >> rho_z2;

			if (var1 == "minuit_improve:")	{if(var2=="true") use_improve = true; else use_improve = false;}
			if (var1 == "minuit_hesse:")	{if(var2=="true") use_hesse = true; else use_hesse = false;}
			if (var1 == "minuit_minos:")	{if(var2=="true") use_minos = true; else use_minos = false;}
			if (var1 == "sidebins:")		{if(var2=="true") use_sidebins = true; else use_sidebins = false;}
			else{
				it = opt.find( var1 );
				if(it != opt.end()){
					stringstream ( var2 ) >> lh_cut[(*it).second.first][(*it).second.second];
				}
			}
		}
	}
	// cout << rho_z1 << " " << rho_z2 << endl;
	data_nb = data_lf_nb-data_ff_nb;
	cout << endl;
	// cout << rpipe << endl;
	stream.close();
}

//     0 - pion
//     1 - kaon
//     2 - proton
//     3 - electron (not available in old versions of RICH code)
//     4 - muon (not available in old versions of RICH code)
//     5 - background

/**********************************************************************/
void create_hist(){
	input = new TFile(Form("%s/%s_%d.root",data_file.c_str(),data_template.c_str(),data_ff_nb),"READ");
	stringstream nn;

	const string chan[8] = {"k0_pip","k0_pim","phi_kp","phi_km","lambda_pip","lambda_pim","rho_pip","rho_pim"};
	const string id[5]   = {"a","pi","k","p","u"};
	const Int_t Nbins[8]   = { 120,  120,    30,    30,   70,   70,  150,  150};
	const Double_t min[8]  = {0.44, 0.44, 0.995, 0.995,  1.1,  1.1, 550., 550.};
	const Double_t max[8]  = {0.56, 0.56, 1.042, 1.042, 1.13, 1.13, 1150.,1150.};

	test_hist=new TH2D("test_hist","",500,-1,1,500,0.,0.3);

	for(int i = 0; i<8; i++){      //8
		for(int j = 0; j<5; j++){  //4
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					nn.str("");
					nn.clear();
					nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
					h[i][j][p][t] = new TH1D(nn.str().c_str(),"",Nbins[i],min[i],max[i]);
					h[i][j][p][t]->Sumw2();

					nn.str("");
					nn.clear();
					nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
					h2[i][j][p][t] = new TH2D(nn.str().c_str(),"",100,-1.,1.,80,0.,0.4);

					nn.str("");
					nn.clear();
					nn << "all: " << p_bins[p]<< " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
					h2[i][j][p][t]->SetTitle(nn.str().c_str());
					h2[i][j][p][t]->GetXaxis()->SetTitle("#alpha");
					h2[i][j][p][t]->GetYaxis()->SetTitle("p_{t}");
				}
			}
		}
	}


	for(int i = 6; i<8; i++){      //8
		for(int j = 5; j<7; j++){
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					nn.str("");
					nn.clear();
					if(j==5)nn << "h_" << chan[i] << "_kk_" << p <<"_" <<t;
					if(j==6)nn << "h_" << chan[i] << "_la_" << p <<"_" <<t;
					if(j==5)h[i][j][p][t] = new TH1D(nn.str().c_str(),"",Nbins[i],min[i],max[i]);
					if(j==6)h[i][j][p][t] = new TH1D(nn.str().c_str(),"",20*Nbins[i],min[i],max[i]);
					h[i][j][p][t]->Sumw2();
				}
			}
		}
	}
}

/**********************************************************************/
void write_hist(){
	stringstream nn;

	const string chan[8] = {"k0_pip","k0_pim","phi_kp","phi_km","lambda_pip","lambda_pim","rho_pip","rho_pim"};
	const string id[5]   = {"a","pi","k","p","u"};


	for(int i =0; i<8; i++){
		cout << std::left
			<< setw(3) << i
			<< setw(15) << h[i][0][0][0]->GetEntries()
			<< setw(15) << h[i][1][0][0]->GetEntries()
			<< setw(15) << h[i][2][0][0]->GetEntries()
			<< setw(15) << h[i][3][0][0]->GetEntries()
			<< setw(15) << h[i][4][0][0]->GetEntries()
			<< setw(15) << h[i][0][0][0]->GetEntries() - h[i][1][0][0]->GetEntries() - h[i][2][0][0]->GetEntries() - h[i][3][0][0]->GetEntries() - h[i][4][0][0]->GetEntries()  << endl;
	}

	if(analysis=="K0L")
	{
		TFile* output = new TFile("hist_k0.root","RECREATE");
		for(int i = 0; i<2; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH1D*)elder->Get(nn.str().c_str()))->Write();
						((TH1D*)h[i][j][p][t])->Write();

						nn.str("");
						nn.clear();
						nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH2D*)elder->Get(nn.str().c_str()))->Write();
						((TH2D*)h2[i][j][p][t])->Write();
					}
				}
			}
		}
		output->Close();
		delete output;

		output = new TFile("hist_lambda.root","RECREATE");
		for(int i = 4; i<6; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH1D*)elder->Get(nn.str().c_str()))->Write();
						((TH1D*)h[i][j][p][t])->Write();

						nn.str("");
						nn.clear();
						nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH2D*)elder->Get(nn.str().c_str()))->Write();
						((TH2D*)h2[i][j][p][t])->Write();
					}
				}
			}
		}
		output->Close();
	}
	if(analysis=="iPhi")
	{
		TFile* output = new TFile("hist_incl_phi2.root","RECREATE");
		for(int i = 2; i<4; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// (TH1D*)elder->Get(nn.str().c_str());
						((TH1D*)h[i][j][p][t])->Write();
						//h->Write();

						nn.str("");
						nn.clear();
						nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH2D*)elder->Get(nn.str().c_str()))->Write();
						((TH2D*)h2[i][j][p][t])->Write();
					}
				}
			}
		}
		output->Close();
	}
	if(analysis=="Rho")
	{
		TFile* output = new TFile("hist_rho.root","RECREATE");
		for(int i = 6; i<8; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH1D*)input->Get(nn.str().c_str()))->Write();
						((TH1D*)h[i][j][p][t])->Write();

						nn.str("");
						nn.clear();
						nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						// ((TH2D*)input->Get(nn.str().c_str()))->Write();
						((TH2D*)h2[i][j][p][t])->Write();
					}
				}
			}
		}

		for(int i = 6; i<8; i++){
			int j = 5;     //8
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					nn.str("");
					nn.clear();
					nn << "h_" << chan[i] << "_kk_" << p <<"_" <<t;
					// ((TH1D*)input->Get(nn.str().c_str()))->Write();
					((TH1D*)h[i][j][p][t])->Write();

					nn.str("");
					nn.clear();
					nn << "h_" << chan[i] << "_la_" << p <<"_" <<t;
					// ((TH1D*)input->Get(nn.str().c_str()))->Write();
					((TH1D*)h[i][j][p][t])->Write();
				}
			}
		}
		output->Close();
	}

	TCanvas* c = new TCanvas("c","",800,600);
	test_hist->Draw("col");
	c->Print("test.C");
	c->Print("test.pdf");

	input->Close();

	delete input;
}

/**********************************************************************/
void get_input_data_t1(){
	TTree *tree = (TTree*)input->Get("tree");
	Long64_t nentries1 = tree->GetEntries();

	tree->SetBranchStatus("Run",0);
	tree->SetBranchStatus("Evt",0);
	tree->SetBranchStatus("TriggerMask",0);
	tree->SetBranchStatus("lv_beam",0);
	tree->SetBranchStatus("lv_scat",0);
	tree->SetBranchStatus("Q2",0);
	tree->SetBranchStatus("xbj",0);
	tree->SetBranchStatus("y",0);
	tree->SetBranchStatus("vx",0);
	tree->SetBranchStatus("vy",0);
	tree->SetBranchStatus("vz",0);
	tree->SetBranchStatus("v2x",0);
	tree->SetBranchStatus("v2y",0);
	tree->SetBranchStatus("v2z",0);
	tree->SetBranchStatus("Chi2",0);
	tree->SetBranchStatus("lv_p",0);
	tree->SetBranchStatus("lv_pi",0);
// 	tree->SetBranchStatus("pt1",0);
	tree->SetBranchStatus("pt2",0);
	tree->SetBranchStatus("pp_mom",0);
// 	tree->SetBranchStatus("pp_theta",0);
	tree->SetBranchStatus("pm_mom",0);
// 	tree->SetBranchStatus("pm_theta",0);

	stringstream cut;

	static TLorentzVector* lv_lambda;
	static TLorentzVector* lv_k0;
	static TLorentzVector* lv_pim;
	static TLorentzVector* lv_pip;

	static double alpha, pt1,pp_theta,pm_theta;
	static double pi_thr, k_thr, p_thr;

	static double pp_lh[7],pm_lh[7];
	for(int i = 0; i<7; i++) {
		pm_lh[i] = -1;
		pp_lh[i] = -1;
	}

	static double pp_x,  pp_y,  pm_x,  pm_y;

	tree->SetBranchAddress("lv_pip",&lv_pip);
	tree->SetBranchAddress("lv_pim",&lv_pim);
	tree->SetBranchAddress("lv_lambda",&lv_lambda);
	tree->SetBranchAddress("lv_k0",&lv_k0);
	tree->SetBranchAddress("alpha",&alpha);
	tree->SetBranchAddress("pp_x",&pp_x);
	tree->SetBranchAddress("pp_y",&pp_y);
	tree->SetBranchAddress("pm_x",&pm_x);
	tree->SetBranchAddress("pm_y",&pm_y);
	tree->SetBranchAddress("pm_lh",&pm_lh);
	tree->SetBranchAddress("pp_lh",&pp_lh);
	tree->SetBranchAddress("k_thr",&k_thr);
	tree->SetBranchAddress("pi_thr",&pi_thr);
	tree->SetBranchAddress("p_thr",&p_thr);
	tree->SetBranchAddress("pt1",&pt1);
	tree->SetBranchAddress("pp_theta",&pp_theta);
	tree->SetBranchAddress("pm_theta",&pm_theta);

	cout << nentries1  << endl;
	// cout << "tree 1" << endl;
	for (Long64_t jentry=0; jentry<nentries1;jentry++) {

		tree->GetEntry(jentry);
		int p_bin_m = -1;
		int p_bin_p = -1;
		int t_bin_m = -1;
		int t_bin_p = -1;

// 		continue;
		for(int i = 0 ; i<Np; i++){
			if(lv_pim->Vect().Mag() >= p_bins[i] && lv_pim->Vect().Mag() < p_bins[i+1]){
				if(lv_pip->Vect().Mag() > pi_thr){
					p_bin_m = i;
				}
			}
			if(lv_pip->Vect().Mag() >= p_bins[i] && lv_pip->Vect().Mag() < p_bins[i+1]){
				if(lv_pim->Vect().Mag() > pi_thr){
					p_bin_p = i;
				}
			}
		}
		for(int i = 0 ; i<Nt; i++){
// 			if(lv_pim->Vect().Theta() >= t_bins[i] && lv_pim->Vect().Theta() < t_bins[i+1]) t_bin_m = i;
// 			if(lv_pip->Vect().Theta() >= t_bins[i] && lv_pip->Vect().Theta() < t_bins[i+1]) t_bin_p = i;
			if(pm_theta >= t_bins[i] && pm_theta < t_bins[i+1]) t_bin_m = i;
			if(pp_theta >= t_bins[i] && pp_theta < t_bins[i+1]) t_bin_p = i;
		}

		if(p_bin_m == -1 && p_bin_p == -1 && t_bin_m == -1 && t_bin_p == -1) continue;

		int id_p=3;
		int id_m=3;


		const int id_lst[5] = {0,1,5,2,3};
		if(lv_pip->Vect().Mag()>pi_thr){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pp_lh[id_lst[i]] > lh_cut[i][0] * pp_lh[0]){
					if(lh_cut[i][1] == -1 || pp_lh[id_lst[i]] > lh_cut[i][1] * pp_lh[1]){
						if(lh_cut[i][2] == -1 || pp_lh[id_lst[i]] > lh_cut[i][2] * pp_lh[2]){
							if(lh_cut[i][3] == -1 || pp_lh[id_lst[i]] > lh_cut[i][3] * pp_lh[3]){
								if(lh_cut[i][4] == -1 || pp_lh[id_lst[i]] > lh_cut[i][4] * pp_lh[4]){
									if(lh_cut[i][5] == -1 || pp_lh[id_lst[i]] > lh_cut[i][5] * pp_lh[5]){
										if(i != 2 || (lv_pip->Vect().Mag()<=p_thr+thr_diff && id_p !=0 && id_p !=1)){
											if(i != 3 || lv_pip->Vect().Mag()>p_thr-thr_diff){
												id_p = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}


// 		if( lv_pip->Vect().Mag() <= p_thr+thr_diff && pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) id_p = 5;
		if( pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) id_p = 3;
		if( lv_pip->Vect().Mag() <= p_thr+thr_diff && ((pp_lh[0] == 0 && pp_lh[1] == 0 && pp_lh[2] == 0
																							 && pp_lh[3] == 0 && pp_lh[4] == 0 && pp_lh[5] == 0)
																						   || (pp_lh[0] < lh_cut[4][0] * pp_lh[5]
																							 && pp_lh[1] < lh_cut[4][1] * pp_lh[5])) ) id_p = 2;

		if(lv_pim->Vect().Mag()>pi_thr){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pm_lh[id_lst[i]] > lh_cut[i][0] * pm_lh[0]){
					if(lh_cut[i][1] == -1 || pm_lh[id_lst[i]] > lh_cut[i][1] * pm_lh[1]){
						if(lh_cut[i][2] == -1 || pm_lh[id_lst[i]] > lh_cut[i][2] * pm_lh[2]){
							if(lh_cut[i][3] == -1 || pm_lh[id_lst[i]] > lh_cut[i][3] * pm_lh[3]){
								if(lh_cut[i][4] == -1 || pm_lh[id_lst[i]] > lh_cut[i][4] * pm_lh[4]){
									if(lh_cut[i][5] == -1 || pm_lh[id_lst[i]] > lh_cut[i][5] * pm_lh[5]){
										if(i != 2 || (lv_pim->Vect().Mag()<=p_thr+thr_diff && id_m !=0 && id_m !=1)){
											if(i != 3 || lv_pim->Vect().Mag()>p_thr-thr_diff){
												id_m = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

// 		if( lv_pim->Vect().Mag() <= p_thr+thr_diff && pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1) id_m = 3;
		if( pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1) id_m = 3;
		if( lv_pim->Vect().Mag() <= p_thr+thr_diff && ((pm_lh[0] == 0 && pm_lh[1] == 0 && pm_lh[2] == 0
																							 && pm_lh[3] == 0 && pm_lh[4] == 0 && pm_lh[5] == 0)
																						   || (pm_lh[0] < lh_cut[4][2] * pm_lh[5]
																							 && pm_lh[1] < lh_cut[4][3] * pm_lh[5])) ) id_m = 2;

		if( (pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) || (pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1)) continue;

		int lambda = 0;
		if( alpha <0) lambda = -1;
		if( alpha >0) lambda =  1;

		bool pipe = false;
		if(rpipe){
			if(pp_x*pp_x + pp_y*pp_y >=25. && pm_x*pm_x + pm_y*pm_y >=25.) pipe = true;
		}

		if( pipe || !rpipe){
			if(p_bin_m != -1 && t_bin_m != -1  && id_p == 0){
				if(TMath::Abs(lv_lambda->Mag()-1.11568)>0.01){
					h[0][0][p_bin_m][t_bin_m]->Fill(lv_k0->Mag());
					h2[0][0][p_bin_m][t_bin_m]->Fill(alpha,pt1);
				}
				if(lambda == -1 && TMath::Abs(lv_k0->Mag()-0.497)>0.02 ){
					h[4][0][p_bin_m][t_bin_m]->Fill(lv_lambda->Mag());
					h2[4][0][p_bin_m][t_bin_m]->Fill(alpha,pt1);
				}
				for(int i = 0; i<5; i++){
					if(id_m == id_lst[i]){
						if(TMath::Abs(lv_lambda->Mag()-1.11568)>0.01){
							if(id_m!=5){
								h[0][id_m+1][p_bin_m][t_bin_m]->Fill(lv_k0->Mag());
								h2[0][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							} else{
								h[0][3][p_bin_m][t_bin_m]->Fill(lv_k0->Mag());
								h2[0][3][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							}
						}

						if(lambda == -1 && TMath::Abs(lv_k0->Mag()-0.497)>0.02 ){
							if(id_m!=5){
								h[4][id_m+1][p_bin_m][t_bin_m]->Fill(lv_lambda->Mag());
								h2[4][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							} else{
								h[4][3][p_bin_m][t_bin_m]->Fill(lv_lambda->Mag());
								h2[4][3][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							}
						}
					}
				}
			}
			if(p_bin_p != -1 && t_bin_p != -1 && id_m == 0){
				if(TMath::Abs(lv_lambda->Mag()-1.11568)>0.01){
					h[1][0][p_bin_p][t_bin_p]->Fill(lv_k0->Mag());
					h2[1][0][p_bin_p][t_bin_p]->Fill(alpha,pt1);
				}
				if(lambda == 1 && TMath::Abs(lv_k0->Mag()-0.497)>0.02 ){
					h[5][0][p_bin_p][t_bin_p]->Fill(lv_lambda->Mag());
					h2[5][0][p_bin_p][t_bin_p]->Fill(alpha,pt1);
				}
				for(int i = 0; i<5; i++){
					if(id_p == id_lst[i]){
						if(TMath::Abs(lv_lambda->Mag()-1.11568)>0.01){
							if(id_p!=5){
								h[1][id_p+1][p_bin_p][t_bin_p]->Fill(lv_k0->Mag());
								h2[1][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							} else{
								h[1][3][p_bin_p][t_bin_p]->Fill(lv_k0->Mag());
								h2[1][3][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							}
						}
						if(lambda == 1 && TMath::Abs(lv_k0->Mag()-0.497)>0.02 ){
							if(id_p!=5){
								h[5][id_p+1][p_bin_p][t_bin_p]->Fill(lv_lambda->Mag());
								h2[5][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							} else{
								h[5][3][p_bin_p][t_bin_p]->Fill(lv_lambda->Mag());
								h2[5][3][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							}
						}

					}
				}
			}
		}

	}

	delete tree;
}

/**********************************************************************/
void get_input_data_t2(){
	TTree *tree2 = (TTree*)input->Get("tree2");
	Long64_t nentries2 = tree2->GetEntries();

	cout << nentries2  << endl;
	// cout << "tree 2" << endl;

	double m_pi=0.139570;

	tree2->SetBranchStatus("Run",0);
// 	tree2->SetBranchStatus("Evt",0);
	tree2->SetBranchStatus("TriggerMask",0);
	tree2->SetBranchStatus("lv_beam",0);
	tree2->SetBranchStatus("lv_scat",0);
	tree2->SetBranchStatus("Q2",0);
	tree2->SetBranchStatus("alpha",0);
	tree2->SetBranchStatus("xbj",0);
	tree2->SetBranchStatus("y",0);
// 	tree2->SetBranchStatus("Emiss",0);
	tree2->SetBranchStatus("vx",0);
	tree2->SetBranchStatus("vy",0);
	tree2->SetBranchStatus("vz",0);
	tree2->SetBranchStatus("pt1",0);
	tree2->SetBranchStatus("pt2",0);
	tree2->SetBranchStatus("pp_mom",0);
// 	tree2->SetBranchStatus("pp_theta",0);
	tree2->SetBranchStatus("pm_mom",0);
// 	tree2->SetBranchStatus("pm_theta",0);

	stringstream cut;

	static Long64_t Evt;

	static TLorentzVector* lv_phi;
	static TLorentzVector* lv_kp;
	static TLorentzVector* lv_km;

	static TLorentzVector* lv_beam;
	static TLorentzVector* lv_scat;

	static double pi_thr2, k_thr2, p_thr2, alpha, pt1;
	static double pm_theta,pp_theta;

	static double pp_lh2[7],pm_lh2[7];
	for(int i = 0; i<7; i++) {
		pm_lh2[i] = -1;
		pp_lh2[i] = -1;
	}
	static int Nout;
	static double pp_x2, pp_y2, pm_x2, pm_y2,Emiss;

	tree2->SetBranchAddress("Evt",&Evt);
	tree2->SetBranchAddress("lv_phi",&lv_phi);
	tree2->SetBranchAddress("lv_kp",&lv_kp);
	tree2->SetBranchAddress("lv_km",&lv_km);
	tree2->SetBranchAddress("lv_beam",&lv_beam);
	tree2->SetBranchAddress("lv_scat",&lv_scat);
	tree2->SetBranchAddress("pp_x",&pp_x2);
	tree2->SetBranchAddress("pp_y",&pp_y2);
	tree2->SetBranchAddress("pm_x",&pm_x2);
	tree2->SetBranchAddress("pm_y",&pm_y2);
	tree2->SetBranchAddress("pm_lh",pm_lh2);
	tree2->SetBranchAddress("pp_lh",pp_lh2);
	tree2->SetBranchAddress("k_thr",&k_thr2);
	tree2->SetBranchAddress("pi_thr",&pi_thr2);
	tree2->SetBranchAddress("p_thr",&p_thr2);
	tree2->SetBranchAddress("pt1",&pt1);
	tree2->SetBranchAddress("alpha",&alpha);
	tree2->SetBranchAddress("pm_theta",&pm_theta);
	tree2->SetBranchAddress("pp_theta",&pp_theta);
	tree2->SetBranchAddress("Nout",&Nout);
	tree2->SetBranchAddress("Emiss",&Emiss);
	/*
	map<Long64_t,int> reject;

	for (Long64_t jentry=0; jentry<nentries2;jentry++) {
		tree2->GetEntry(jentry);

		if(lv_phi->Mag() >  1.042 || lv_phi->Mag() < 0.995) continue;

		map<Long64_t,int>::iterator it2;
		it2 = reject.find( Evt );

		if(it2 != reject.end()) (*it2).second++;
		else reject[Evt] = 1;
	}
	*/

	for (Long64_t jentry=0; jentry<nentries2;jentry++) {
		tree2->GetEntry(jentry);
		/*
		map<Long64_t,int>::iterator it2;
		it2 = reject.find( Evt );

		if(it2 != reject.end()){
			if((*it2).second>1){
				cout << (*it2).second<< endl;
				continue;
			}
		}
		*/
//		double real_Emiss = Emiss + lv_beam->E() - lv_scat->E() - lv_phi->E()+0.9382723/2;
		if( Nout == 3 && TMath::Abs(Emiss)<2.5 ) continue;
		TLorentzVector lv_pip;
		TLorentzVector lv_pim;
		TLorentzVector lv_ks1;
		TLorentzVector lv_ks2;

		lv_pip.SetVectM(lv_kp->Vect(),m_pi);
		lv_pim.SetVectM(lv_km->Vect(),m_pi);

		lv_ks1 = lv_pim + *lv_kp;
		lv_ks2 = lv_pip + *lv_km;
//		cout << lv_pip.Mag() << " " << lv_pim.Mag() << " " << lv_ks.Mag() << endl;
		if(TMath::Abs(lv_ks1.Mag()-0.89166)>0.05){
			if(TMath::Abs(lv_ks2.Mag()-0.89166)>0.05){
				test_hist->Fill(alpha,pt1);
			}
		}

		int p_bin_m = -1;
		int p_bin_p = -1;
		int t_bin_m = -1;
		int t_bin_p = -1;

		for(int i = 0 ; i<Np; i++){
			if(lv_km->Vect().Mag() >= p_bins[i] && lv_km->Vect().Mag() < p_bins[i+1]){
				if(lv_kp->Vect().Mag() > k_thr2){
					p_bin_m = i;
				}
			}
			if(lv_kp->Vect().Mag() >= p_bins[i] && lv_kp->Vect().Mag() < p_bins[i+1]){
				if(lv_km->Vect().Mag() > k_thr2){
					p_bin_p = i;
				}
			}
		}
		for(int i = 0 ; i<Nt; i++){
// 			if(lv_km->Vect().Theta() >= t_bins[i] && lv_km->Vect().Theta() < t_bins[i+1]) t_bin_m = i;
// 			if(lv_kp->Vect().Theta() >= t_bins[i] && lv_kp->Vect().Theta() < t_bins[i+1]) t_bin_p = i;
			if(pm_theta >= t_bins[i] && pm_theta < t_bins[i+1]) t_bin_m = i;
			if(pp_theta >= t_bins[i] && pp_theta < t_bins[i+1]) t_bin_p = i;
		}

		if(p_bin_m == -1 && p_bin_p == -1 && t_bin_m == -1 && t_bin_p == -1) continue;

		int id_p=3;
		int id_m=3;

		const int id_lst[5] = {0,1,5,2,3};

		if(lv_kp->Vect().Mag()>pi_thr2){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][0] * pp_lh2[0]){
					if(lh_cut[i][1] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][1] * pp_lh2[1]){
						if(lh_cut[i][2] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][2] * pp_lh2[2]){
							if(lh_cut[i][3] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][3] * pp_lh2[3]){
								if(lh_cut[i][4] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][4] * pp_lh2[4]){
									if(lh_cut[i][5] == -1 || pp_lh2[id_lst[i]] > lh_cut[i][5] * pp_lh2[5]){
										if(i != 2 || (lv_kp->Vect().Mag()<=p_thr2+thr_diff && id_p !=0 && id_p !=1)){
											if(i != 3 || lv_kp->Vect().Mag()>p_thr2-thr_diff){
												id_p = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
// 		if( lv_kp->Vect().Mag() <= p_thr2+thr_diff && pp_lh2[0] == -1 && pp_lh2[1] == -1 && pp_lh2[2] == -1 && pp_lh2[3] == -1 && pp_lh2[4] == -1 && pp_lh2[5] == -1) id_p = 5;
		if( pp_lh2[0] == -1 && pp_lh2[1] == -1 && pp_lh2[2] == -1 && pp_lh2[3] == -1 && pp_lh2[4] == -1 && pp_lh2[5] == -1) id_p = 3;
		if( lv_kp->Vect().Mag() <= p_thr2+thr_diff && ((pp_lh2[0] == 0 && pp_lh2[1] == 0 && pp_lh2[2] == 0
																							 && pp_lh2[3] == 0 && pp_lh2[4] == 0 && pp_lh2[5] == 0)
																						   || (pp_lh2[0] < lh_cut[4][0] * pp_lh2[5]
																							 && pp_lh2[1] < lh_cut[4][1] * pp_lh2[5])) ) id_p = 2;
// 		int des[4][8];
// 		for(int kk = 0; kk<4; kk++){
// 			for(int hh = 0 ; hh<8; hh++){
// 				des[kk][hh] = -9;
// 			}
// 		}
		if(lv_km->Vect().Mag()>pi_thr2){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][0] * pm_lh2[0]){
					if(lh_cut[i][1] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][1] * pm_lh2[1]){
						if(lh_cut[i][2] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][2] * pm_lh2[2]){
							if(lh_cut[i][3] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][3] * pm_lh2[3]){
								if(lh_cut[i][4] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][4] * pm_lh2[4]){
									if(lh_cut[i][5] == -1 || pm_lh2[id_lst[i]] > lh_cut[i][5] * pm_lh2[5]){
										if(i != 2 || (lv_km->Vect().Mag()<=p_thr2+thr_diff && id_m !=0 && id_m !=1)){
											if(i != 3 || lv_km->Vect().Mag()>p_thr2-thr_diff){
												id_m = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

// 		if( lv_km->Vect().Mag() <= p_thr2+thr_diff && pm_lh2[0] == -1 && pm_lh2[1] == -1 && pm_lh2[2] == -1 && pm_lh2[3] == -1 && pm_lh2[4] == -1 && pm_lh2[5] == -1) id_m = 5;
		if( pm_lh2[0] == -1 && pm_lh2[1] == -1 && pm_lh2[2] == -1 && pm_lh2[3] == -1 && pm_lh2[4] == -1 && pm_lh2[5] == -1) id_m = 3;
		if( lv_km->Vect().Mag() <= p_thr2+thr_diff && ((pm_lh2[0] == 0 && pm_lh2[1] == 0 && pm_lh2[2] == 0
																							 && pm_lh2[3] == 0 && pm_lh2[4] == 0 && pm_lh2[5] == 0)
																						   || (pm_lh2[0] < lh_cut[4][2] * pm_lh2[5]
																							 && pm_lh2[1] < lh_cut[4][3] * pm_lh2[5])) ) id_m = 2;

		if( (pp_lh2[0] == -1 && pp_lh2[1] == -1 && pp_lh2[2] == -1 && pp_lh2[3] == -1 && pp_lh2[4] == -1 && pp_lh2[5] == -1) || (pm_lh2[0] == -1 && pm_lh2[1] == -1 && pm_lh2[2] == -1 && pm_lh2[3] == -1 && pm_lh2[4] == -1 && pm_lh2[5] == -1)) continue;

		bool pipe = false;
		if(rpipe){
			if(pp_x2*pp_x2 + pp_y2*pp_y2 >=25. && pm_x2*pm_x2 + pm_y2*pm_y2 >=25.) pipe = true;
		}

		if( pipe || !rpipe){
			if(p_bin_m != -1 && t_bin_m != -1  && id_p == 1){
// 				if(TMath::Abs(lv_ks1.Mag()-0.89166)>0.05){
// 				if(lv_ks1.Mag()-0.89166<-0.05){
//					if(TMath::Abs(lv_ks2.Mag()-0.89166)>0.03){
						h[2][0][p_bin_m][t_bin_m]->Fill(lv_phi->Mag());
						h2[2][0][p_bin_m][t_bin_m]->Fill(alpha,pt1);
						for(int i = 0; i<5; i++){
							if(id_m == id_lst[i]){
								if(id_m!=5){
									h[2][id_m+1][p_bin_m][t_bin_m]->Fill(lv_phi->Mag());
									h2[2][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pt1);
								} else{
									h[2][3][p_bin_m][t_bin_m]->Fill(lv_phi->Mag());
									h2[2][3][p_bin_m][t_bin_m]->Fill(alpha,pt1);
								}
							}
						}
// 					}
// 				}
			}
			if(p_bin_p != -1 && t_bin_p != -1 && id_m == 1){
//				if(TMath::Abs(lv_ks1.Mag()-0.89166)>0.03){
// 					if(TMath::Abs(lv_ks2.Mag()-0.89166)>0.05){
// 					if(lv_ks2.Mag()-0.89166<-0.05){
						h[3][0][p_bin_p][t_bin_p]->Fill(lv_phi->Mag());
						h2[3][0][p_bin_p][t_bin_p]->Fill(alpha,pt1);
						for(int i = 0; i<5; i++){
							if(id_p == id_lst[i]){
								if(id_p!=5){
									h[3][id_p+1][p_bin_p][t_bin_p]->Fill(lv_phi->Mag());
									h2[3][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pt1);
								} else{
									h[3][3][p_bin_p][t_bin_p]->Fill(lv_phi->Mag());
									h2[3][3][p_bin_p][t_bin_p]->Fill(alpha,pt1);
								}
							}
						}
// 					}
//				}
			}
		}


	}

	delete tree2;

}

/**********************************************************************/
void get_input_data_rho(){
	TTree *tree = (TTree*)input->Get("tree");
	Long64_t nentries1 = tree->GetEntries();

	// tree->SetBranchStatus("nu",0);
	// tree->SetBranchStatus("Emiss",0);
	// tree->SetBranchStatus("Nout",0);
	tree->SetBranchStatus("pp_mom",0);
	tree->SetBranchStatus("pm_mom",0);

	stringstream cut;

	static TLorentzVector* lv_pim;
	static TLorentzVector* lv_pip;

	static double pp_theta,pm_theta;
	static double pi_thr, k_thr, p_thr, Emiss;
	static TLorentzVector* lv_beam;
	static TLorentzVector* lv_scat;
	static double m_rho = 0.776;
	static int Nout;

	static double pp_lh[7],pm_lh[7];
	for(int i = 0; i<7; i++) {
		pm_lh[i] = -1;
		pp_lh[i] = -1;
	}

	static double pp_x,  pp_y,  pm_x,  pm_y;
	tree->SetBranchAddress("lv_beam",&lv_beam);
	tree->SetBranchAddress("lv_scat",&lv_scat);
	tree->SetBranchAddress("k_thr",&k_thr);
	tree->SetBranchAddress("pi_thr",&pi_thr);
	tree->SetBranchAddress("p_thr",&p_thr);
	tree->SetBranchAddress("lv_pip",&lv_pip);
	tree->SetBranchAddress("lv_pim",&lv_pim);
	// tree->SetBranchAddress("m_rho",&m_rho);
	tree->SetBranchAddress("pp_x",&pp_x);
	tree->SetBranchAddress("pp_y",&pp_y);
	tree->SetBranchAddress("pp_theta",&pp_theta);
	tree->SetBranchAddress("pp_lh",&pp_lh);
	// tree->SetBranchAddress("pt1",&pt1);
	// tree->SetBranchAddress("alpha",&alpha);
	tree->SetBranchAddress("pm_x",&pm_x);
	tree->SetBranchAddress("pm_y",&pm_y);
	tree->SetBranchAddress("pm_theta",&pm_theta);
	tree->SetBranchAddress("pm_lh",&pm_lh);
	tree->SetBranchAddress("Nout",&Nout);
	tree->SetBranchAddress("Emiss",&Emiss);

	cout << "Entries: "<< nentries1  << endl;
	// cout << 2 << endl;

	// TH1D* hz[4];
	// hz[0] = new TH1D("hzep","",100,0,1);
	// hz[1] = new TH1D("hzem","",100,0,1);
	// hz[2] = new TH1D("hzip","",100,0,1);
	// hz[3] = new TH1D("hzim","",100,0,1);
	//
	// for(int i =0;i<4;i++){
		// hz[i]->GetXaxis()->SetTitle("z");
		// hz[i]->GetXaxis()->SetTitleOffset(1.01);
		// hz[i]->GetYaxis()->SetTitleOffset(1.01);
		// hz[i]->GetXaxis()->SetTitleSize(0.045);
		// hz[i]->GetXaxis()->SetLabelSize(0.045);
		// hz[i]->GetXaxis()->SetNdivisions(505);
		//
		// hz[i]->GetYaxis()->SetTitleSize(0.045);
		// hz[i]->GetYaxis()->SetTitleOffset(1.1);
		// hz[i]->GetYaxis()->SetLabelSize(0.045);
		// hz[i]->GetYaxis()->SetNdivisions(505);
	// }

	// TH1D* hde1 = new TH1D("hde1","",2000,-1,1);
	// TH1D* hde2 = new TH1D("hde2","",2000,-1,1);

	// TH1D* hde3 = new TH1D("hde3","",2000,-1,1);
	// TH1D* hde4 = new TH1D("hde4","",2000,-1,1);


	for (Long64_t jentry=0; jentry<nentries1;jentry++) {
		tree->GetEntry(jentry);
		int p_bin_m = -1;
		int p_bin_p = -1;
		int t_bin_m = -1;
		int t_bin_p = -1;

		// if(Nout != 3) continue;
		// if(TMath::Abs(Emiss) > 2.5) continue;

		if(m_rho < 0.55 || m_rho >= 0.95) continue;
		for(int i = 0 ; i<Np; i++){
			if(lv_pim->Vect().Mag() >= p_bins[i] && lv_pim->Vect().Mag() < p_bins[i+1]){
				if(lv_pip->Vect().Mag() > pi_thr){
					p_bin_m = i;
				}
			}
			if(lv_pip->Vect().Mag() >= p_bins[i] && lv_pip->Vect().Mag() < p_bins[i+1]){
				if(lv_pim->Vect().Mag() > pi_thr){
					p_bin_p = i;
				}
			}
		}
		for(int i = 0 ; i<Nt; i++){
// 			if(lv_pim->Vect().Theta() >= t_bins[i] && lv_pim->Vect().Theta() < t_bins[i+1]) t_bin_m = i;
// 			if(lv_pip->Vect().Theta() >= t_bins[i] && lv_pip->Vect().Theta() < t_bins[i+1]) t_bin_p = i;
			if(pm_theta >= t_bins[i] && pm_theta < t_bins[i+1]) t_bin_m = i;
			if(pp_theta >= t_bins[i] && pp_theta < t_bins[i+1]) t_bin_p = i;
		}

		if(p_bin_m == -1 && p_bin_p == -1 && t_bin_m == -1 && t_bin_p == -1) continue;

		int id_p=3;
		int id_m=3;

		const int id_lst[5] = {0,1,5,2,3};
		if(lv_pip->Vect().Mag()>pi_thr){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pp_lh[id_lst[i]] > lh_cut[i][0] * pp_lh[0]){
					if(lh_cut[i][1] == -1 || pp_lh[id_lst[i]] > lh_cut[i][1] * pp_lh[1]){
						if(lh_cut[i][2] == -1 || pp_lh[id_lst[i]] > lh_cut[i][2] * pp_lh[2]){
							if(lh_cut[i][3] == -1 || pp_lh[id_lst[i]] > lh_cut[i][3] * pp_lh[3]){
								if(lh_cut[i][4] == -1 || pp_lh[id_lst[i]] > lh_cut[i][4] * pp_lh[4]){
									if(lh_cut[i][5] == -1 || pp_lh[id_lst[i]] > lh_cut[i][5] * pp_lh[5]){
										if(i != 2 || (lv_pip->Vect().Mag()<=p_thr+thr_diff && id_p !=0 && id_p !=1)){
											if(i != 3 || lv_pip->Vect().Mag()>p_thr-thr_diff){
												id_p = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

// 		if( lv_pip->Vect().Mag() <= p_thr+thr_diff && pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) id_p = 5;
		if( pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) id_p = 3;
		if( lv_pip->Vect().Mag() <= p_thr+thr_diff && ((pp_lh[0] == 0 && pp_lh[1] == 0 && pp_lh[2] == 0
																							 && pp_lh[3] == 0 && pp_lh[4] == 0 && pp_lh[5] == 0)
																						   || (pp_lh[0] < lh_cut[4][0] * pp_lh[5]
																							 && pp_lh[1] < lh_cut[4][1] * pp_lh[5])) ) id_p = 2;

		if(lv_pim->Vect().Mag()>pi_thr){
			for(int i = 0; i<4; i++){
				if(lh_cut[i][0] == -1 || pm_lh[id_lst[i]] > lh_cut[i][0] * pm_lh[0]){
					if(lh_cut[i][1] == -1 || pm_lh[id_lst[i]] > lh_cut[i][1] * pm_lh[1]){
						if(lh_cut[i][2] == -1 || pm_lh[id_lst[i]] > lh_cut[i][2] * pm_lh[2]){
							if(lh_cut[i][3] == -1 || pm_lh[id_lst[i]] > lh_cut[i][3] * pm_lh[3]){
								if(lh_cut[i][4] == -1 || pm_lh[id_lst[i]] > lh_cut[i][4] * pm_lh[4]){
									if(lh_cut[i][5] == -1 || pm_lh[id_lst[i]] > lh_cut[i][5] * pm_lh[5]){
										if(i != 2 || (lv_pim->Vect().Mag()<=p_thr+thr_diff && id_m !=0 && id_m !=1)){
											if(i != 3 || lv_pim->Vect().Mag()>p_thr-thr_diff){
												id_m = id_lst[i];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

// 		if( lv_pim->Vect().Mag() <= p_thr+thr_diff && pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1) id_m = 3;
		if( pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1) id_m = 3;
		if( lv_pim->Vect().Mag() <= p_thr+thr_diff && ((pm_lh[0] == 0 && pm_lh[1] == 0 && pm_lh[2] == 0
																							 && pm_lh[3] == 0 && pm_lh[4] == 0 && pm_lh[5] == 0)
																						   || (pm_lh[0] < lh_cut[4][2] * pm_lh[5]
																							 && pm_lh[1] < lh_cut[4][3] * pm_lh[5])) ) id_m = 2;

		if( (pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) || (pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1)) continue;

		bool pipe = false;
		if(rpipe){
			if(pp_x*pp_x + pp_y*pp_y >=25. && pm_x*pm_x + pm_y*pm_y >=25.) pipe = true;
		}


		TLorentzVector lv_pm, lv_pp;
		lv_pp.SetVectM(lv_pip->Vect(),0.9382723);
		lv_pm.SetVectM(lv_pim->Vect(),0.9382723);

		TLorentzVector lv_km, lv_kp;
		lv_kp.SetVectM(lv_pip->Vect(),0.493677);
		lv_km.SetVectM(lv_pim->Vect(),0.493677);

		double z_pim = lv_pim->E()/(lv_beam->E()-lv_scat->E());
		double z_pip = lv_pip->E()/(lv_beam->E()-lv_scat->E());

		double m_kstar_kp = (lv_kp + *lv_pim).Mag();
		double m_kstar_km = (lv_km + *lv_pip).Mag();

		double m_lam_pp = (lv_pp + *lv_pim).Mag();
		double m_lam_pm = (lv_pm + *lv_pip).Mag();

		// double de_kp =  -TMath::Power((*lv_pip + *lv_pim).E(),2) + TMath::Power((lv_kp + *lv_pim).E(),2);
		// double de_km =  -TMath::Power((*lv_pip + *lv_pim).E(),2) + TMath::Power((lv_km + *lv_pip).E(),2);

		// de_kp *=1000;
		// de_kp *=1000;
		// de_km *=1000;
		// de_km *=1000;
		// if(TMath::Abs((lv_pp + *lv_pim).Mag()-1.11568)<0.02) continue;
		// if(TMath::Abs((lv_pm + *lv_pip).Mag()-1.11568)<0.02) continue;

		// if(TMath::Abs((lv_kp + *lv_pim).Mag()-0.89166)<0.05) continue;
		// if(TMath::Abs((lv_km + *lv_pip).Mag()-0.89166)<0.05) continue;

		// if((lv_kp + lv_km).Mag()<1.05) continue;


		if(m_rho <0.550 || m_rho >=0.950) continue;
		if( pipe || !rpipe){
			if(p_bin_m != -1 && t_bin_m != -1  && id_p == 0){
				if(z_pim >= rho_z1 && z_pim < rho_z2){
					h[6][0][p_bin_m][t_bin_m]->Fill(1000.*m_rho);
					// if(p_bin_m==11 && t_bin_m==1 && m_rho > 0.65 && m_rho<0.82){
						// hz[3]->Fill(lv_pim->E()/nu);
						// if(TMath::Abs(Emiss) < 2.5 && Nout ==3){
							// hz[1]->Fill(lv_pim->E()/nu);
						// }
					// }
					// h2[0][0][p_bin_m][t_bin_m]->Fill(alpha,pt1);
					for(int i = 0; i<5; i++){
						if(id_m == id_lst[i]){
							if(id_m==1){
								h[6][5][p_bin_m][t_bin_m]->Fill(1000.*m_kstar_km);
							}
							if(id_m==2 || id_m==5){
								h[6][6][p_bin_m][t_bin_m]->Fill(1000.*m_lam_pm);
							}
							if(id_m!=5){
								h[6][id_m+1][p_bin_m][t_bin_m]->Fill(1000.*m_rho);
								// h2[0][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							} else{
								h[6][3][p_bin_m][t_bin_m]->Fill(1000.*m_rho);
								// h2[0][3][p_bin_m][t_bin_m]->Fill(alpha,pt1);
							}
						}
					}
				}
			}
			if(p_bin_p != -1 && t_bin_p != -1 && id_m == 0){
				if(z_pip >= rho_z1 && z_pip < rho_z2){
					// if(de_kp>0.4){
					h[7][0][p_bin_p][t_bin_p]->Fill(1000.*m_rho);
					// if(p_bin_p==11 && t_bin_p==1 && m_rho > 0.65 && m_rho<0.82){
						// hz[2]->Fill(lv_pip->E()/nu);
						// if(TMath::Abs(Emiss) < 2.5 && Nout ==3){
							// hz[0]->Fill(lv_pip->E()/nu);
						// }
					// }
					// h2[1][0][p_bin_p][t_bin_p]->Fill(alpha,pt1);
					for(int i = 0; i<5; i++){
						if(id_p == id_lst[i]){
							if(id_p==1){
								h[7][5][p_bin_p][t_bin_p]->Fill(1000.*m_kstar_kp);
								// if(p_bin_p== 10 && t_bin_p == 1 && m_kstar_kp> 0.87&& m_kstar_kp <0.92) hde1->Fill(de_kp);
								// if(p_bin_p==  4 && t_bin_p == 1 && m_kstar_kp> 0.87&& m_kstar_kp <0.92) hde2->Fill(de_kp);
							}
							if(id_p==2 || id_p==5){
								h[7][6][p_bin_p][t_bin_p]->Fill(1000.*m_lam_pp);
							}
							// if(id_p==0){
								// if(p_bin_p== 10 && t_bin_p == 1 && m_rho> 0.7&& m_rho <0.8) hde3->Fill(de_kp);
								// if(p_bin_p==  4 && t_bin_p == 1 && m_rho> 0.7&& m_rho <0.8) hde4->Fill(de_kp);
							// }
							if(id_p!=5){
								h[7][id_p+1][p_bin_p][t_bin_p]->Fill(1000.*m_rho);
								// h2[1][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							} else{
								h[7][3][p_bin_p][t_bin_p]->Fill(1000.*m_rho);
								// h2[1][3][p_bin_p][t_bin_p]->Fill(alpha,pt1);
							}
						}
					}
					// }
				}
			}
		}
	}
	// TCanvas * cc = new TCanvas("cc","",800,600);
	// gStyle->SetOptStat(0);
	//
	// TLegend* leg= new TLegend(0.12,0.75,0.27,0.89);
	// leg->SetLineColor(0);
	// leg->SetNColumns(1);
	// leg->SetFillStyle(0);
	// leg->SetTextSize(0.04);
	// leg->AddEntry(hz[0], "#pi+",	"l");
	// leg->AddEntry(hz[1], "#pi-",	"l");
	//
	// hz[0]->SetLineColor(kRed+1);
	// hz[0]->SetLineWidth(2);
	// hz[0]->SetTitle("Exclusive #rho");
	//
	// hz[2]->SetLineColor(kRed+1);
	// hz[2]->SetLineWidth(2);
	// hz[2]->SetTitle("Inclusive #rho");
	//
	// hz[1]->SetLineColor(kAzure+2);
	// hz[1]->SetLineWidth(2);
	// hz[1]->SetTitle("Exclusive #rho");
	//
	// hz[3]->SetLineColor(kAzure+2);
	// hz[3]->SetLineWidth(2);
	// hz[3]->SetTitle("Inclusive #rho");
	//
	// hz[1]->Draw();
	// hz[0]->Draw("same");
	// leg->Draw();
	// cc->Print("z.pdf(");
	// cc->Clear();
	//
	// hz[2]->Draw();
	// hz[3]->Draw("same");
	// leg->Draw();
	// cc->Print("z.pdf)");





	// hde1->Draw();
	// cc->Print("de_1.root");
	// cc->Clear();

	// hde2->Draw();
	// cc->Print("de_2.root");
	// cc->Clear();

	// hde3->Draw();
	// cc->Print("de_3.root");
	// cc->Clear();

	// hde4->Draw();
	// cc->Print("de_4.root");
	// cc->Clear();

	delete tree;
}

/**********************************************************************/
void get_plots(){
	const string chan[8] = {"k0_pip","k0_pim","phi_kp","phi_km","lambda_pip","lambda_pim","rho_pip","rho_pim"};
	const string id[5]   = {"a","pi","k","p","u"};
	stringstream nn;
	// cout << 1 << endl;
	if(start<=0){
		// cout << "k0"<< endl;
		input_k0 = new TFile(hist_file_k0.c_str(),"READ");
		for(int i = 0; i<2; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						h[i][j][p][t] = (TH1D*)input_k0->Get(nn.str().c_str());
					}
				}
			}
		}
	}
	if(start<=2){
		// cout << "phi0"<< endl;
		input_phi = new TFile(hist_file_phi.c_str(),"READ");
		for(int i = 2; i<4; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						h[i][j][p][t] = (TH1D*)input_phi->Get(nn.str().c_str());
					}
				}
			}
		}
	}
	if(start<=4){
		// cout << "lambda"<< endl;
		input_lam = new TFile(hist_file_lam.c_str(),"READ");
		for(int i = 4; i<6; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						h[i][j][p][t] = (TH1D*)input_lam->Get(nn.str().c_str());
					}
				}
			}
		}
	}
	if(start<=6){
		// cout << "rho0"<< endl;
		input_rho = new TFile(hist_file_rho.c_str(),"READ");
		for(int i = 6; i<8; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						nn.str("");
						nn.clear();
						nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
						h[i][j][p][t] = (TH1D*)input_rho->Get(nn.str().c_str());
					}
				}
			}
			int j = 5;  //4
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					nn.str("");
					nn.clear();
					nn << "h_" << chan[i] << "_kk_" << p <<"_" <<t;
					h[i][j][p][t] = (TH1D*)input_rho->Get(nn.str().c_str());
				}
			}
			j = 6;
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					nn.str("");
					nn.clear();
					nn << "h_" << chan[i] << "_la_" << p <<"_" <<t;
					h[i][j][p][t] = (TH1D*)input_rho->Get(nn.str().c_str());
				}
			}
		}
	}
}

/**********************************************************************/
void rm_plots(){
	if(start<=0){
		for(int i = 0; i<2; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						delete h[i][j][p][t];
					}
				}
			}
		}
	}
	if(start<=2){
		for(int i = 2; i<4; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						delete h[i][j][p][t];
					}
				}
			}
		}
	}
	if(start<=4){
		for(int i = 4; i<6; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						delete h[i][j][p][t];
					}
				}
			}
		}
	}
	if(start<=6){
		for(int i = 6; i<8; i++){      //6
			for(int j = 0; j<5; j++){  //4
				for(int p = 0; p<Np;p++){
					for(int t = 0; t<Nt; t++){
						delete h[i][j][p][t];
					}
				}
			}
			int j = 5;  //4
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					delete h[i][j][p][t];
				}
			}
			j = 6;
			for(int p = 0; p<Np;p++){
				for(int t = 0; t<Nt; t++){
					delete h[i][j][p][t];
				}
			}
		}
	}
}

/**********************************************************************/
void fit_table_k0(int cc){
	// Model for K0
	if( cc != 0 && cc != 1) return;
	cout << endl;
	if(cc == 0) cout << "Fits of K0 sample for pi- efficiency:" << endl;
	if(cc == 1) cout << "Fits of K0 sample for pi+ efficiency:" << endl;
	for(int t = 0; t<Nt; t++){
		for(int p = 0; p<Np; p++){
// 	for(int p = 3; p<4; p++){
// 		for(int t = 1; t<2; t++){
			if(t==3 && p>6) continue;
			cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << endl;
			//Signal
			RooRealVar x("x","M",0.44,0.56,"GeV");
			x.setBins(120);

			RooRealVar mean("mean","mean",0.5,0.47,0.52,"GeV") ;
			RooRealVar sigma1("sigma1","sigma1",0.01,0.,0.1,"GeV");
			RooRealVar sigma2("sigma2","sigma2",0.005,0.,0.1,"GeV");

			RooGaussian gauss1("gauss1","gauss1",x,mean,sigma1) ;
			RooGaussian gauss2("gauss2","gauss2",x,mean,sigma2) ;

			//Background
			RooRealVar k0a("k0a","k0a",	-0.2,	-1.,	1.) ;
			RooRealVar k1a("k1a","k1a",	-0.3,	-1.,	1.) ;
			RooRealVar k2a("k2a","k2a",	0.,		-1.,	1.) ;
// 			RooRealVar k3a("k3a","k3a",	0.2,	-1.,	1.) ;
// 			RooRealVar k4a("k4a","k4a",	0.07,	-1.,	1.) ;
			RooChebychev bgna("bgna","bgna",x,RooArgSet(k0a,k1a,k2a));//,k3a,k4a)) ;


// 			RooRealVar k0pi("k0pi","k0pi",	-0.2,	-1.,	1.) ;
// 			RooRealVar k1pi("k1pi","k1pi",	-0.3,	-1.,	1.) ;
// 			RooRealVar k2pi("k2pi","k2pi",	0.,		-1.,	1.) ;
// 			RooRealVar k3pi("k3pi","k3pi",	0.2,	-1.,	1.) ;
// 			RooRealVar k4pi("k4pi","k4pi",	0.07,	-1.,	1.) ;
// 			RooChebychev bgnpi("bgnpi","bgnpi",x,RooArgSet(k0pi,k1pi));//,k2pi,k3pi,k4pi)) ;

// 			RooRealVar k0k("k0k","k0k",	-0.2,	-1.,	1.) ;
// 			RooRealVar k1k("k1k","k1k",	-0.3,	-1.,	1.) ;
// 			RooRealVar k2k("k2k","k2k",	0.,		-1.,	1.) ;
// 			RooRealVar k3k("k3k","k3k",	0.2,	-1.,	1.) ;
// 			RooRealVar k4k("k4k","k4k",	0.07,	-1.,	1.) ;
// 			RooChebychev bgnk("bgnk","bgnk",x,RooArgSet(k0k,k1k));//,k2k));//,k3k,k4k)) ;

			RooRealVar k0p("k0p","k0p",	-0.6,	-1.,	1.) ;
			RooRealVar k1p("k1p","k1p",	-0.4,	-1.,	1.) ;
			RooRealVar k2p("k2p","k2p",	0.08,	-1.,	1.) ;
// 			RooRealVar k3p("k3p","k3p",	0.2,	-1.,	1.) ;
// 			RooRealVar k4p("k4p","k4p",	0.07,	-1.,	1.) ;
			RooChebychev bgnp("bgnp","bgnp",x,RooArgSet(k0p,k1p,k2p));//,k3p,k4p)) ;

// 			RooRealVar k0u("k0u","k0u",	-0.2,	-1.,	1.) ;
// 			RooRealVar k1u("k1u","k1u",	-0.3,	-1.,	1.) ;
// 			RooRealVar k2u("k2u","k2u",	0.,		-1.,	1.) ;
// 			RooRealVar k3u("k3u","k3u",	0.2,	-1.,	1.) ;
// 			RooRealVar k4u("k4u","k4u",	0.07,	-1.,	1.) ;
// 			RooChebychev bgnu("bgnu","bgnu",x,RooArgSet(k0u,k1u,k2u));//,k3u,k4u)) ;

			//Construct composite pdf
			Int_t ent = h[0+cc][0][p][t]->GetEntries(); //[p][t]
// 			RooRealVar N_a_s("N_a_s","N_a_s",0.95*ent,0.,1.05*ent) ;
			RooRealVar N_a_b("N_a_b","N_a_b",0.05*ent,0.,1.05*ent) ;

			ent = h[0+cc][1][p][t]->GetEntries(); //[p][t]
			RooRealVar N_pi_s("N_pi_s","N_pi_s",0.9*ent,0.,1.05*ent) ;
			RooRealVar N_pi_b("N_pi_b","N_pi_b",0.1*ent,0.,1.05*ent) ;

			ent = h[0+cc][2][p][t]->GetEntries(); //[p][t]
			RooRealVar N_k_s("N_k_s","N_k_s",0.3*ent,0.,1.05*ent) ;
			RooRealVar N_k_b("N_k_b","N_k_b",0.7*ent,0.,1.05*ent) ;

			ent = h[0+cc][3][p][t]->GetEntries(); //[p][t]
			RooRealVar N_p_s("N_p_s","N_p_s",0.2*ent,0.,1.05*ent) ;
			RooRealVar N_p_b("N_p_b","N_p_b",0.8*ent,0.,1.05*ent) ;

			ent = h[0+cc][4][p][t]->GetEntries(); //[p][t]
			RooRealVar N_u_s("N_u_s","N_u_s",0.8*ent,0.,1.05*ent) ;
			RooRealVar N_u_b("N_u_b","N_u_b",0.05*ent,0.,1.05*ent) ;

			if(cc == 1 && t==2 && p >Np-4){
				cout << cc << " " << t << " " << p << endl;
				N_p_s.setVal(0.);
			}

			RooRealVar frac_s("frac_s","frac_s",0.65,0.5,0.7) ; //
			RooAddPdf sig("sig","sig",RooArgList(gauss1,gauss2),frac_s) ;

			RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));

			RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgna),RooArgList(N_a_s,N_a_b)) ;
			RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgna),RooArgList(N_pi_s,N_pi_b)) ; // checked with bgnpi
			RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgna),RooArgList(N_k_s,N_k_b)) ; // checked with bgnk
			RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgnp),RooArgList(N_p_s,N_p_b)) ;
			RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgna),RooArgList(N_u_s,N_u_b)) ;

			RooDataHist data_all("data_all","data_all",x,Import(*h[0+cc][0][p][t])); //[p][t]
			RooDataHist data_pi("data_pi","data_pi",x,Import(*h[0+cc][1][p][t])); //[p][t]
			RooDataHist data_k("data_k","data_k",x,Import(*h[0+cc][2][p][t])); //[p][t]
			RooDataHist data_p("data_p","data_p",x,Import(*h[0+cc][3][p][t])); //[p][t]
			RooDataHist data_unk("data_unk","data_unk",x,Import(*h[0+cc][4][p][t])); //[p][t]

			RooCategory sample("sample","sample") ;
			sample.defineType("all") ;
			sample.defineType("pi") ;
			sample.defineType("k") ;
			sample.defineType("p") ;
			sample.defineType("unk") ;

			RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));

			RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
			simPdf.addPdf(model_all,"all") ;
			simPdf.addPdf(model_pi,"pi") ;
			simPdf.addPdf(model_k,"k") ;
			simPdf.addPdf(model_p,"p") ;
			simPdf.addPdf(model_unk,"unk") ;

			if(!use_sidebins){
				RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
				RooAbsReal* nll = simPdf.createNLL(combData,Extended(true)) ;
				RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
				RooMinuit minu(nll_r) ;
				minu.setPrintLevel(-1);
				minu.setNoWarn();
				int i = 0;

				do{
// 					double delta = N_a_s.getVal() - (N_pi_s.getVal() + N_k_s.getVal() + N_p_s.getVal() + N_u_s.getVal());
					if(i>0){
// 						N_a_s.setVal( N_a_s.getVal()  - delta/2.);
						N_pi_s.setVal(1.01*N_pi_s.getVal());
						N_k_s.setVal( 0.98*N_k_s.getVal());
						N_p_s.setVal( 0.98*N_p_s.getVal());
						N_u_s.setVal( 0.98*N_u_s.getVal());
					}

					minu.hesse();
					minu.simplex();
					minu.migrad();
					minu.improve();
					minu.hesse();

					r[0+cc][p][t] = minu.save();

// 					delta = 1 - (N_pi_s.getVal() + N_k_s.getVal() + N_p_s.getVal() + N_u_s.getVal()) / N_a_s.getVal();
					i++;

				}while( r[0+cc][p][t]->covQual() != 3 && i<= retry);
				delete nll;
			}else {

				N_id[0+cc][0][p][t] = h[0+cc][0][p][t]->Integral(33,85) - h[0+cc][0][p][t]->Integral(3,29) - h[0+cc][0][p][t]->Integral(89,115);
				N_id[0+cc][1][p][t] = ( h[0+cc][1][p][t]->Integral(33,85) - h[0+cc][1][p][t]->Integral(3,29) - h[0+cc][1][p][t]->Integral(89,115) )/N_id[0+cc][0][p][t];
				N_id[0+cc][2][p][t] = ( h[0+cc][2][p][t]->Integral(33,85) - h[0+cc][2][p][t]->Integral(3,29) - h[0+cc][2][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];
				N_id[0+cc][3][p][t] = ( h[0+cc][3][p][t]->Integral(33,85) - h[0+cc][3][p][t]->Integral(3,29) - h[0+cc][3][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];
				N_id[0+cc][4][p][t] = ( h[0+cc][4][p][t]->Integral(33,85) - h[0+cc][4][p][t]->Integral(3,29) - h[0+cc][4][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];

			}

			N_id[0+cc][0][p][t] = N_a_s.getVal();
			N_id[0+cc][1][p][t] = N_pi_s.getVal();
			N_id[0+cc][2][p][t] = N_k_s.getVal() ;
			N_id[0+cc][3][p][t] = N_p_s.getVal() ;
			N_id[0+cc][4][p][t] = N_u_s.getVal() ;

			TGaxis::SetMaxDigits(3);
			stringstream nn;
			nn.str("");
			nn << "all: " << p_bins[p]<< " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
			RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
// 			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
			// simPdf.paramOn(frame1,Layout(0.8,0.99,0.85));
			// frame1->getAttText()->SetTextSize(0.02);

			RooPlot* frame2 = x.frame(Title("pi")) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;
// 			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;

			RooPlot* frame3 = x.frame(Title("k")) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;
// 			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;

			RooPlot* frame4 = x.frame(Title("p")) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;
// 			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgnp"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;

			RooPlot* frame5 = x.frame(Title("noID")) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;
// 			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;

			TCanvas* c = new TCanvas("c","",3200,600) ;

			frame1->GetXaxis()->SetTitleOffset(1.01);
			frame1->GetXaxis()->SetTitleSize(0.045);
			frame1->GetXaxis()->SetLabelSize(0.045);
			frame1->GetXaxis()->SetNdivisions(505);
			frame1->GetYaxis()->SetTitleOffset(1.3);
			frame1->GetYaxis()->SetTitleSize(0.045);
			frame1->GetYaxis()->SetLabelSize(0.045);
			frame1->GetYaxis()->SetNdivisions(505);

			frame2->GetXaxis()->SetTitleOffset(1.01);
			frame2->GetXaxis()->SetTitleSize(0.045);
			frame2->GetXaxis()->SetLabelSize(0.045);
			frame2->GetXaxis()->SetNdivisions(505);
			frame2->GetYaxis()->SetTitleOffset(1.3);
			frame2->GetYaxis()->SetTitleSize(0.045);
			frame2->GetYaxis()->SetLabelSize(0.045);
			frame2->GetYaxis()->SetNdivisions(505);

			frame3->GetXaxis()->SetTitleOffset(1.01);
			frame3->GetXaxis()->SetTitleSize(0.045);
			frame3->GetXaxis()->SetLabelSize(0.045);
			frame3->GetXaxis()->SetNdivisions(505);
			frame3->GetYaxis()->SetTitleOffset(1.3);
			frame3->GetYaxis()->SetTitleSize(0.045);
			frame3->GetYaxis()->SetLabelSize(0.045);
			frame3->GetYaxis()->SetNdivisions(505);

			frame4->GetXaxis()->SetTitleOffset(1.01);
			frame4->GetXaxis()->SetTitleSize(0.045);
			frame4->GetXaxis()->SetLabelSize(0.045);
			frame4->GetXaxis()->SetNdivisions(505);
			frame4->GetYaxis()->SetTitleOffset(1.3);
			frame4->GetYaxis()->SetTitleSize(0.045);
			frame4->GetYaxis()->SetLabelSize(0.045);
			frame4->GetYaxis()->SetNdivisions(505);

			frame5->GetXaxis()->SetTitleOffset(1.01);
			frame5->GetXaxis()->SetTitleSize(0.045);
			frame5->GetXaxis()->SetLabelSize(0.045);
			frame5->GetXaxis()->SetNdivisions(505);
			frame5->GetYaxis()->SetTitleOffset(1.3);
			frame5->GetYaxis()->SetTitleSize(0.045);
			frame5->GetYaxis()->SetLabelSize(0.045);
			frame5->GetYaxis()->SetNdivisions(505);

			c->Clear();
			TPad* pa[5];
			pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
			pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
			pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
			pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
			pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
			for(int dr = 0; dr<5;dr++){
				pa[dr]->SetLeftMargin(0.12);
				pa[dr]->SetRightMargin(0.08);
				pa[dr]->Draw();
			}
			pa[0]->cd() ; frame1->Draw() ;
			// TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
			// bb->SetY1(0.3);
			pa[1]->cd() ; frame2->Draw() ;
			pa[2]->cd() ; frame3->Draw() ;
			pa[3]->cd() ; frame4->Draw() ;
			pa[4]->cd() ; frame5->Draw() ;

			if(p==0 && t == 0){
				if(cc == 0) c->Print("test_k0_0.pdf[");
				else if(cc == 1) c->Print("test_k0_1.pdf[");
			}



			if(cc == 0) c->Print("test_k0_0.pdf");
			if(cc == 1) c->Print("test_k0_1.pdf");


			if(p==Np-1 && t == Nt-1){
			// if(p== 6 && t == Nt-1){
				if(cc == 0) c->Print("test_k0_0.pdf]");
				else if(cc == 1) c->Print("test_k0_1.pdf]");
			}
			delete frame1;
			delete frame2;
			delete frame3;
			delete frame4;
			delete frame5;
			for(int l=0;l<5;l++) delete pa[l];
			delete c;
		}
	}
}

/**********************************************************************/
void fit_table_phi(int cc){
	// Model for Phi
// 	ofstream res_out;
// 	if(cc == 0){
// 		res_out.open ("fit_out_0.dat");
// 	} else {
// 		res_out.open ("fit_out_1.dat");
// 	}
// 	res_out << std::left;
	if( cc != 0 && cc != 1) return;
	cout << endl;
	if(cc == 0) cout << "Fits of PHI sample for K- efficiency:" << endl;
	if(cc == 1) cout << "Fits of PHI sample for K+ efficiency:" << endl;

	for(int t = 0; t<Nt; t++){
// 	for(int t = 2; t<3; t++){
		for(int p = 0; p<Np; p++){

			if(t==3 && p>6) continue;
			cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << endl;
// 		for(int p = 0; p<2; p++){
 		//Signal
			RooRealVar x("x","M",0.995,1.042,"GeV");

			RooRealVar mean("mean","mean",1.01945,1.,1.04,"GeV") ;
			RooRealVar sigma("sigma","sigma",0.0031,0.0025,0.005,"GeV");
			RooRealVar width("width","width",0.00426);
			RooRealVar d1("d1","d1",1.);
			RooRealVar d2("d2","d2",0.4937);

			RooRelBreitWigner sb("sb","relBW",x,mean,width,d1,d2,d2,d1);
			RooGaussian sg("sg","gauss",x,mean,sigma);
//
//
			RooFFTConvPdf sig("sig","sig",x,sb,sg);
			sig.setShift(0,0);

			//Background
			RooRealVar thr("thr","threshold",0.987354,0.985,0.99,"GeV");
			RooRealVar thr2("thr2","threshold",0.987354,0.985,0.99,"GeV");
			RooRealVar thr3("thr3","threshold",0.987354,0.985,0.99,"GeV");
			RooRealVar thr4("thr4","threshold",0.987354,0.985,0.99,"GeV");
			RooRealVar thr5("thr5","threshold",0.987354,0.985,0.99,"GeV");

			RooRealVar n("n","n",  0.4,  0.,   1.5) ;
			RooRealVar a("a","a",  4.,   0.,    100.) ;
			RooGenericPdf bgn1("bgn1","x<thr ? 0 :TMath::Power(x - thr,n)*TMath::Exp(-a*(x - thr))",RooArgList(x,thr,n,a));

			RooRealVar n2("n2","n2",  0.4,  0.,   1.5) ;
			RooRealVar a2("a2","a2",  4.,   0.,    100.) ;
			RooGenericPdf bgn2("bgn2","x<thr2 ? 0 :TMath::Power(x - thr2,n2)*TMath::Exp(-a2*(x - thr2))",RooArgList(x,thr2,n2,a2));

			RooRealVar n3("n3","n3",  0.4,  0.,   1.5) ;
			RooRealVar a3("a3","a3",  4.,   0.,    100.) ;
			RooGenericPdf bgn3("bgn3","x<thr3 ? 0 :TMath::Power(x - thr3,n3)*TMath::Exp(-a3*(x - thr3))",RooArgList(x,thr3,n3,a3));

			RooRealVar n4("n4","n4",  0.4,  0.,   1.5) ;
			RooRealVar a4("a4","a4",  4.,   0.,    100.) ;
			RooGenericPdf bgn4("bgn4","x<thr4 ? 0 :TMath::Power(x - thr4,n4)*TMath::Exp(-a4*(x - thr4))",RooArgList(x,thr4,n4,a4));

			RooRealVar n5("n5","n5",  0.4,  0.,   1.5) ;
			RooRealVar a5("a5","a5",  4.,   0.,    100.) ;
			RooGenericPdf bgn5("bgn5","x<thr5 ? 0 :TMath::Power(x - thr5,n5)*TMath::Exp(-a5*(x - thr5))",RooArgList(x,thr5,n5,a5));


			//Construct composite pdf
			Int_t ent = h[2+cc][0][p][t]->GetEntries(); //[p][t]
// 			RooRealVar N_a_s("N_a_s","N_a_s1",0.8*ent,0.,1.05*ent) ;
			RooRealVar N_a_b("N_a_b","N_a_b",0.2*ent,0.,1.05*ent) ;

			ent = h[2+cc][1][p][t]->GetEntries(); //[p][t]
			RooRealVar N_pi_s("N_pi_s","N_pi_s",0.12*ent,0.,1.05*ent) ;
			RooRealVar N_pi_b("N_pi_b","N_pi_b",0.5*ent,0.,1.05*ent) ;
			RooRealVar N_pi_b2("N_pi_b2","N_pi_b2",0.38*ent,0.,1.05*ent) ;

			ent = h[2+cc][2][p][t]->GetEntries(); //[p][t]
			RooRealVar N_k_s("N_k_s","N_k_s",0.87*ent,0.,1.05*ent) ;
			RooRealVar N_k_b("N_k_b","N_k_b",0.13*ent,0.,1.05*ent) ;

			ent = h[2+cc][3][p][t]->GetEntries(); //[p][t]
			RooRealVar N_p_s("N_p_s","N_p_s",0.87*ent,0.,1.05*ent) ;
			RooRealVar N_p_b("N_p_b","N_p_b",0.13*ent,0.,1.05*ent) ;

			ent = h[2+cc][4][p][t]->GetEntries(); //[p][t]
			RooRealVar N_u_s("N_u_s","N_u_s",0.77*ent,0.,1.05*ent) ;
			RooRealVar N_u_b("N_u_b","N_u_b",0.23*ent,0.,1.05*ent) ;
			RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));
// 			RooFormulaVar N_a_b("N_a_b","N_pi_b + N_k_b + N_p_b + N_u_b",RooArgSet(N_pi_b,N_k_b,N_p_b,N_u_b));

			//bedingung nur für t = 2

			if(p>0){
				double offset1 = 1;
				double offset2 = 1;
// 				if(p==1||p==2){
// 					offset1=1.05;
// 					offset2=0.95;
// 				}
// 				if(t==1 && cc == 1 && (p==1||p==2)){
// 					offset1=1.05;
// 					offset2=0.95;
// 				}
				thr.setVal(last_para[0]);
				thr2.setVal(last_para[1]);
				thr3.setVal(last_para[2]);
				thr4.setVal(last_para[3]);
				thr5.setVal(last_para[4]);
				n.setVal(last_para[5]);
				n2.setVal(last_para[6]);
				n3.setVal(last_para[7]);
				n4.setVal(last_para[8]);
				n5.setVal(last_para[9]);
				a.setVal(last_para[10]);
				a2.setVal(last_para[11]);
				a3.setVal(last_para[12]);
				a4.setVal(last_para[13]);
				a5.setVal(last_para[14]);
				N_pi_s.setVal(offset2*last_para[15]);
				N_pi_b.setVal(offset2*last_para[16]);
				N_k_s.setVal(offset1*last_para[17]);
				N_k_b.setVal(offset1*last_para[18]);
				N_p_s.setVal(offset2*last_para[19]);
				N_p_b.setVal(offset2*last_para[20]);
				N_u_s.setVal(offset2*last_para[21]);
				N_u_b.setVal(offset2*last_para[22]);
				N_a_b.setVal(last_para[23]);
				mean.setVal(last_para[24]);
				sigma.setVal(last_para[25]);
				mean.setVal(last_para[26]);
			}



			RooArgList* lst_k_pdf = new RooArgList;
			lst_k_pdf->add(sig);

// 			if(t == 0 && p < 3) {
// 				lst_k_pdf->add(bgn1);
// 			} else {
				lst_k_pdf->add(bgn3);
// 			}

			RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgn1),RooArgList(N_a_s,N_a_b)) ;
			RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgn2),RooArgList(N_pi_s,N_pi_b)) ;
			RooAddPdf model_k("model_k","model_k",*lst_k_pdf,RooArgList(N_k_s,N_k_b)) ;
			RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgn1),RooArgList(N_p_s,N_p_b)) ;
			RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgn1),RooArgList(N_u_s,N_u_b)) ;

			x.setBins(30);
// 			x.setBins(70);
// 			x.setBins(140);
// 			x.setRange("fit_range",0.99,1.06);

// 			if(!use_sidebins) for(int i = 0 ; i< 5; i++) h[2+cc][i][p][t]->Rebin();

			RooDataHist data_all("data_all","data_all",x,Import(*h[2+cc][0][p][t])); //[p][t]
			RooDataHist data_pi("data_pi","data_pi",x,Import(*h[2+cc][1][p][t])); //[p][t]
			RooDataHist data_k("data_k","data_k",x,Import(*h[2+cc][2][p][t])); //[p][t]
			RooDataHist data_p("data_p","data_p",x,Import(*h[2+cc][3][p][t])); //[p][t]
			RooDataHist data_unk("data_unk","data_unk",x,Import(*h[2+cc][4][p][t])); //[p][t]

			RooCategory sample("sample","sample") ;
			sample.defineType("all") ;
			sample.defineType("pi") ;
			sample.defineType("k") ;
			sample.defineType("p") ;

			RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));
			RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
			simPdf.addPdf(model_all,"all") ;
			simPdf.addPdf(model_pi,"pi") ;
			simPdf.addPdf(model_k,"k") ;
			simPdf.addPdf(model_p,"p") ;
			simPdf.addPdf(model_unk,"unk") ;

			if(!use_sidebins){

// 				RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
				RooFormulaVar restriction("restriction","0",RooArgSet());

				RooAbsReal* nll = simPdf.createNLL(combData,Extended(true));
				RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
				RooMinuit minu(nll_r) ;
				minu.setPrintLevel(-1);
				minu.setNoWarn();
				int i = 0;

				do{
					if(i>0){
						if(t == 0){
							N_pi_s.setVal(0.95*N_pi_s.getVal());
							N_k_s.setVal( 1.05*N_k_s.getVal() );
							N_p_s.setVal( 0.90*N_p_s.getVal() );
							N_u_s.setVal( 0.95*N_u_s.getVal() );
						}else if(t == 1){
							N_pi_s.setVal(0.94*N_pi_s.getVal());
							N_k_s.setVal( 1.06*N_k_s.getVal() );
							N_p_s.setVal( 0.91*N_p_s.getVal() );
							N_u_s.setVal( 0.94*N_u_s.getVal() );
						}else if(t == 2){
							N_pi_s.setVal(0.91*N_pi_s.getVal());
							N_k_s.setVal( 1.09*N_k_s.getVal() );
							N_p_s.setVal( 0.91*N_p_s.getVal() );
							N_u_s.setVal( 0.92*N_u_s.getVal() );
						}else if(t == 3){
							N_pi_s.setVal(0.90*N_pi_s.getVal());
							N_k_s.setVal( 1.08*N_k_s.getVal() );
							N_p_s.setVal( 0.90*N_p_s.getVal() );
							N_u_s.setVal( 0.90*N_u_s.getVal() );
						}
					}

					minu.hesse();
// 					if( (p> Np-5 || p<4) && !(cc==1 && p ==2) && !(cc == 1 && p == Np-2) )
						minu.simplex();
					minu.migrad();
					minu.improve();
					minu.hesse();

					r[2+cc][p][t] = minu.save();

					i++;

// 				}while( TMath::Abs(1 - (N_pi_s.getVal() + N_k_s.getVal() + N_p_s.getVal() + N_u_s.getVal())/ N_a_s.getVal() ) >1e-4 && i<= retry);
				}while( r[2+cc][p][t]->covQual() != 3 && i<= retry);
// 					-1 "Unknown, matrix was externally provided"
// 					 0 "Not calculated at all"
// 					 1 "Approximation only, not accurate"
// 					 2 "Full matrix, but forced positive-definite"
// 					 3 "Full, accurate covariance matrix"
				delete nll;
			}else {

				N_id[2+cc][0][p][t] = h[2+cc][0][p][t]->Integral(26,55) - h[2+cc][0][p][t]->Integral(10,24) - h[2+cc][0][p][t]->Integral(57,72);
				N_id[2+cc][1][p][t] = ( h[2+cc][1][p][t]->Integral(26,55) - h[2+cc][1][p][t]->Integral(10,24) - h[2+cc][1][p][t]->Integral(57,72) )/N_id[2+cc][0][p][t];
				N_id[2+cc][2][p][t] = ( h[2+cc][2][p][t]->Integral(26,55) - h[2+cc][2][p][t]->Integral(10,24) - h[2+cc][2][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];
				N_id[2+cc][3][p][t] = ( h[2+cc][3][p][t]->Integral(26,55) - h[2+cc][3][p][t]->Integral(10,24) - h[2+cc][3][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];
				N_id[2+cc][4][p][t] = ( h[2+cc][4][p][t]->Integral(26,55) - h[2+cc][4][p][t]->Integral(10,24) - h[2+cc][4][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];

			}


// 			res_out << "t bin: " << t << " p bin: " << p << " Hesse0: " << hes0 <<" Migrad0: " << mig0<<" Migrad: " << mig << " Hesse: " << hes << " Minos: " << min  << "     " << h[2+cc][0][p][t]->GetEntries() -h[2+cc][1][p][t]->GetEntries() -h[2+cc][2][p][t]->GetEntries() -h[2+cc][3][p][t]->GetEntries() -h[2+cc][4][p][t]->GetEntries() <<"\n";
			N_id[2+cc][0][p][t] = N_a_s.getVal();
			N_id[2+cc][1][p][t] = N_pi_s.getVal();
			N_id[2+cc][2][p][t] = N_k_s.getVal();
			N_id[2+cc][3][p][t] = N_p_s.getVal();
			N_id[2+cc][4][p][t] = N_u_s.getVal();

			// if(p>0){
				last_para[0] = thr.getVal();
				last_para[1] = thr2.getVal();
				last_para[2] = thr3.getVal();
				last_para[3] = thr4.getVal();
				last_para[4] = thr5.getVal();
				last_para[5] = n.getVal();
				last_para[6] = n2.getVal();
				last_para[7] = n3.getVal();
				last_para[8] = n4.getVal();
				last_para[9] = n5.getVal();
				last_para[10] = a.getVal();
				last_para[11] = a2.getVal();
				last_para[12] = a3.getVal();
				last_para[13] = a4.getVal();
				last_para[14] = a5.getVal();
				last_para[15] = N_pi_s.getVal();
				last_para[16] = N_pi_b.getVal();
				last_para[17] = N_k_s.getVal();
				last_para[18] = N_k_b.getVal();
				last_para[19] = N_p_s.getVal();
				last_para[20] = N_p_b.getVal();
				last_para[21] = N_u_s.getVal();
				last_para[22] = N_u_b.getVal();
				last_para[23] = N_a_b.getVal();
				last_para[24] = mean.getVal();
				last_para[25] = sigma.getVal();
				last_para[26] = mean.getVal();
			// }
// 			N_id[2+cc][0][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_a_s"))->getVal();
// 			N_id[2+cc][1][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_pi"))->getVal();
// 			N_id[2+cc][2][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_k"))->getVal();
// 			N_id[2+cc][3][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_p"))->getVal();

// 			N_id[2+cc][0][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_a_s"))->getVal();
// 			N_id[2+cc][1][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_pi_s"))->getVal();
// 			N_id[2+cc][2][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_k_s"))->getVal();
// 			N_id[2+cc][3][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_p_s"))->getVal();

//			R_id[2+cc][0][p][t] = TMath::Sqrt(N_a_s1.getError()*N_a_s1.getError() + N_a_s2.getError()*N_a_s2.getError());
//			R_id[2+cc][1][p][t] = TMath::Sqrt(N_pi_s1.getError()*N_pi_s1.getError() + N_pi_s2.getError()*N_pi_s2.getError());
//			R_id[2+cc][2][p][t] = TMath::Sqrt(N_k_s1.getError()*N_k_s1.getError() + N_k_s2.getError()*N_k_s2.getError());
//			R_id[2+cc][3][p][t] = TMath::Sqrt(N_p_s1.getError()*N_p_s1.getError() + N_p_s2.getError()*N_p_s2.getError());
			stringstream nn;
			TGaxis::SetMaxDigits(3);
			nn.str("");
			nn << "all: " << p_bins[p]<< " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
			RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
// 			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
			// simPdf.paramOn(frame1,Layout(0.7,0.99,0.97));
			// frame1->getAttText()->SetTextSize(0.02);

			RooPlot* frame2 = x.frame(Title("pi")) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;
// 			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgn2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;

			RooPlot* frame3 = x.frame(Title("k")) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;
// 			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgn3"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;

			RooPlot* frame4 = x.frame(Title("p")) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;
// 			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;

			RooPlot* frame5 = x.frame(Title("noID")) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;
// 			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;


			TCanvas* c = new TCanvas("c","",3200,600) ;

			frame1->GetXaxis()->SetTitleOffset(1.01);
			frame1->GetXaxis()->SetTitleSize(0.045);
			frame1->GetXaxis()->SetLabelSize(0.045);
			frame1->GetXaxis()->SetNdivisions(505);
			frame1->GetYaxis()->SetTitleOffset(1.3);
			frame1->GetYaxis()->SetTitleSize(0.045);
			frame1->GetYaxis()->SetLabelSize(0.045);
			frame1->GetYaxis()->SetNdivisions(505);

			frame2->GetXaxis()->SetTitleOffset(1.01);
			frame2->GetXaxis()->SetTitleSize(0.045);
			frame2->GetXaxis()->SetLabelSize(0.045);
			frame2->GetXaxis()->SetNdivisions(505);
			frame2->GetYaxis()->SetTitleOffset(1.3);
			frame2->GetYaxis()->SetTitleSize(0.045);
			frame2->GetYaxis()->SetLabelSize(0.045);
			frame2->GetYaxis()->SetNdivisions(505);

			frame3->GetXaxis()->SetTitleOffset(1.01);
			frame3->GetXaxis()->SetTitleSize(0.045);
			frame3->GetXaxis()->SetLabelSize(0.045);
			frame3->GetXaxis()->SetNdivisions(505);
			frame3->GetYaxis()->SetTitleOffset(1.3);
			frame3->GetYaxis()->SetTitleSize(0.045);
			frame3->GetYaxis()->SetLabelSize(0.045);
			frame3->GetYaxis()->SetNdivisions(505);

			frame4->GetXaxis()->SetTitleOffset(1.01);
			frame4->GetXaxis()->SetTitleSize(0.045);
			frame4->GetXaxis()->SetLabelSize(0.045);
			frame4->GetXaxis()->SetNdivisions(505);
			frame4->GetYaxis()->SetTitleOffset(1.3);
			frame4->GetYaxis()->SetTitleSize(0.045);
			frame4->GetYaxis()->SetLabelSize(0.045);
			frame4->GetYaxis()->SetNdivisions(505);

			frame5->GetXaxis()->SetTitleOffset(1.01);
			frame5->GetXaxis()->SetTitleSize(0.045);
			frame5->GetXaxis()->SetLabelSize(0.045);
			frame5->GetXaxis()->SetNdivisions(505);
			frame5->GetYaxis()->SetTitleOffset(1.3);
			frame5->GetYaxis()->SetTitleSize(0.045);
			frame5->GetYaxis()->SetLabelSize(0.045);
			frame5->GetYaxis()->SetNdivisions(505);

			c->Clear();
			TPad* pa[5];
			pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
			pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
			pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
			pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
			pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
			for(int dr = 0; dr<5;dr++){
				pa[dr]->SetLeftMargin(0.12);
				pa[dr]->SetRightMargin(0.08);
				pa[dr]->Draw();
			}
			pa[0]->cd() ; frame1->Draw() ;
			// TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
			// bb->SetY1(0.3);
			pa[1]->cd() ; frame2->Draw() ;
			pa[2]->cd() ; frame3->Draw() ;
			pa[3]->cd() ; frame4->Draw() ;
			pa[4]->cd() ; frame5->Draw() ;

			if(p==0 && t == 0){
				if(cc == 0) c->Print("test_phi_0.pdf[");
				else if(cc == 1) c->Print("test_phi_1.pdf[");
			}

			if(cc == 0) c->Print("test_phi_0.pdf");
			if(cc == 1) c->Print("test_phi_1.pdf");


			if(p==Np-1 && t == Nt-1){
			// if(p== 6 && t == Nt-1){
				if(cc == 0) c->Print("test_phi_0.pdf]");
				else if(cc == 1) c->Print("test_phi_1.pdf]");
			}

// 			delete bb;
			delete lst_k_pdf;
			delete frame1;
			delete frame2;
			delete frame3;
			delete frame4;
			delete frame5;
			for(int l=0;l<5;l++) delete pa[l];
			delete c;

		}
	}
	first_phi = false;
// 	res_out.close();
}

/**********************************************************************/
void fit_table_lambda(int cc){
	// Model for Lambda
	retry = 1;
	if( cc != 0 && cc != 1) return;
	cout << endl;
	if(cc == 0) cout << "Fits of LAMBDA sample for pbar efficiency:" << endl;
	if(cc == 1) cout << "Fits of LAMBDA sample for p efficiency:" << endl;

	for(int t = 0; t<Nt;t++){
		for(int p = 0; p<Np;p++){
			for(int k = 0; k<30;k++){
				last_para2[p][t][k] = 0.;
			}
		}
	}


// 	for(int t = 0; t<Nt; t++){
	for(int t = 0; t<Nt; t++){
		for(int p = 0; p<Np; p++){
// 	for(int t = 1; t<2; t++){
// 		for(int p = 3; p<4; p++){
//Signal
			if(t==3 && p>6) continue;
			cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) <<p << endl;
			RooRealVar x("x","M",1.1,1.13,"GeV");
			x.setBins(70);
			RooRealVar mean("mean","mean",1.115,1.11,1.12,"GeV") ;
// 			RooRealVar mean("mean","mean",1.116) ;
// 			RooRealVar sigma2("sigma2","sigma2",0.002,0.0005,0.0021,"GeV");
			RooRealVar sigma1("sigma1","sigma1",0.004,0.0020,0.005,"GeV");
// 			RooRealVar sigma1("sigma1","sigma1",0.003);
			RooRealVar frac_s("frac_s","frac_s",0.15,0.1,0.7) ;

			RooFormulaVar sigma2("sigma2","sigma1*frac_s",RooArgSet(sigma1,frac_s));

			RooGaussian gauss1("gauss1","gauss1",x,mean,sigma1) ;
			RooGaussian gauss2("gauss2","gauss2",x,mean,sigma2) ;


			RooAddPdf sig("sig","sig",RooArgList(gauss1,gauss2),frac_s) ;
// 			RooVoigtian sig("sig","sig",x,mean,sigma2,sigma1);
			//Background
			RooRealVar n("n","n",  4.1,  0.,   7.5) ;
			RooRealVar a("a","a",  60.,   0.,    300.) ;
// 			RooRealVar thr("thr","threshold",1.0778423,1.,1.09,"GeV");
			RooRealVar thr("thr","threshold",1.0778423);
			RooGenericPdf bgna("bgna","x<thr ? 0 :TMath::Power(x - thr,n)*TMath::Exp(-a*(x - thr))",RooArgList(x,thr,n,a));

			RooRealVar n2("n2","n2",  2.7,  0.,   7.5) ;
			RooRealVar a2("a2","a2",  19.,   0.,    300.) ;
// 			RooRealVar thr2("thr2","threshold2",1.0778423,1.,1.09,"GeV");
			RooGenericPdf bgn2("bgn2","x<thr ? 0 :TMath::Power(x - thr,n2)*TMath::Exp(-a2*(x - thr))",RooArgList(x,thr,n2,a2));

			RooRealVar n3("n3","n3",  1.7,  0.,   7.5) ;
			RooRealVar a3("a3","a3",  31.,   0.,    300.) ;
// 			RooRealVar thr3("thr3","threshold3",1.0778423,1.,1.09,"GeV");
			RooGenericPdf bgn3("bgn3","x<thr ? 0 :TMath::Power(x - thr,n3)*TMath::Exp(-a3*(x - thr))",RooArgList(x,thr,n3,a3));

			RooRealVar n4("n4","n4",  1.7,  0.,   7.5) ;
			RooRealVar a4("a4","a4",  31.,   0.,    300.) ;
// 			RooRealVar thr4("thr4","threshold4",1.0778423,1.,1.09,"GeV");
			RooGenericPdf bgn4("bgn4","x<thr ? 0 :TMath::Power(x - thr,n4)*TMath::Exp(-a4*(x - thr))",RooArgList(x,thr,n4,a4));

			RooRealVar n5("n5","n5",  1.7,  0.,   7.5) ;
			RooRealVar a5("a5","a5",  31.,   0.,    300.) ;
// 			RooRealVar thr5("thr5","threshold5",1.0778423,1.,1.09,"GeV");
			RooGenericPdf bgn5("bgn5","x<thr ? 0 :TMath::Power(x - thr,n5)*TMath::Exp(-a5*(x - thr))",RooArgList(x,thr,n5,a5));

			//Construct composite pdf
			Int_t ent = h[4+cc][0][p][t]->GetEntries(); //[p][t]
// 			RooRealVar N_a_s("N_a_s","N_a_s",0.4*ent,0.,1.05*ent) ;
			RooRealVar N_a_b("N_a_b","N_a_b",0.45*ent,0.,1.05*ent) ;

			ent = h[4+cc][1][p][t]->GetEntries(); //[p][t]
			RooRealVar N_pi_s("N_pi_s","N_pi_s",0.13*ent,0.,1.05*ent) ;
			RooRealVar N_pi_b("N_pi_b","N_pi_b",0.82*ent,0.,1.05*ent) ;

			ent = h[4+cc][2][p][t]->GetEntries(); //[p][t]
			RooRealVar N_k_s("N_k_s","N_k_s",0.45*ent,0.,1.05*ent) ;
			RooRealVar N_k_b("N_k_b","N_k_b",0.5*ent,0.,1.05*ent) ;

			ent = h[4+cc][3][p][t]->GetEntries(); //[p][t]
			RooRealVar N_p_s("N_p_s","N_p_s",0.72*ent,0.,1.05*ent) ;
			RooRealVar N_p_b("N_p_b","N_p_b",0.11*ent,0.,1.05*ent) ;

			ent = h[4+cc][4][p][t]->GetEntries(); //[p][t]
			RooRealVar N_u_s("N_u_s","N_u_s",0.77*ent,0.,1.05*ent) ;
			RooRealVar N_u_b("N_u_b","N_u_b",0.11*ent,0.,1.05*ent) ;

			RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));
// 			RooFormulaVar N_a_b("N_a_b","N_pi_b + N_k_b + N_p_b + N_u_b",RooArgSet(N_pi_b,N_k_b,N_p_b,N_u_b));

			if(t == 1 && p == Np-2) {
				N_k_s.setVal(0);
				N_pi_s.setVal(0);
			}



			RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgna),RooArgList(N_a_s,N_a_b)) ;
			RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgna),RooArgList(N_pi_s,N_pi_b)) ;
			RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgna),RooArgList(N_k_s,N_k_b)) ;
			RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgna),RooArgList(N_p_s,N_p_b)) ;
			RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgna),RooArgList(N_u_s,N_u_b)) ;


			RooDataHist data_all("data_all","data_all",x,Import(*h[4+cc][0][p][t])); //[p][t]
			RooDataHist data_pi("data_pi","data_pi",x,Import(*h[4+cc][1][p][t])); //[p][t]
			RooDataHist data_k("data_k","data_k",x,Import(*h[4+cc][2][p][t])); //[p][t]
			RooDataHist data_p("data_p","data_p",x,Import(*h[4+cc][3][p][t])); //[p][t]
			RooDataHist data_unk("data_unk","data_unk",x,Import(*h[4+cc][4][p][t])); //[p][t]



			RooCategory sample("sample","sample") ;
			sample.defineType("all") ;
			sample.defineType("pi") ;
			sample.defineType("k") ;
			sample.defineType("p") ;
			sample.defineType("unk") ;

			RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));
			RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
			simPdf.addPdf(model_all,"all") ;
			simPdf.addPdf(model_pi,"pi") ;
			simPdf.addPdf(model_k,"k") ;
			simPdf.addPdf(model_p,"p") ;
			simPdf.addPdf(model_unk,"unk") ;


			if(p>0 && t==0){
				thr.setVal(last_para2[t][p-1][0]);
// 				thr2.setVal(last_para2[t][p-1][1]);
// 				thr3.setVal(last_para2[t][p-1][2]);
// 				thr4.setVal(last_para2[t][p-1][3]);
// 				thr5.setVal(last_para2[t][p-1][4]);
				n.setVal(last_para2[t][p-1][5]);
				n2.setVal(last_para2[t][p-1][6]);
				n3.setVal(last_para2[t][p-1][7]);
				n4.setVal(last_para2[t][p-1][8]);
				n5.setVal(last_para2[t][p-1][9]);
				a.setVal(last_para2[t][p-1][10]);
				a2.setVal(last_para2[t][p-1][11]);
				a3.setVal(last_para2[t][p-1][12]);
				a4.setVal(last_para2[t][p-1][13]);
				a5.setVal(last_para2[t][p-1][14]);
				N_pi_s.setVal(last_para2[t][p-1][15]);
				N_pi_b.setVal(last_para2[t][p-1][16]);
				N_k_s.setVal(last_para2[t][p-1][17]);
				N_k_b.setVal(last_para2[t][p-1][18]);
				N_p_s.setVal(last_para2[t][p-1][19]);
				N_p_b.setVal(last_para2[t][p-1][20]);
				N_u_s.setVal(last_para2[t][p-1][21]);
				N_u_b.setVal(last_para2[t][p-1][22]);
				N_a_b.setVal(last_para2[t][p-1][23]);
				mean.setVal(last_para2[t][p-1][24]);
				sigma1.setVal(last_para2[t][p-1][25]);
				frac_s.setVal(last_para2[t][p-1][26]);
			}else if(t>0){
				thr.setVal(last_para2[t-1][p][0]);
// 				thr2.setVal(last_para2[t][p-1][1]);
// 				thr3.setVal(last_para2[t][p-1][2]);
// 				thr4.setVal(last_para2[t][p-1][3]);
// 				thr5.setVal(last_para2[t][p-1][4]);
				n.setVal(last_para2[t-1][p][5]);
				n2.setVal(last_para2[t-1][p][6]);
				n3.setVal(last_para2[t-1][p][7]);
				n4.setVal(last_para2[t-1][p][8]);
				n5.setVal(last_para2[t-1][p][9]);
				a.setVal(last_para2[t-1][p][10]);
				a2.setVal(last_para2[t-1][p][11]);
				a3.setVal(last_para2[t-1][p][12]);
				a4.setVal(last_para2[t-1][p][13]);
				a5.setVal(last_para2[t-1][p][14]);
				N_pi_s.setVal(last_para2[t-1][p][15]);
				N_pi_b.setVal(last_para2[t-1][p][16]);
				N_k_s.setVal(last_para2[t-1][p][17]);
				N_k_b.setVal(last_para2[t-1][p][18]);
				N_p_s.setVal(last_para2[t-1][p][19]);
				N_p_b.setVal(last_para2[t-1][p][20]);
				N_u_s.setVal(last_para2[t-1][p][21]);
				N_u_b.setVal(last_para2[t-1][p][22]);
				N_a_b.setVal(last_para2[t-1][p][23]);
				mean.setVal(last_para2[t-1][p][24]);
				sigma1.setVal(last_para2[t-1][p][25]);
				frac_s.setVal(last_para2[t-1][p][26]);
			}

			if(!use_sidebins){
				RooFormulaVar restriction("restriction","0",RooArgSet());
// 				RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
				RooAbsReal* nll = simPdf.createNLL(combData,Extended(true)) ;
				RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
				RooMinuit minu(nll_r) ;

				minu.setPrintLevel(-1);
				minu.setNoWarn();
				int i = 0;

				do{
					mean.setVal(1.115) ;
					sigma1.setVal(0.004);
					frac_s.setVal(0.15) ;

					if(i>0){
						if(cc==0){
							N_pi_s.setVal(0.41*N_pi_s.getVal());
							N_k_s.setVal( 0.57*N_k_s.getVal() );
							N_p_s.setVal( 1.72*N_p_s.getVal() );
							N_u_s.setVal( 0.34*N_u_s.getVal() );

							N_pi_b.setVal( 1.02*N_pi_b.getVal() );
							N_k_b.setVal(  1.02*N_k_b.getVal() );
							N_p_b.setVal(  0.98*N_p_b.getVal() );
							N_u_b.setVal(  1.02*N_u_b.getVal() );
						}else{
							N_pi_s.setVal(0.33*N_pi_s.getVal());
							N_k_s.setVal( 0.42*N_k_s.getVal() );
							N_p_s.setVal( 1.84*N_p_s.getVal() );
							N_u_s.setVal( 0.27*N_u_s.getVal() );

							N_pi_b.setVal( 1.01*N_pi_b.getVal() );
							N_k_b.setVal(  1.01*N_k_b.getVal() );
							N_p_b.setVal(  0.99*N_p_b.getVal() );
							N_u_b.setVal(  1.01*N_u_b.getVal() );
						}

						// if(i%2 == 1){
							// N_pi_s.setVal((1.-i/((double)retry*1.07))*N_pi_s.getVal());
							// N_k_s.setVal( (1.-i/((double)retry*1.07))*N_k_s.getVal() );
							// N_p_s.setVal( (1.+i/((double)retry*1.07))*N_p_s.getVal() );
							// N_u_s.setVal( (1.-i/((double)retry*1.07))*N_u_s.getVal() );
						// }else if(i%2 == 0){
							// N_pi_s.setVal((1.+i/((double)retry*1.05))*N_pi_s.getVal());
							// N_k_s.setVal( (1.+i/((double)retry*1.05))*N_k_s.getVal() );
							// N_p_s.setVal( (1.-i/((double)retry*1.05))*N_p_s.getVal() );
							// N_u_s.setVal( (1.+i/((double)retry*1.05))*N_u_s.getVal() );
						// }


						/*
						if(t == 0){
							if(i%2 == 0){
								N_pi_s.setVal((1.-i/((double)retry*1.3))*N_pi_s.getVal());
								N_k_s.setVal( (1.-i/((double)retry*1.3))*N_k_s.getVal() );
								N_p_s.setVal( (1.+i/((double)retry*1.3))*N_p_s.getVal() );
								N_u_s.setVal( (1.-i/((double)retry*1.3))*N_u_s.getVal() );

// 								N_pi_s.setVal(0.90*N_pi_s.getVal());
// 								N_k_s.setVal( 0.90*N_k_s.getVal() );
// 								N_p_s.setVal( 1.10*N_p_s.getVal() );
// 								N_u_s.setVal( 0.90*N_u_s.getVal() );
							}else if(i%2 == 1){
// 								N_pi_s.setVal(1.10*N_pi_s.getVal());
// 								N_k_s.setVal( 1.10*N_k_s.getVal() );
// 								N_p_s.setVal( 0.90*N_p_s.getVal() );
// 								N_u_s.setVal( 1.10*N_u_s.getVal() );

								N_pi_s.setVal((1.+i/((double)retry*1.35))*N_pi_s.getVal());
								N_k_s.setVal( (1.+i/((double)retry*1.35))*N_k_s.getVal() );
								N_p_s.setVal( (1.-i/((double)retry*1.35))*N_p_s.getVal() );
								N_u_s.setVal( (1.+i/((double)retry*1.35))*N_u_s.getVal() );
							}
						}else if(t == 1){
							if(i%2 == 0){
								N_pi_s.setVal((1.-i/((double)retry*1.3))*N_pi_s.getVal());
								N_k_s.setVal( (1.-i/((double)retry*1.3))*N_k_s.getVal() );
								N_p_s.setVal( (1.+i/((double)retry*1.3))*N_p_s.getVal() );
								N_u_s.setVal( (1.-i/((double)retry*1.3))*N_u_s.getVal() );

// 								N_pi_s.setVal(0.90*N_pi_s.getVal());
// 								N_k_s.setVal( 0.90*N_k_s.getVal() );
// 								N_p_s.setVal( 1.10*N_p_s.getVal() );
// 								N_u_s.setVal( 0.90*N_u_s.getVal() );
							}else if(i%2 == 1){
// 								N_pi_s.setVal(1.10*N_pi_s.getVal());
// 								N_k_s.setVal( 1.10*N_k_s.getVal() );
// 								N_p_s.setVal( 0.90*N_p_s.getVal() );
// 								N_u_s.setVal( 1.10*N_u_s.getVal() );

								N_pi_s.setVal((1.+i/((double)retry*1.35))*N_pi_s.getVal());
								N_k_s.setVal( (1.+i/((double)retry*1.35))*N_k_s.getVal() );
								N_p_s.setVal( (1.-i/((double)retry*1.35))*N_p_s.getVal() );
								N_u_s.setVal( (1.+i/((double)retry*1.35))*N_u_s.getVal() );
							}
						}else if(t == 2){
							if(i%2 == 0){
								N_pi_s.setVal((1.-i/((double)retry*1.4))*N_pi_s.getVal());
								N_k_s.setVal( (1.-i/((double)retry*1.4))*N_k_s.getVal() );
								N_p_s.setVal( (1.+i/((double)retry*1.4))*N_p_s.getVal() );
								N_u_s.setVal( (1.-i/((double)retry*1.4))*N_u_s.getVal() );

// 								N_pi_s.setVal(0.50*N_pi_s.getVal());
// 								N_k_s.setVal( 0.50*N_k_s.getVal() );
// 								N_p_s.setVal( 1.50*N_p_s.getVal() );
// 								N_u_s.setVal( 0.50*N_u_s.getVal() );
							}else if(i%2 == 1){
								N_pi_s.setVal((1.+i/((double)retry*1.2))*N_pi_s.getVal());
								N_k_s.setVal( (1.+i/((double)retry*1.2))*N_k_s.getVal() );
								N_p_s.setVal( (1.-i/((double)retry*1.2))*N_p_s.getVal() );
								N_u_s.setVal( (1.+i/((double)retry*1.2))*N_u_s.getVal() );

// 								N_pi_s.setVal(1.50*N_pi_s.getVal());
// 								N_k_s.setVal( 1.50*N_k_s.getVal() );
// 								N_p_s.setVal( 0.50*N_p_s.getVal() );
// 								N_u_s.setVal( 1.50*N_u_s.getVal() );
							}
						}else if(t == 3){
							if(i%2 == 0){
								N_pi_s.setVal((1.-i/((double)retry*1.4))*N_pi_s.getVal());
								N_k_s.setVal( (1.-i/((double)retry*1.4))*N_k_s.getVal() );
								N_p_s.setVal( (1.+i/((double)retry*1.4))*N_p_s.getVal() );
								N_u_s.setVal( (1.-i/((double)retry*1.4))*N_u_s.getVal() );

// 								N_pi_s.setVal(0.50*N_pi_s.getVal());
// 								N_k_s.setVal( 0.50*N_k_s.getVal() );
// 								N_p_s.setVal( 1.50*N_p_s.getVal() );
// 								N_u_s.setVal( 0.50*N_u_s.getVal() );
							}else if(i%2 == 1){
								N_pi_s.setVal((1.+i/((double)retry*1.3))*N_pi_s.getVal());
								N_k_s.setVal( (1.+i/((double)retry*1.3))*N_k_s.getVal() );
								N_p_s.setVal( (1.-i/((double)retry*1.3))*N_p_s.getVal() );
								N_u_s.setVal( (1.+i/((double)retry*1.3))*N_u_s.getVal() );

// 								N_pi_s.setVal(1.50*N_pi_s.getVal());
// 								N_k_s.setVal( 1.50*N_k_s.getVal() );
// 								N_p_s.setVal( 0.50*N_p_s.getVal() );
// 								N_u_s.setVal( 1.50*N_u_s.getVal() );
							}
						}*/
					}
					minu.hesse();
					minu.simplex();
					minu.migrad();
					minu.improve();
					minu.hesse();
					r[4+cc][p][t] = minu.save();

					i++;

				}while( r[4+cc][p][t]->covQual() != 3 && i<= retry);
				delete nll;
			}else {

				N_id[4+cc][0][p][t] = h[4+cc][0][p][t]->Integral(18,37) - h[4+cc][0][p][t]->Integral(6,15) - h[4+cc][0][p][t]->Integral(40,50);
				N_id[4+cc][1][p][t] = ( h[4+cc][1][p][t]->Integral(18,37) - h[4+cc][1][p][t]->Integral(6,15) - h[4+cc][1][p][t]->Integral(40,50) )/N_id[4+cc][0][p][t];
				N_id[4+cc][2][p][t] = ( h[4+cc][2][p][t]->Integral(18,37) - h[4+cc][2][p][t]->Integral(6,15) - h[4+cc][2][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];
				N_id[4+cc][3][p][t] = ( h[4+cc][3][p][t]->Integral(18,37) - h[4+cc][3][p][t]->Integral(6,15) - h[4+cc][3][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];
				N_id[4+cc][4][p][t] = ( h[4+cc][4][p][t]->Integral(18,37) - h[4+cc][4][p][t]->Integral(6,15) - h[4+cc][4][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];

			}



			N_id[4+cc][0][p][t] = N_a_s.getVal();
			N_id[4+cc][1][p][t] = N_pi_s.getVal();
			N_id[4+cc][2][p][t] = N_k_s.getVal() ;
			N_id[4+cc][3][p][t] = N_p_s.getVal() ;
			N_id[4+cc][4][p][t] = N_u_s.getVal() ;


			if(p>-1){
				last_para2[p][t][0] = thr.getVal();
// 				last_para2[p][t][1] = thr2.getVal();
// 				last_para2[p][t][2] = thr3.getVal();
// 				last_para2[p][t][3] = thr4.getVal();
// 				last_para2[p][t][4] = thr5.getVal();
				last_para2[p][t][5] = n.getVal();
				last_para2[p][t][6] = n2.getVal();
				last_para2[p][t][7] = n3.getVal();
				last_para2[p][t][8] = n4.getVal();
				last_para2[p][t][9] = n5.getVal();
				last_para2[p][t][10] = a.getVal();
				last_para2[p][t][11] = a2.getVal();
				last_para2[p][t][12] = a3.getVal();
				last_para2[p][t][13] = a4.getVal();
				last_para2[p][t][14] = a5.getVal();
				last_para2[p][t][15] = N_pi_s.getVal();
				last_para2[p][t][16] = N_pi_b.getVal();
				last_para2[p][t][17] = N_k_s.getVal();
				last_para2[p][t][18] = N_k_b.getVal();
				last_para2[p][t][19] = N_p_s.getVal();
				last_para2[p][t][20] = N_p_b.getVal();
				last_para2[p][t][21] = N_u_s.getVal();
				last_para2[p][t][22] = N_u_b.getVal();
				last_para2[p][t][23] = N_a_b.getVal();
				last_para2[p][t][24] = mean.getVal();
				last_para2[p][t][25] = sigma1.getVal();
				last_para2[p][t][26] = frac_s.getVal();
			}


			stringstream nn;
			TGaxis::SetMaxDigits(3);
			nn.str("");
			nn << "all: " << p_bins[p]<< " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
			RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
// 			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
// 			simPdf.plotOn(frame1,Slice(sample,"all"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
			// simPdf.paramOn(frame1,Layout(0.7,0.99,0.97));
			// frame1->getAttText()->SetTextSize(0.02);

			RooPlot* frame2 = x.frame(Title("pi")) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;
// 			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
// 			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;

			RooPlot* frame3 = x.frame(Title("k")) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;
// 			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
// 			simPdf.plotOn(frame3,Slice(sample,"k"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;

			RooPlot* frame4 = x.frame(Title("p")) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;
// 			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
// 			simPdf.plotOn(frame4,Slice(sample,"p"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;


			RooPlot* frame5 = x.frame(Title("noID")) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;
// 			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
// 			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;

			TCanvas* c = new TCanvas("c","",3200,600) ;
			frame1->GetXaxis()->SetTitleOffset(1.01);
			frame1->GetXaxis()->SetTitleSize(0.045);
			frame1->GetXaxis()->SetLabelSize(0.045);
			frame1->GetXaxis()->SetNdivisions(505);
			frame1->GetYaxis()->SetTitleOffset(1.3);
			frame1->GetYaxis()->SetTitleSize(0.045);
			frame1->GetYaxis()->SetLabelSize(0.045);
			frame1->GetYaxis()->SetNdivisions(505);

			frame2->GetXaxis()->SetTitleOffset(1.01);
			frame2->GetXaxis()->SetTitleSize(0.045);
			frame2->GetXaxis()->SetLabelSize(0.045);
			frame2->GetXaxis()->SetNdivisions(505);
			frame2->GetYaxis()->SetTitleOffset(1.3);
			frame2->GetYaxis()->SetTitleSize(0.045);
			frame2->GetYaxis()->SetLabelSize(0.045);
			frame2->GetYaxis()->SetNdivisions(505);

			frame3->GetXaxis()->SetTitleOffset(1.01);
			frame3->GetXaxis()->SetTitleSize(0.045);
			frame3->GetXaxis()->SetLabelSize(0.045);
			frame3->GetXaxis()->SetNdivisions(505);
			frame3->GetYaxis()->SetTitleOffset(1.3);
			frame3->GetYaxis()->SetTitleSize(0.045);
			frame3->GetYaxis()->SetLabelSize(0.045);
			frame3->GetYaxis()->SetNdivisions(505);

			frame4->GetXaxis()->SetTitleOffset(1.01);
			frame4->GetXaxis()->SetTitleSize(0.045);
			frame4->GetXaxis()->SetLabelSize(0.045);
			frame4->GetXaxis()->SetNdivisions(505);
			frame4->GetYaxis()->SetTitleOffset(1.3);
			frame4->GetYaxis()->SetTitleSize(0.045);
			frame4->GetYaxis()->SetLabelSize(0.045);
			frame4->GetYaxis()->SetNdivisions(505);

			frame5->GetXaxis()->SetTitleOffset(1.01);
			frame5->GetXaxis()->SetTitleSize(0.045);
			frame5->GetXaxis()->SetLabelSize(0.045);
			frame5->GetXaxis()->SetNdivisions(505);
			frame5->GetYaxis()->SetTitleOffset(1.3);
			frame5->GetYaxis()->SetTitleSize(0.045);
			frame5->GetYaxis()->SetLabelSize(0.045);
			frame5->GetYaxis()->SetNdivisions(505);

			c->Clear();
			TPad* pa[5];
			pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
			pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
			pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
			pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
			pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
			for(int dr = 0; dr<5;dr++){
				pa[dr]->SetLeftMargin(0.12);
				pa[dr]->SetRightMargin(0.08);
				pa[dr]->Draw();
			}
			pa[0]->cd() ; frame1->Draw() ;
			// TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
			// bb->SetY1(0.3);
			pa[1]->cd() ; frame2->Draw() ;
			pa[2]->cd() ; frame3->Draw() ;
			pa[3]->cd() ; frame4->Draw() ;
			pa[4]->cd() ; frame5->Draw() ;

			if(p==0 && t == 0){
				if(cc == 0) c->Print("test_lambda_0.pdf[");
				else if(cc == 1) c->Print("test_lambda_1.pdf[");
			}

			if(cc == 0) c->Print("test_lambda_0.pdf");
			if(cc == 1) c->Print("test_lambda_1.pdf");

			if(p==Np-1 && t == Nt-1){
			// if(p== 6 && t == Nt-1){
				if(cc == 0) c->Print("test_lambda_0.pdf]");
				else if(cc == 1) c->Print("test_lambda_1.pdf]");
			}

// 			delete bb;
			delete frame1;
			delete frame2;
			delete frame3;
			delete frame4;
			delete frame5;
			for(int l=0;l<5;l++) delete pa[l];
			delete c;
		}
	}
}

/**********************************************************************/
void fit_table_rho(int cc){
	// Model for Rho0
	if( cc != 0 && cc != 1) return;

	if(cc == 0) cout << "Fits of Rho0 sample for pi- efficiency:" << endl;
	if(cc == 1) cout << "Fits of Rho0 sample for pi+ efficiency:" << endl;
	for(int t = 0; t<Nt; t++){
		for(int p = 0; p<Np; p++){

			// if(p!= 10) continue;
			// if(t!=1 ) continue;
			// 
			// if(t==3 && p>6) continue;

			cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << endl;
			//Signal

			RooRealVar x("x","M",550.,1150,"MeV");
			x.setBins(150);
			x.setRange("fitRange_all",550.,950.);
			x.setRange("fitRange_pi",550.,950.);
			x.setRange("fitRange_k",550.,950.);
			x.setRange("fitRange_p",550.,950.);
			x.setRange("fitRange_unk",550.,950.);
			x.setRange("fitRange_ks",780.,1050.);
			x.setRange("fitRange_la",950.,1150.);

			// RooRealVar x("x","M",550.,950,"MeV");
			// x.setBins(100);

			/************************/
			/*   test to fit K*     */
			/************************/

			RooRealVar k_mean("k_mean","mean",892,890,894,"MeV") ;
			RooRealVar k_sigma1("k_sigma1","sigma1",34,25,50,"MeV");
			RooRealVar k_sigma2("k_sigma2","sigma2",16,10,20,"MeV");
			RooRealVar k_frac_s("k_frac_s","frac_s",0.67,0.6,0.75) ;

			RooGaussian k_gauss1("k_gauss1","gauss1",x,k_mean,k_sigma1) ;
			RooGaussian k_gauss2("k_gauss2","gauss2",x,k_mean,k_sigma2) ;

			RooAddPdf k_sig("k_sig","sig",RooArgList(k_gauss1,k_gauss2),k_frac_s) ;

			RooRealVar k_n("k_n","k_n",  4.1,  0.,   7.5) ;
			RooRealVar k_a("k_a","k_a",  0.015,   0.,    300.) ;
			RooRealVar k_thr("k_thr","threshold",633.247);
			RooGenericPdf k_bgn("k_bgn","x<k_thr ? 0 :TMath::Power(x - k_thr,k_n)*TMath::Exp(-k_a*(x - k_thr))",RooArgList(x,k_thr,k_n,k_a));

			// RooRealVar k_k0("k_k0","k_k0",	 0.2,	-1.,	1.) ;
			// RooRealVar k_k1("k_k1","k_k1",	-0.3,	-1.,	1.) ;
			// RooChebychev k_bgn("k_bgn","k_bgn",x,RooArgSet(k_k0,k_k1));

			//Construct composite pdf
			Int_t k_ent = h[6+cc][5][p][t]->GetEntries(); //[p][t]
			RooRealVar k_N_s("k_N_s","N_a_s",0.1*k_ent,0.,1.05*k_ent) ;
			RooRealVar k_N_b("k_N_b","N_a_b",0.99*k_ent,0.,1.05*k_ent) ;

			// RooRealVar N2("N2","N2",N_k_s.getVal()) ;
			RooRealVar b2mean("b2mean","mean",850,550,1050,"MeV");
			RooRealVar b2sigma("b2sigma","sigma",34,15,50,"MeV");
			RooRealVar b2a("b2a","b2a",1,-5,5,"");
			RooRealVar b2n("b2n","b2n",0,0,50,"");
			// RooCBShape b2gauss("b2gauss","gauss",x,b2mean,b2sigma,b2a,b2n);
			RooGaussian b2gauss("b2gauss","gauss",x,b2mean,b2sigma);


			// RooAddPdf model_ks("model_ks","model_ks",RooArgList(k_sig,k_bgn),RooArgList(k_N_s,k_N_b)) ;


			/***************/
			/*  Rho fit    */
			/***************/

			RooRealVar mean("mean","mean",770,720,820,"MeV") ;
			// RooRealVar sigma1("sigma1","sigma1",40.,0.,100.,"MeV");
			RooRealVar sigma1("sigma1","sigma1",42.,0.,100.,"MeV");
			// RooRealVar sigma2("sigma2","sigma2",5.,0.,100.,"MeV");
			RooRealVar sigma2("sigma2","sigma2",70.,0.,100.,"MeV");

			RooGaussian gauss1("gauss1","gauss1",x,mean,sigma1) ;
			RooGaussian gauss2("gauss2","gauss2",x,mean,sigma2) ;

			//Background
			RooRealVar k0a("k0a","k0a",	-0.2,	-1.,	1.) ;
			RooRealVar k1a("k1a","k1a",	-0.2,	-1.,	1.) ;
			RooChebychev bgna("bgna","bgna",x,RooArgSet(k0a,k1a));//,k3a,k4a)) ;



			RooRealVar k0k("k0k","k0k",	-0.2,	-1.,	1.) ;
			RooRealVar k1k("k1k","k1k",	-0.2,	-1.,	1.) ;
			RooChebychev bgnk("bgnk","bgnk",x,RooArgSet(k0k,k1k));//,k2k));//,k3k,k4k)) ;

			RooRealVar k0p("k0p","k0p",	-0.2,	-1.,	1.) ;
			RooRealVar k1p("k1p","k1p",	-0.2,	-1.,	1.) ;
			// RooRealVar k2p("k2p","k2p",	-0.2,	-1.,	1.) ;
			RooChebychev bgnp("bgnp","bgnp",x,RooArgSet(k0p,k1p));//,k3p,k4p));

			RooRealVar bmean("bmean","mean",650,550,950,"MeV");
			RooRealVar bsigma("bsigma","sigma",34,15,50,"MeV");
			RooRealVar ba("ba","ba",1,-5,5,"");
			RooRealVar bn("bn","bn",4,0,50,"");
			RooGaussian bgauss("bgauss","gauss",x,bmean,bsigma);
			// RooCBShape bgauss("bgauss","gauss",x,bmean,bsigma,ba,bn);


			//Construct composite pdf

			Int_t ent = h[6+cc][0][p][t]->GetEntries(); //[p][t]
			RooRealVar N_a_b("N_a_b","N_a_b",0.05*ent,0.,1.05*ent) ;

			ent = h[6+cc][1][p][t]->GetEntries(); //[p][t]
			RooRealVar N_pi_s("N_pi_s","N_pi_s",0.9*ent,0.,1.05*ent) ;
			RooRealVar N_pi_b("N_pi_b","N_pi_b",0.1*ent,0.,1.05*ent) ;

			ent = h[6+cc][2][p][t]->GetEntries(); //[p][t]
			RooRealVar N_k_s("N_k_s","N_k_s",0.2*ent,0.,1.05*ent) ;
			RooRealVar N_k_b("N_k_b","N_k_b",0.8*ent,0.,1.05*ent) ;

			ent = h[6+cc][3][p][t]->GetEntries(); //[p][t]
			RooRealVar N_p_s("N_p_s","N_p_s",0.2*ent,0.,1.05*ent) ;
			RooRealVar N_p_b("N_p_b","N_p_b",0.8*ent,0.,1.05*ent) ;

			ent = h[6+cc][4][p][t]->GetEntries(); //[p][t]
			RooRealVar N_u_s("N_u_s","N_u_s",0.8*ent,0.,1.05*ent) ;
			RooRealVar N_u_b("N_u_b","N_u_b",0.05*ent,0.,1.05*ent) ;









			if(cc == 1 && t==2 && p >Np-4){
				N_p_s.setVal(0.);
			}

			RooRealVar frac_s("frac_s","frac_s",0.65,0.5,0.7); //
			RooAddPdf sig("sig","sig",RooArgList(gauss1,gauss2),frac_s);

			RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));

			RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgna),RooArgList(N_a_s,N_a_b));
			RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgna),RooArgList(N_pi_s,N_pi_b)); // checked with bgnpi
			RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgna,bgauss),RooArgList(N_k_s,N_k_b, k_N_s)); // checked with bgnk
			// RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgnk),RooArgList(N_k_s,N_k_b)); // checked with bgnk
			RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgna),RooArgList(N_p_s,N_p_b));
			RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgna),RooArgList(N_u_s,N_u_b));

			RooAddPdf model_ks("model_ks","model_ks",RooArgList(k_sig,k_bgn,b2gauss),RooArgList(k_N_s,k_N_b,N_k_s));
			// RooAddPdf model_ks("model_ks","model_ks",RooArgList(k_sig,k_bgn),RooArgList(k_N_s,k_N_b));

			RooCategory sample("sample","sample") ;
			sample.defineType("all") ;
			sample.defineType("pi") ;
			sample.defineType("k") ;
			sample.defineType("p") ;
			sample.defineType("unk") ;
			sample.defineType("ks") ;

			RooDataHist combData("combData","combined data",x,Index(sample),
				Import("all",*h[6+cc][0][p][t]),
				Import("pi", *h[6+cc][1][p][t]),
				Import("k",  *h[6+cc][2][p][t]),
				Import("p",  *h[6+cc][3][p][t]),
				Import("unk",*h[6+cc][4][p][t]),
				Import("ks", *h[6+cc][5][p][t])
			);


			RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
			simPdf.addPdf(model_all,"all");
			simPdf.addPdf(model_pi,"pi");
			simPdf.addPdf(model_k,"k");
			simPdf.addPdf(model_p,"p");
			simPdf.addPdf(model_unk,"unk");
			simPdf.addPdf(model_ks,"ks");


			if(!use_sidebins){

				RooAbsReal* nll = simPdf.createNLL(combData,Extended(true),Range("fitRange"),SplitRange()) ;
				RooMinuit minu(*nll) ;
				minu.setPrintLevel(-1);
				minu.setNoWarn();

				int i = 0;
				do{
					if(i>0){
						N_pi_s.setVal(1.01*N_pi_s.getVal());
						N_k_s.setVal( 0.98*N_k_s.getVal());
						N_p_s.setVal( 0.98*N_p_s.getVal());
						N_u_s.setVal( 0.98*N_u_s.getVal());
						k_N_s.setVal( 0.98*k_N_s.getVal());
					}

					minu.hesse();
					minu.simplex();
					minu.migrad();
					minu.improve();
					minu.hesse();

					r[6+cc][p][t] = minu.save();

					i++;
				}while( r[6+cc][p][t]->covQual() != 3 && i<= retry);
				delete nll;
			}else {
				N_id[6+cc][0][p][t] =   h[6+cc][0][p][t]->Integral(33,85) - h[6+cc][0][p][t]->Integral(3,29) - h[6+cc][0][p][t]->Integral(89,115);
				N_id[6+cc][1][p][t] = ( h[6+cc][1][p][t]->Integral(33,85) - h[6+cc][1][p][t]->Integral(3,29) - h[6+cc][1][p][t]->Integral(89,115) ) /N_id[6+cc][0][p][t];
				N_id[6+cc][2][p][t] = ( h[6+cc][2][p][t]->Integral(33,85) - h[6+cc][2][p][t]->Integral(3,29) - h[6+cc][2][p][t]->Integral(89,115) ) /N_id[6+cc][0][p][t];
				N_id[6+cc][3][p][t] = ( h[6+cc][3][p][t]->Integral(33,85) - h[6+cc][3][p][t]->Integral(3,29) - h[6+cc][3][p][t]->Integral(89,115) ) /N_id[6+cc][0][p][t];
				N_id[6+cc][4][p][t] = ( h[6+cc][4][p][t]->Integral(33,85) - h[6+cc][4][p][t]->Integral(3,29) - h[6+cc][4][p][t]->Integral(89,115) ) /N_id[6+cc][0][p][t];
			}
			//
			//
			r[6+cc][p][t]->Print();
			//
			// N_id[6+cc][0][p][t] = N_a_s.getVal();
			// N_id[6+cc][1][p][t] = N_pi_s.getVal();
			// N_id[6+cc][2][p][t] = N_k_s.getVal() ;
			// N_id[6+cc][3][p][t] = N_p_s.getVal() ;
			// N_id[6+cc][4][p][t] = N_u_s.getVal() ;

			TGaxis::SetMaxDigits(4);
			stringstream nn;
			nn.str("");
			nn << "all: " << p_bins[p]<< " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
			RooPlot* frame1 = x.frame(Title(nn.str().c_str()),Range("fitRange_all"));
			combData.plotOn(frame1,Cut("sample==sample::all")) ;
// 			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame1,Range("fitRange_all"),Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame1,Range("fitRange_all"),Slice(sample,"all"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame1,Range("fitRange_all"),Slice(sample,"all"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame1,Range("fitRange_all"),Slice(sample,"all"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame1,Cut("sample==sample::all")) ;

			RooPlot* frame2 = x.frame(Title("pi"),Range("fitRange_pi")) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;
// 			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame2,Range("fitRange_pi"),Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame2,Range("fitRange_pi"),Slice(sample,"pi"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame2,Range("fitRange_pi"),Slice(sample,"pi"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame2,Range("fitRange_pi"),Slice(sample,"pi"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame2,Cut("sample==sample::pi")) ;

			RooPlot* frame3 = x.frame(Title("k"),Range("fitRange_k")) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;
// 			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),Components("bgauss"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGray)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			simPdf.plotOn(frame3,Range("fitRange_k"),Slice(sample,"k"),Components("bgauss"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGray)) ;
			combData.plotOn(frame3,Cut("sample==sample::k")) ;

			RooPlot* frame4 = x.frame(Title("p"),Range("fitRange_p")) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;
// 			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame4,Range("fitRange_p"),Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame4,Range("fitRange_p"),Slice(sample,"p"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame4,Range("fitRange_p"),Slice(sample,"p"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame4,Range("fitRange_p"),Slice(sample,"p"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame4,Cut("sample==sample::p")) ;

			RooPlot* frame5 = x.frame(Title("noID"),Range("fitRange_unk")) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;
// 			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
			simPdf.plotOn(frame5,Range("fitRange_unk"),Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame5,Range("fitRange_unk"),Slice(sample,"unk"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame5,Range("fitRange_unk"),Slice(sample,"unk"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame5,Range("fitRange_unk"),Slice(sample,"unk"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			combData.plotOn(frame5,Cut("sample==sample::unk")) ;


			RooPlot* frame6 = x.frame(Title("K^{*}"),Range("fitRange_ks")) ;
			combData.plotOn(frame6,Cut("sample==sample::ks")) ;
			simPdf.plotOn(frame6,Range("fitRange_ks"),Slice(sample,"ks"),ProjWData(sample,combData),LineWidth(lw)) ;
			simPdf.plotOn(frame6,Range("fitRange_ks"),Slice(sample,"ks"),Components("k_gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
			simPdf.plotOn(frame6,Range("fitRange_ks"),Slice(sample,"ks"),Components("k_gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
			simPdf.plotOn(frame6,Range("fitRange_ks"),Slice(sample,"ks"),Components("k_bgn"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
			simPdf.plotOn(frame6,Range("fitRange_ks"),Slice(sample,"ks"),Components("b2gauss"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGray)) ;
			combData.plotOn(frame6,Cut("sample==sample::ks")) ;


			TCanvas* c = new TCanvas("c","",4000,600) ;

			frame1->GetXaxis()->SetTitleOffset(1.01);
			frame1->GetXaxis()->SetTitleSize(0.045);
			frame1->GetXaxis()->SetLabelSize(0.045);
			frame1->GetXaxis()->SetNdivisions(505);
			frame1->GetYaxis()->SetTitleOffset(1.3);
			frame1->GetYaxis()->SetTitleSize(0.045);
			frame1->GetYaxis()->SetLabelSize(0.045);
			frame1->GetYaxis()->SetNdivisions(505);

			frame2->GetXaxis()->SetTitleOffset(1.01);
			frame2->GetXaxis()->SetTitleSize(0.045);
			frame2->GetXaxis()->SetLabelSize(0.045);
			frame2->GetXaxis()->SetNdivisions(505);
			frame2->GetYaxis()->SetTitleOffset(1.3);
			frame2->GetYaxis()->SetTitleSize(0.045);
			frame2->GetYaxis()->SetLabelSize(0.045);
			frame2->GetYaxis()->SetNdivisions(505);

			frame3->GetXaxis()->SetTitleOffset(1.01);
			frame3->GetXaxis()->SetTitleSize(0.045);
			frame3->GetXaxis()->SetLabelSize(0.045);
			frame3->GetXaxis()->SetNdivisions(505);
			frame3->GetYaxis()->SetTitleOffset(1.3);
			frame3->GetYaxis()->SetTitleSize(0.045);
			frame3->GetYaxis()->SetLabelSize(0.045);
			frame3->GetYaxis()->SetNdivisions(505);

			frame4->GetXaxis()->SetTitleOffset(1.01);
			frame4->GetXaxis()->SetTitleSize(0.045);
			frame4->GetXaxis()->SetLabelSize(0.045);
			frame4->GetXaxis()->SetNdivisions(505);
			frame4->GetYaxis()->SetTitleOffset(1.3);
			frame4->GetYaxis()->SetTitleSize(0.045);
			frame4->GetYaxis()->SetLabelSize(0.045);
			frame4->GetYaxis()->SetNdivisions(505);

			frame5->GetXaxis()->SetTitleOffset(1.01);
			frame5->GetXaxis()->SetTitleSize(0.045);
			frame5->GetXaxis()->SetLabelSize(0.045);
			frame5->GetXaxis()->SetNdivisions(505);
			frame5->GetYaxis()->SetTitleOffset(1.3);
			frame5->GetYaxis()->SetTitleSize(0.045);
			frame5->GetYaxis()->SetLabelSize(0.045);
			frame5->GetYaxis()->SetNdivisions(505);

			frame6->GetXaxis()->SetTitleOffset(1.01);
			frame6->GetXaxis()->SetTitleSize(0.045);
			frame6->GetXaxis()->SetLabelSize(0.045);
			frame6->GetXaxis()->SetNdivisions(505);
			frame6->GetYaxis()->SetTitleOffset(1.3);
			frame6->GetYaxis()->SetTitleSize(0.045);
			frame6->GetYaxis()->SetLabelSize(0.045);
			frame6->GetYaxis()->SetNdivisions(505);

			c->Clear();
			TPad* pa[6];
			pa[0] = new TPad("p1","",0./6.,	0.,	1./6. ,1.);
			pa[1] = new TPad("p2","",1./6.,	0.,	2./6. ,1.);
			pa[2] = new TPad("p3","",2./6.,	0.,	3./6. ,1.);
			pa[3] = new TPad("p4","",3./6.,	0.,	4./6. ,1.);
			pa[4] = new TPad("p5","",4./6.,	0.,	5./6. ,1.);
			pa[5] = new TPad("p6","",5./6.,	0.,	6./6. ,1.);
			for(int dr = 0; dr<6;dr++){
				pa[dr]->SetLeftMargin(0.12);
				pa[dr]->SetRightMargin(0.08);
				pa[dr]->Draw();
			}
			pa[0]->cd() ; frame1->Draw() ;
			// TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
			// bb->SetY1(0.3);
			pa[1]->cd() ; frame2->Draw() ;
			pa[2]->cd() ; frame3->Draw() ;
			pa[3]->cd() ; frame4->Draw() ;
			pa[4]->cd() ; frame5->Draw() ;
			pa[5]->cd() ; frame6->Draw() ;

			c->Print("test.root");

			if(p==0 && t == 0){
			// if(p== 10 && t == 1){
				if(cc == 0) c->Print("test_rho0_pim.pdf[");
				else if(cc == 1) c->Print("test_rho0_pip.pdf[");
			}

			if(cc == 0) c->Print("test_rho0_pim.pdf");
			if(cc == 1) c->Print("test_rho0_pip.pdf");


			if(p==Np-1 && t == Nt-1){
			// if(p== 6 && t == Nt-1){
			// if(p== 11 && t == 2){
				if(cc == 0) c->Print("test_rho0_pim.pdf]");
				else if(cc == 1) c->Print("test_rho0_pip.pdf]");
			}
			delete frame1;
			delete frame2;
			delete frame3;
			delete frame4;
			delete frame5;
			delete frame6;
			for(int ll=0;ll<6;ll++) delete pa[ll];
			delete c;
		}
	}
}

/**********************************************************************/
void print_table(){
	return;
	if(start<=0 && stop >2){
		input_k0->Close();
		delete input_k0;
	}else if(start<=2 && stop >4){
		input_lam->Close();
		delete input_lam;
	}else if(start<=4 && stop >6){
		input_phi->Close();
		delete input_phi;
	}else if(start<=6 && stop >8){
		input_rho->Close();
		delete input_rho;
	}
	TGraphErrors* gr[8][4][Nt];

	TGraph* grC[8][Nt];

	double ppp[Np];
	double aa[Np];
	stringstream nn;

	for(int i = 0; i< Np; i++){
		ppp[i] = p_bins[i]/2. + p_bins[i+1]/2.;
		aa[i] =0.;
	}

	TFile* ff= new TFile(out_file.c_str(),"RECREATE");
	string p_name[8] = {"pi_m","pi_p","k_m","k_p","p_m","p_p","rpi_m","rpi_p"};
	string id[4] = {"pi","k","p","u"};
	string id2[4] = {"pi","k","p","rpi"};
	Int_t color[4] = {kOrange+7,kAzure+4,kTeal+4,kYellow+2};

	ofstream ofs_matrix("rich_mat.txt", std::ofstream::out | std::ofstream::trunc);

	int color2[7] = {kRed, kOrange+7, kYellow+2, kSpring-6, kCyan-6, kAzure+3, kViolet-1};
	double shift[4] = {0.,0.1,0.2,0.3};

	int cov_elem[4] = {6,2,4,8}; // pi,k,p,u

	for(int p = 0; p< Np; p++)
	{
  	for(int t = 1; t< 3; t++)
		{
				ofs_matrix << p << "\t" << t;
				for(int i = start; i<stop; i++)
				{
					for(int j = 1; j<4; j++)
					{
						double val;
						double aaa = N_id[i][j][p][t];
						double ggg = N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t];
						val = (ggg ? aaa/ggg : 0);
						ofs_matrix << "\t" << val;
					}
				}
				ofs_matrix << endl;
		}
	}

	ofs_matrix.close();

	for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
		for(int t = 0; t< Nt; t++){

			if(t!=1 && t!=2) continue;

			double diff[Np];
			for(int p = 0; p< Np; p++){
				if(t==3 && p>6) continue;

				if(p!= 10 && p!=11) continue;


// 				diff[p] = (N_id[i][1][p][t] + N_id[i][2][p][t] + N_id[i][3][p][t] +N_id[i][4][p][t]) + shift[t];
				if(r[i][p][t]) diff[p] = (r[i][p][t]->covQual()/3.) + shift[t];
				else diff[p] = -1;
			}
			if(t!=3)grC[i][t] = new TGraph(Np,ppp,diff);
			else if(t==3)grC[i][t] = new TGraph(7,ppp,diff);
			grC[i][t]->SetMarkerStyle(33);
			grC[i][t]->SetMarkerSize(2);
			grC[i][t]->SetMarkerColor(color2[t]);
			grC[i][t]->SetLineColor(color2[t]);
		}
		for(int j = 1; j<5; j++) {  // identified as (pi,k,p)
			for(int t = 0; t< Nt; t++){
// 			for(int t = 1; t< 2; t++){

				if(t!=1 &&t!=2) continue;

				double val[Np],err[Np];
				for(int p = 0; p< Np; p++){

					if(p!= 10 && p!=11) continue;

// 				for(int p = 0; p< 2; p++){
					if(t==3 && p>6) continue;
					double aaa = N_id[i][j][p][t];
					double ggg = N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t];
					val[p] = aaa/ggg;

					double tmp = 0.;

					double tmp_jak[4];

					for(int ll=0; ll<4;ll++){
						tmp_jak[ll] = -aaa/(ggg*ggg);
						if(ll == j-1) tmp_jak[ll] += 1./ggg;
					}

					for(int ll = 0; ll<4;ll++) tmp += tmp_jak[ll]*tmp_jak[ll]*r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[ll]);

					for(int ll = 0; ll<3; ll++){
						for(int hh = ll+1; hh<4; hh++){
							tmp += tmp_jak[ll]*tmp_jak[hh]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[hh]);
						}
					}

					err[p] = TMath::Sqrt(tmp);
					cout << i << " " << j << "    " << t << " " << p << "  " << val[p] << " +/- " << err[p] << endl;
					// err[p] = TMath::Sqrt((aaa+1)*(ggg-aaa+1)/( (ggg+2)*(ggg+2)*(ggg+3) ));
				}
				if(t != 3)gr[i][j-1][t] = new TGraphErrors(Np,ppp,val,aa,err);
				else if(t==3)gr[i][j-1][t] = new TGraphErrors(7,ppp,val,aa,err);
				gr[i][j-1][t]->SetMarkerStyle(33);
				gr[i][j-1][t]->SetMarkerColor(color[j-1]);
				gr[i][j-1][t]->SetLineColor(color[j-1]);

				nn.str("");
				nn.clear();
				nn << p_name[i] << "_" << id[j-1] << "_" << t;
				gr[i][j-1][t]->SetName(nn.str().c_str());

				gr[i][j-1][t]->Write();
			}
		}
	}

	for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
		for(int t = 0; t< Nt; t++){

			if(t!=1 && t!=2) continue;

			for(int p = 0; p< Np; p++){
				if(t==3 && p>6) continue;

				if(p!= 10 && p!=11) continue;

				nn.str("");
				nn.clear();
				nn <<"fit_" << p_name[i] << "_" << t << "_" << p;
				r[i][p][t]->Write(nn.str().c_str());
				delete r[i][p][t];
			}
		}
	}
	gStyle->SetOptStat(0);

	for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
		for(int j = 1; j<5; j++) {  // identified as (pi,k,p)
			for(int t = 0; t< Nt; t++){

				if(t!=1 && t!=2) continue;

				gr[i][j-1][t]->SetMarkerColor(color2[t]);
				gr[i][j-1][t]->SetLineColor(color2[t]);
			}
		}
	}

	TCanvas* c = new TCanvas("c","",800,600);
	TLegend *leg = new TLegend(0.13,0.25,0.33,0.55);
	leg->SetLineColor(0);
	leg->SetNColumns(1);
	leg->SetFillStyle(0);

	for(int t = 0; t< Nt; t++){

		if(t!=1 && t!=2) continue;

		nn.str("");
		nn << t_bins[t] <<"#leq #theta < " << t_bins[t+1];
		leg->AddEntry(gr[start][0][t],nn.str().c_str(),"pl");
	}

	string lable[8] = {"#pi^{-}","#pi^{+}","K^{-}","K^{+}","#bar{p}","p","unknown","unknown"};

	double min[4][4];
	double max[4][4];


	max[0][0]=1.0;
	max[0][1]=1.0;
	max[0][2]=1.0;
	max[0][3]=1.0;

	max[1][0]=1.0;
	max[1][1]=1.0;
	max[1][2]=1.0;
	max[1][3]=1.0;

	max[2][0]=1.0;
	max[2][1]=1.0;
	max[2][2]=1.0;
	max[2][3]=1.0;

	max[3][0]=1.0;
	max[3][1]=1.0;
	max[3][2]=1.0;
	max[3][3]=1.0;

	min[0][0]=0.;
	min[0][1]=0.;
	min[0][2]=0.;
	min[0][3]=0.;

	min[1][0]=0.;
	min[1][1]=0.;
	min[1][2]=0.;
	min[1][3]=0.;

	min[2][0]=0.;
	min[2][1]=0.;
	min[2][2]=0.;
	min[2][3]=0.;

	min[3][0]=0.;
	min[3][1]=0.;
	min[3][2]=0.;
	min[3][3]=0.;

	TLine* line = new TLine(0,1,50,1);
	line->SetLineStyle(2);
	line->SetLineColor(kGray+1);

	TLine* li[4];
// 	li[0] = new TLine(0,0.995,50,0.995);
// 	li[1] = new TLine(0,1.000,50,1.000);
// 	li[2] = new TLine(0,1.005,50,1.005);
	li[0] = new TLine(0,1.0,50,1.0);
	li[1] = new TLine(0,1.1,50,1.1);
	li[2] = new TLine(0,1.2,50,1.2);
	li[3] = new TLine(0,1.3,50,1.3);
	for(int l = 0 ; l<4; l++){
		li[l]->SetLineStyle(2);
		li[l]->SetLineColor(kGray+1);
	}

	for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
		for(int j = 1; j<5; j++) {  // identified as (pi,k,p)
			c->Clear();
			TH1F* hr = c->DrawFrame(0.,min[i/2][j-1],50.,max[i/2][j-1]);
			hr->GetXaxis()->SetTitle("p (GeV/c)");
			nn.str("");
			nn << lable[i] <<" #rightarrow " << lable[2*(j-1)+i%2];//;id2[j-1];
			hr->SetTitle(nn.str().c_str());
			leg->Draw();
			for(int t = 0 ; t< Nt; t++){

				if(t!=1 && t!=2) continue;

				gr[i][j-1][t]->Draw("pl0");
			}
			nn.str("");
			if(i%2 == 0) nn << "table/" << id2[i/2] <<"/" << id[i/2] << "m_" << id[j-1] <<".pdf";
			else nn << "table/" << id2[i/2] <<"/" << id[i/2] << "p_" << id[j-1] <<".pdf";
			c->Print(nn.str().c_str());
			for(int t = 0 ; t< Nt; t++){
				if(t!=1 && t!=2) continue;
				delete gr[i][j-1][t];
			}
			delete hr;
		}
	}


	c->Print("check_sum.pdf[");
	for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
		c->Clear();
// 		TH1F* hr = c->DrawFrame(0.,0.98,50.,1.02);
		TH1F* hr = c->DrawFrame(0.,0.,50.,1.37);
		hr->GetXaxis()->SetTitle("p (GeV/c)");
		hr->GetXaxis()->SetLabelSize(0.045);
		hr->GetXaxis()->SetTitleSize(0.045);
		hr->GetXaxis()->SetTitleOffset(0.9);
		hr->GetYaxis()->SetLabelSize(0.045);
		nn.str("");
		nn << "Check " << lable[i];
		hr->SetTitle(nn.str().c_str());
		for(int t = 0 ; t< Nt; t++){

			if(t!=1 && t!=2) continue;

			li[t]->Draw("same");
			grC[i][t]->Draw("pl");
		}
		c->Print("check_sum.pdf");

		for(int t = 0 ; t< Nt; t++){
			if(t!=1 && t!=2) continue;
			delete grC[i][t];
		}

		delete hr;
	}
	c->Print("check_sum.pdf]");


	/*
	string charge[2] = {"n","p"};
	for(int c = 0; c<2;c++){
		for(int t = 0; t<Nt;t++){
			for(int p=0; p<Np;p++){
				TMatrixD eff3(3,3);
				nn.str("");
				nn << "eff3_"<<charge[c] <<"_t"<<t<<"_p"<<TMath::Nint(p_bins[p]);
				//0 pi
				//1 K
				//2 p

				eff3[0][2] = N_id[0+c][3][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
				eff3[0][1] = N_id[0+c][2][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
				eff3[0][0] = N_id[0+c][1][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);

				eff3[1][2] = N_id[2+c][3][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
				eff3[1][1] = N_id[2+c][2][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
				eff3[1][0] = N_id[2+c][1][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);

				eff3[2][2] = N_id[4+c][3][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff3[2][1] = N_id[4+c][2][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff3[2][0] = N_id[4+c][1][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);

				eff3.Write(nn.str().c_str());
				TMatrixD eff4(3,4);
				nn.str("");
				nn << "eff4_"<<charge[c] <<"_t"<<t<<"_p"<<TMath::Nint(p_bins[p]);
				//0 pi
				//1 K
				//2 p
				//3 u

				eff4[0][3] = N_id[0+c][4][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
				eff4[0][2] = N_id[0+c][3][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
				eff4[0][1] = N_id[0+c][2][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
				eff4[0][0] = N_id[0+c][1][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);

				eff4[1][3] = N_id[2+c][4][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
				eff4[1][2] = N_id[2+c][3][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
				eff4[1][1] = N_id[2+c][2][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
				eff4[1][0] = N_id[2+c][1][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);

				eff4[2][3] = N_id[4+c][4][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff4[2][2] = N_id[4+c][3][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff4[2][1] = N_id[4+c][2][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff4[2][0] = N_id[4+c][1][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
				eff4.Write(nn.str().c_str());
			}
		}
	}
	*/


	ff->Close();


	delete leg;
	delete line;
	for(int i=0;i<4;i++)delete li[i];
	delete c;
	delete ff;
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
