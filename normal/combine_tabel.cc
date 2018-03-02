void combine_tabel(){
	TFile* f[3];
	f[0] = new TFile("rich_pi.root","read");
	f[1] = new TFile("rich_k.root","read");
	f[2] = new TFile("rich_p.root","read");
	
	TFile* f_out = new TFile("rich_2011_more_theta.root","recreate");

	string p_name[6] = {"pi_m","pi_p","k_m","k_p","p_m","p_p"};
	string id[4] = {"pi","k","p","u"};

	TGraphErrors* gr;

	stringstream nn;

	for(int i = 0; i<6; i++) {  // particle (pi,k,p) 0,6
		for(int j = 1; j<5; j++) {  // identified as (pi,k,p)
			for(int t = 0; t< 4; t++){
				nn.str("");
				nn.clear();
				nn << p_name[i] << "_" << id[j-1] << "_" << t;
				((TGraphErrors*)f[i/2]->Get(nn.str().c_str()))->Write();
			}	
		}
	}

	f_out->Close();
	f[0]->Close();
	f[1]->Close();
	f[2]->Close();


}
