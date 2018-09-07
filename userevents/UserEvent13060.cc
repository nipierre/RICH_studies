#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <sstream>
#include "TH1.h"
#include "TH2.h"
#include "G3part.h"
#include "TROOT.h"

#include "Phast.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "PaParticle.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaDigit.h"
#include "PaAlgo.h"
#include "PaPid.h"
#include "PaMetaDB.h"

#include "TargetCell.h"


/****************************/
/* RICH Matrix				*/
/* inclusive phi 			*/
/* selection				*/
/****************************/


enum trigger{
	IT 		= 1 	<<0,
	MT 		= 1  	<<1,
	LT 		= 1 	<<2,
	OT		= 1 	<<3,
	CT		= 1 	<<4,
	IV 		= 1		<<5,
	HaloT	= 1 	<<6,
	BT 		= 1 	<<7,
	MTinc 	= 1 	<<8,
	LAST 	= 1		<<9,
	NRand 	= 1 	<<10,
	TRand 	= 1 	<<11
};



void UserEvent13060(PaEvent& e)
{
	static TTree* tree2(NULL);
	static Long64_t Evt;
	static int Run, TriggerMask;
	static float vx, vy ,vz;
	static double Q2, xbj, y;

	static TLorentzVector lv_beam, lv_scat, lv_gamma;
	static TLorentzVector lv_kp, lv_km, lv_phi;
	static double Emiss;
	static double pt1, pt2, alpha;
	static double pi_thr, k_thr, p_thr;

	const float m_pi = G3partMass[8];
	const float m_p = G3partMass[14];
	const float m_k = G3partMass[11];
	const float m_phi = G3partMass[91];
//	const float m_lambda = G3partMass[18];
//	const float m_k0 = G3partMass[10];
	const float m_mu = G3partMass[5];

	static double pp_lh[7],pm_lh[7];

	for(int i = 0; i<7; i++) {
		pm_lh[i] = -1;
		pp_lh[i] = -1;
	}

	static double pp_x, pp_y, pp_mom, pp_theta, pm_x, pm_y, pm_mom, pm_theta;
	static int Nout,charge1, charge2;
//	static TH2D* h2;

	TargetCell* fTcell;
	fTcell = new TargetCell();

	static bool first(true);
	if(first){

		tree2 = new TTree("tree2","tree2");
		tree2->Branch("Run",    			&Run,    			"Run/I");
		tree2->Branch("Evt",    			&Evt,    			"Evt/L");
		tree2->Branch("TriggerMask", 		&TriggerMask, 		"TriggerMask/I");
		tree2->Branch("lv_beam",			"TLorentzVector", 	&lv_beam);
		tree2->Branch("lv_scat",			"TLorentzVector", 	&lv_scat);
		tree2->Branch("Q2",					&Q2,				"Q2/D");
		tree2->Branch("xbj",				&xbj,				"xbj/D");
		tree2->Branch("y",					&y,					"y/D");
		tree2->Branch("Emiss",				&Emiss,				"Emiss/D");
		tree2->Branch("vx",					&vx,				"vx/F");
		tree2->Branch("vy",					&vy,				"vy/F");
		tree2->Branch("vz",					&vz,				"vz/F");
		tree2->Branch("lv_kp",				"TLorentzVector", 	&lv_kp);
		tree2->Branch("lv_km",				"TLorentzVector", 	&lv_km);
		tree2->Branch("lv_phi",				"TLorentzVector", 	&lv_phi);
		tree2->Branch("pt1",				&pt1,				"pt1/D");
		tree2->Branch("pt2",				&pt2,				"pt2/D");
		tree2->Branch("alpha",				&alpha,				"alpha/D");
		tree2->Branch("pp_x",				&pp_x,				"pp_x/D");
		tree2->Branch("pp_y",				&pp_y,				"pp_y/D");
		tree2->Branch("pp_mom",				&pp_mom,			"pp_mom/D");
		tree2->Branch("pp_theta",			&pp_theta,			"pp_theta/D");
		tree2->Branch("pm_x",				&pm_x,				"pm_x/D");
		tree2->Branch("pm_y",				&pm_y,				"pm_y/D");
		tree2->Branch("pm_mom",				&pm_mom,			"pm_mom/D");
		tree2->Branch("pm_theta",			&pm_theta,			"pm_theta/D");
		tree2->Branch("pm_lh",				pm_lh,				"pm_lh[7]/D");
		tree2->Branch("pp_lh",				pp_lh,				"pp_lh[7]/D");
		tree2->Branch("k_thr",				&k_thr,				"k_thr/D");
		tree2->Branch("pi_thr",				&pi_thr,			"pi_thr/D");
		tree2->Branch("p_thr",				&p_thr,				"p_thr/D");
		tree2->Branch("Nout", 				&Nout,		 		"Nout/I");
		tree2->Branch("charge1", 			&charge1,		 	"charge1/I");
		tree2->Branch("charge2", 			&charge2,		 	"charge2/I");

		first=false;


	}

	double n_rich = 1 + PaMetaDB::Ref().NminusOne(Run)/1e6;
	if (PaMetaDB::Ref().NminusOne(Run) < 0) n_rich = 1 +  PaSetup::Ref().NminusOne()/1e6;
	k_thr = (m_k/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));
	pi_thr = (m_pi/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));
	p_thr = (m_p/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));

	Run        =  e.RunNum();
	Evt        =  e.UniqueEvNum();
	TriggerMask = (e.TrigMask() & 0xffff);

	const PaSetup& setup = PaSetup::Ref();
	const double z_rich = setup.Rich().DetPos(0).Z();// - 1.65;

	PaPid pid;

	int Bvtx   =  e.iBestPrimaryVertex();
// 	cout << 1 << endl;
	for(int iv = 0; iv < e.NVertex(); iv++){ // loop over reconstructed vertices
		const PaVertex& v = e.vVertex(iv);
		if(! v.IsPrimary()) continue; // not primary. Skip.
		if(Bvtx !=-1 && Bvtx != iv) continue;

		vx = v.X();
		vy = v.Y();
		vz = v.Z();

		int imu0 = v.InParticle();
		int imu1 = v.iMuPrim();
		if(imu0 == -1) continue;
		if(imu1 == -1) continue;
// 		cout << 2<< endl;
		const PaParticle& pa_beam = e.vParticle(imu0);
		const PaParticle& pa_scat = e.vParticle(imu1);

		const PaTPar& par_beam = pa_beam.ParInVtx(iv);
		const PaTPar& par_scat = pa_scat.ParInVtx(iv);

    if( !fTcell->TargetCell::CrossCells(par_beam, Run)) continue;
		lv_beam = par_beam.LzVec(m_mu);
		lv_scat = par_scat.LzVec(m_mu);
		lv_gamma = lv_beam - lv_scat;

		Q2 = PaAlgo::Q2(lv_beam,lv_scat);
		xbj = PaAlgo::xbj(lv_beam, lv_scat);
		y = (lv_beam.E()-lv_scat.E())/lv_beam.E();
		if(y<0.1 || y>0.9) continue;
// 		cout << 3 << endl;

// 		if(v.NOutParticles() !=3 ) continue;
		Nout = v.NOutParticles();
		if(Nout < 3) continue;
		for(int l = 0; l< v.NOutParticles();l++){
			int i_part1 = v.iOutParticle(l);
			if(i_part1 == imu1) continue;
			for(int m = l; m< v.NOutParticles();m++){
				int i_part2 = v.iOutParticle(m);
				if(i_part2 == imu1) continue;
				if(i_part2 == i_part1) continue;


				PaParticle p1 = e.vParticle(i_part1);
				PaParticle p2 = e.vParticle(i_part2);
				if(p1.Q() == p2.Q()) continue;
				if(p1.Q() < p2.Q()) swap(p1,p2);   // now first particle is positive

				charge1 = p1.Q();
				charge2 = p2.Q();

				if( p1.ParInVtx(Bvtx).Mom()<9. || p1.ParInVtx(Bvtx).Mom()>55.) continue;
				if( p2.ParInVtx(Bvtx).Mom()<9. || p2.ParInVtx(Bvtx).Mom()>55.) continue;

				lv_kp = p1.ParInVtx(Bvtx).LzVec(m_k);
				lv_km = p2.ParInVtx(Bvtx).LzVec(m_k);
				lv_phi = lv_kp + lv_km;

				TLorentzVector lv_ppp;
				lv_ppp.SetXYZM(0.,0.,0.,m_p);
				double MX2 = (lv_beam-lv_scat+lv_ppp-lv_phi).Mag2();
				Emiss = (MX2 - m_p*m_p)/(2*m_p);
//				if( TMath::Abs(Emiss) > 2.5 ) continue;

				if(fabs(lv_phi.Mag() - m_phi) >= 0.12 )  continue;

				PaTrack tr_p1= e.vTrack(p1.iTrack());
				PaTrack tr_p2= e.vTrack(p2.iTrack());

				if(tr_p1.ZLast() < 350) continue;
				if(tr_p2.ZLast() < 350) continue;

				if(!tr_p1.HasMom()) continue;
				if(!tr_p2.HasMom()) continue;

				for(int i = 0; i<6; i++) pp_lh[i] = pid.GetLike(i, tr_p1);
				for(int i = 0; i<6; i++) pm_lh[i] = pid.GetLike(i, tr_p2);
				pm_lh[6] = pid.SecondLike(tr_p1);
				pp_lh[6] = pid.SecondLike(tr_p2);

				PaTPar par1, par2;
				if( !tr_p1.Extrapolate(z_rich,par1) ) continue;
				if( !tr_p2.Extrapolate(z_rich,par2) ) continue;;

				const PaTPar& p1_par = p1.ParInVtx(iv);
				const PaTPar& p2_par = p2.ParInVtx(iv);

				TVector3 p1_mom = p1_par.Mom3();
				TVector3 p2_mom = p2_par.Mom3();

				TVector3 p3_mom = p1_mom + p2_mom;
				if(p1_mom.Pt(p3_mom) <0.023) continue;

				pt1 = p1_mom.Pt(p3_mom);
				pt2 = p2_mom.Pt(p3_mom);
				alpha = (TMath::Sqrt(p1_mom.Mag2() - pt1*pt1) - TMath::Sqrt(p2_mom.Mag2() - pt2*pt2))/ (TMath::Sqrt(p1_mom.Mag2() - pt1*pt1) + TMath::Sqrt(p2_mom.Mag2() - pt2*pt2));

				pp_x     = par1.X();
				pp_y     = par1.Y();
				pp_mom   = par1.Mom();
				pp_theta = par1.Theta();
				pm_x     = par2.X();
				pm_y     = par2.Y();
				pm_mom   = par2.Mom();
				pm_theta = par2.Theta();

				tree2->Fill();
			}
		}

	}

}

void UserJobEnd13060(){
	const Phast& ph = Phast::Ref();

	if(ph.print) cout<<"[ UserJobEnd0 has been called ]"<<endl;
	ph.h_file->cd();

}
