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
/* K0, lambda	 			*/
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

void UserEvent13010(PaEvent& e)
{

	static TTree* tree(NULL);
	static TTree* tree2(NULL);
//	static TTree* tree3(NULL);

	static float t_t, t_p, t_x, t_y;
	static double t_n;

	static int Run, Evt, TriggerMask;
	static float vx, vy ,vz, v2x, v2y ,v2z, Chi2;
	static double Q2, xbj, y;

	static TLorentzVector lv_beam, lv_scat, lv_gamma;
	static TLorentzVector lv_p, lv_pi, lv_lambda;
	static TLorentzVector lv_pip, lv_pim, lv_k0;
	static TLorentzVector lv_kp, lv_km, lv_phi;
	static double Emiss;
	static double pt1, pt2, alpha;
	static double pi_thr, k_thr, p_thr;

	const float m_pi = G3partMass[8];
	const float m_p = G3partMass[14];
	const float m_k = G3partMass[11];
	const float m_phi = G3partMass[91];
	const float m_lambda = G3partMass[18];
	const float m_k0 = G3partMass[10];
	const float m_mu = G3partMass[5];

	static double pp_lh[7],pm_lh[7];

	for(int i = 0; i<7; i++) {
		pm_lh[i] = -1;
		pp_lh[i] = -1;
	}

	static double pp_x, pp_y, pp_mom, pp_theta, pm_x, pm_y, pm_mom, pm_theta, pm_ch, pp_ch;

//	static TH2D* h2;

	TargetCell* fTcell;
	fTcell = new TargetCell();

	static bool first(true);
	if(first){

		tree = new TTree("tree","tree");
		tree->Branch("Run",    				&Run,    			"Run/I");
		tree->Branch("Evt",    				&Evt,    			"Evt/I");
		tree->Branch("TriggerMask", 		&TriggerMask, 		"TriggerMask/I");
		tree->Branch("lv_beam",				"TLorentzVector", 	&lv_beam);
		tree->Branch("lv_scat",				"TLorentzVector", 	&lv_scat);
		tree->Branch("Q2",					&Q2,				"Q2/D");
		tree->Branch("xbj",					&xbj,				"xbj/D");
		tree->Branch("y",					&y,					"y/D");
		tree->Branch("vx",					&vx,				"vx/F");
		tree->Branch("vy",					&vy,				"vy/F");
		tree->Branch("vz",					&vz,				"vz/F");
		tree->Branch("v2x",					&v2x,				"v2x/F");
		tree->Branch("v2y",					&v2y,				"v2y/F");
		tree->Branch("v2z",					&v2z,				"v2z/F");
		tree->Branch("Chi2",				&Chi2,				"Chi2/D");
		tree->Branch("lv_p",				"TLorentzVector", 	&lv_p);
		tree->Branch("lv_pi",				"TLorentzVector", 	&lv_pi);
		tree->Branch("lv_lambda",			"TLorentzVector", 	&lv_lambda);

		tree->Branch("lv_pip",				"TLorentzVector", 	&lv_pip);
		tree->Branch("lv_pim",				"TLorentzVector", 	&lv_pim);
		tree->Branch("lv_k0",				"TLorentzVector", 	&lv_k0);
		tree->Branch("pt1",					&pt1,				"pt1/D");
		tree->Branch("pt2",					&pt2,				"pt2/D");
		tree->Branch("alpha",				&alpha,				"alpha/D");
		tree->Branch("pp_x",				&pp_x,				"pp_x/D");
		tree->Branch("pp_y",				&pp_y,				"pp_y/D");
		tree->Branch("pp_mom",				&pp_mom,			"pp_mom/D");
		tree->Branch("pp_ch",				&pp_ch,				"pp_ch/D");
		tree->Branch("pp_theta",			&pp_theta,			"pp_theta/D");
		tree->Branch("pm_x",				&pm_x,				"pm_x/D");
		tree->Branch("pm_y",				&pm_y,				"pm_y/D");
		tree->Branch("pm_mom",				&pm_mom,			"pm_mom/D");
		tree->Branch("pm_ch",				&pm_ch,				"pm_ch/D");
		tree->Branch("pm_theta",			&pm_theta,			"pm_theta/D");
		tree->Branch("pm_lh",				pm_lh,				"pm_lh[7]/D");
		tree->Branch("pp_lh",				pp_lh,				"pp_lh[7]/D");
		tree->Branch("k_thr",				&k_thr,				"k_thr/D");
		tree->Branch("pi_thr",				&pi_thr,			"pi_thr/D");
		tree->Branch("p_thr",				&p_thr,				"p_thr/D");

		/* old exclusive phi selection */
		/*
		tree2 = new TTree("tree2","tree2");
		tree2->Branch("Run",    			&Run,    			"Run/I");
		tree2->Branch("Evt",    			&Evt,    			"Evt/I");
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
		tree2->Branch("pp_ch",				&pp_ch,				"pp_ch/D");
		tree2->Branch("pp_theta",			&pp_theta,			"pp_theta/D");
		tree2->Branch("pm_x",				&pm_x,				"pm_x/D");
		tree2->Branch("pm_y",				&pm_y,				"pm_y/D");
		tree2->Branch("pm_mom",				&pm_mom,			"pm_mom/D");
		tree2->Branch("pm_ch",				&pm_ch,				"pm_ch/D");
		tree2->Branch("pm_theta",			&pm_theta,			"pm_theta/D");
		tree2->Branch("pm_lh",				pm_lh,				"pm_lh[7]/D");
		tree2->Branch("pp_lh",				pp_lh,				"pp_lh[7]/D");
		tree2->Branch("k_thr",				&k_thr,				"k_thr/D");
		tree2->Branch("pi_thr",				&pi_thr,			"pi_thr/D");
		tree2->Branch("p_thr",				&p_thr,				"p_thr/D");
		*/
		/*
		tree3 = new TTree("tree3","tree3");
		tree3->Branch("theta",				&t_t,				"t_t/F");
		tree3->Branch("p",					&t_p,				"t_p/F");
		tree3->Branch("n",					&t_n,				"t_n/D");
		tree3->Branch("x",					&t_x,				"t_x/F");
		tree3->Branch("y",					&t_y,				"t_y/F");
		*/
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

	/*
	for( int pp = 0; pp< e.NTrack();pp++){
		const PaTrack & ttr = e.vTrack(pp);
		if( ttr.ZLast()<350 ) continue;
		if( ttr.ZFirst()>350 ) continue;
		if( ttr.iParticle()==-1 ) continue;
		if( !ttr.HasMom() ) continue;
		PaPid ppii;
		if(!ppii.CheckRichInfo(ttr)) continue;
		PaTPar tpar;
		if( !ttr.Extrapolate(z_rich,tpar) ) continue;
		if(tpar.X()*tpar.X() + tpar.Y()*tpar.Y() < 25) continue;
		t_x = tpar.X();
		t_y = tpar.Y();
		t_t=ppii.Theta_Ch(ttr);
		t_p=tpar.Mom();
		t_n=n_rich;
		if(t_t != 0 && t_p <= 60)tree3->Fill();
	}
	*/
	PaPid pid;

	int Bvtx   =  e.iBestPrimaryVertex();
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
		const PaParticle& pa_beam = e.vParticle(imu0);
		const PaParticle& pa_scat = e.vParticle(imu1);

		const PaTPar& par_beam = pa_beam.ParInVtx(iv);
		const PaTPar& par_scat = pa_scat.ParInVtx(iv);

		//if( !fTcell->TargetCell::InTarget(par_beam,Run) ) continue;  THIS FUNCTION IS FUCKED UP. TODO SEE WHY
    if( !fTcell->TargetCell::CrossCells(par_beam, Run)) continue;
		lv_beam = par_beam.LzVec(m_mu);
		lv_scat = par_scat.LzVec(m_mu);
		lv_gamma = lv_beam - lv_scat;

		Q2 = PaAlgo::Q2(lv_beam,lv_scat);
		xbj = PaAlgo::xbj(lv_beam, lv_scat);
		y = (lv_beam.E()-lv_scat.E())/lv_beam.E();
		if(y<0.1 || y>0.9) continue;
		for(Int_t iv2 = 0; iv2 < e.NVertex(); iv2++){
			const PaVertex& v2 = e.vVertex(iv2);
			if(v2.IsPrimary()) continue;
			if(v2.NOutParticles() != 2) continue;
			v2x    = v2.X();
			v2y    = v2.Y();
			v2z    = v2.Z();

			Chi2  = v2.Chi2();

			int ip1 = v2.iOutParticle(0);
			int ip2 = v2.iOutParticle(1);
			PaParticle p1 = e.vParticle(ip1);
			PaParticle p2 = e.vParticle(ip2);
			if(p1.NFitPar() == 0 || p2.NFitPar() == 0) continue;
			if(p1.Q() == p2.Q())                       continue;
			if(p1.Q() < p2.Q()) swap(p1,p2);

			const PaTrack& tr_p1 = e.vTrack(p1.iTrack());
			const PaTrack& tr_p2 = e.vTrack(p2.iTrack());

			if(tr_p1.XX0() >10) continue;
			if(tr_p2.XX0() >10) continue;

			const PaTPar& p1_par = p1.ParInVtx(iv2);
			const PaTPar& p2_par = p2.ParInVtx(iv2);

			if(p1_par.Mom() < 1) continue;
			if(p2_par.Mom() < 1) continue;

			bool pp = false;
			for(int ipv = 0; ipv < p1.NVertex(); ipv++){
			   if(e.vVertex(p1.iVertex(ipv)).IsPrimary()) pp = true;
			}

			if(pp) continue;

			pp = false;
			for(int ipv = 0; ipv < p2.NVertex(); ipv++){
			   if(e.vVertex(p2.iVertex(ipv)).IsPrimary()) pp = true;
			}

			if(pp) continue;

			if(tr_p1.ZLast() < 350) continue;
			if(tr_p2.ZLast() < 350) continue;

			TVector3 p1_mom = p1_par.Mom3();
			TVector3 p2_mom = p2_par.Mom3();

			TVector3 p3_mom = p1_mom + p2_mom;

			TVector3 pv(v2x-vx,v2y-vy,v2z-vz);
			double ang = TMath::ACos( p3_mom*pv/(p3_mom.Mag() * pv.Mag()) ) ;
			if(ang >0.01) continue;
			if(p1_mom.Pt(p3_mom) <0.023) continue;

			if ( TMath::Abs( pv.Mag() ) / TMath::Sqrt( v.Cov(5) + v2.Cov(5) + v.Cov(2) + v2.Cov(2) + v.Cov(1) + v2.Cov(1) ) < 2.0) continue;


			pt1 = p1_mom.Pt(p3_mom);
			pt2 = p2_mom.Pt(p3_mom);
			alpha = (TMath::Sqrt(p1_mom.Mag2() - pt1*pt1) - TMath::Sqrt(p2_mom.Mag2() - pt2*pt2))/ (TMath::Sqrt(p1_mom.Mag2() - pt1*pt1) + TMath::Sqrt(p2_mom.Mag2() - pt2*pt2));

			lv_pip = p1.ParInVtx(iv2).LzVec(m_pi);
			lv_pim = p2.ParInVtx(iv2).LzVec(m_pi);
			lv_k0 = lv_pip + lv_pim;

			PaTPar par1, par2;
			if( !tr_p1.Extrapolate(z_rich,par1) ) continue;
			if( !tr_p2.Extrapolate(z_rich,par2) ) continue;

			for(int i = 0; i<6; i++) pp_lh[i] = pid.GetLike(i, tr_p1);
			for(int i = 0; i<6; i++) pm_lh[i] = pid.GetLike(i, tr_p2);
			pm_lh[6] = pid.SecondLike(tr_p1);
			pp_lh[6] = pid.SecondLike(tr_p2);
			pm_ch = pid.Theta_Ch(tr_p1);
			pp_ch = pid.Theta_Ch(tr_p2);

// 			if( par1.X()*par1.X() + par1.Y()*par1.Y() < 25.  ) continue;
// 			if( par2.X()*par2.X() + par2.Y()*par2.Y() < 25.  ) continue;

			pp_x     = par1.X();
			pp_y     = par1.Y();
			pp_mom   = par1.Mom();
			pp_theta = par1.Theta();
			pm_x     = par2.X();
			pm_y     = par2.Y();
			pm_mom   = par2.Mom();
			pm_theta = par2.Theta();

			if(alpha > 0){
				lv_p  = p1.ParInVtx(iv2).LzVec(m_p);
				lv_pi = p2.ParInVtx(iv2).LzVec(m_pi);
				lv_lambda = lv_p + lv_pi;
			} else if (alpha < 0){
				lv_p  = p2.ParInVtx(iv2).LzVec(m_p);
				lv_pi = p1.ParInVtx(iv2).LzVec(m_pi);
				lv_lambda = lv_p + lv_pi;
			}else continue;


			if(fabs(lv_lambda.Mag() - m_lambda) >= 0.15  && fabs(lv_k0.Mag() - m_k0) >= 0.15 )   continue;

			tree->Fill();
		}




		/* old exclusive phi selection */
		/*
		if(v.NOutParticles() !=3 ) continue;
		int i1=-1, i2=-1;

		for(int l = 0; l< v.NOutParticles();l++){
			int i_part = v.iOutParticle(l);
			if(i_part == imu1) continue;
			if(i1 ==-1) i1 = i_part;
			else i2 = i_part;
		}
		PaParticle p1 = e.vParticle(i1);
		PaParticle p2 = e.vParticle(i2);

		if(p1.Q() == p2.Q()) continue;
		if(p1.Q() < p2.Q()) swap(p1,p2);   // now first particle is positive

		if( p1.ParInVtx(Bvtx).Mom()<2. || p1.ParInVtx(Bvtx).Mom()>70.) continue;
		if( p2.ParInVtx(Bvtx).Mom()<2. || p2.ParInVtx(Bvtx).Mom()>70.) continue;

		lv_kp = p1.ParInVtx(Bvtx).LzVec(m_k);
		lv_km = p2.ParInVtx(Bvtx).LzVec(m_k);
		lv_phi = lv_kp + lv_km;

		TLorentzVector lv_ppp;
		lv_ppp.SetXYZM(0.,0.,0.,m_p);
		double MX2 = (lv_beam-lv_scat+lv_p-lv_phi).Mag2();
		Emiss = (MX2 - m_p*m_p)/(2*m_p);
		if( TMath::Abs(Emiss) > 2.5 ) continue;

		if(fabs(lv_phi.Mag() - m_phi) >= 0.15 )  continue;

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
		pm_ch = pid.Theta_Ch(tr_p1);
		pp_ch = pid.Theta_Ch(tr_p2);

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
		*/
	}

}

void UserJobEnd13010(){
	const Phast& ph = Phast::Ref();

	if(ph.print) cout<<"[ UserJobEnd0 has been called ]"<<endl;
	ph.h_file->cd();

}
