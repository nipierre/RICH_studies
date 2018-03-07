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

//****************************/
/* RICH Matrix				*/
/* inclusive rho0 			*/
/* selection				*/
/****************************/



// Only for 2006 data
void RescaleMom7102(PaEvent & e, bool faster_mode)
{
  static bool if_first=1;

  static float x_sm2=0;
  static float y_sm2=0;
  static float z_sm2=0;
  static float z_sm1=0;
  static float z_tar=0;
  static float detdatcorr=1.0000;
  static float sm1corr=1.0000;
  static float sm2corr=1.0000;
  static float sm1sm2corr=1.0000;


  if(if_first)
  {
    const PaSetup& setup = PaSetup::Ref();
    PaMagInfo* m = PaSetup::Ref().PtrMagField()->getMagInfo();

    const int MagnetsNumber=PaSetup::Ref().PtrMagField()->getNumOfMags();

    x_sm2 = m[MagnetsNumber-1].xcm/10.; // SM2
    y_sm2 = m[MagnetsNumber-1].ycm/10.; // SM2
    z_sm2 = m[MagnetsNumber-1].zcm/10.; // SM2
    z_sm1 = m[MagnetsNumber-2].zcm/10.; // SM1
    z_tar = setup.TargetCenterZ();      // target

    float xx=x_sm2-67.5;
    float yy=y_sm2-45.0;
    float zz=z_sm2+1.5;
    float bx,by,bz;
    setup.MagField(xx,yy,zz,bx,by,bz);


    int usedfield=0;
    if(fabs(by)>1.75) usedfield=5000;
    if(fabs(by)>1.45&&fabs(by)<1.75) usedfield=4000;





    if(usedfield==5000&&fabs(fabs(by/1.91221)-1)>0.003) return; // corr done in det.dat?
    if(usedfield==4000&&fabs(fabs(by/1.57905)-1)>0.003) return; //corr done in det.dat?

    switch(usedfield)
    {
      case 4000:
        detdatcorr=fabs(by/1.57905);
        sm2corr=1.0088/detdatcorr;
        sm1sm2corr=1.0088/detdatcorr;
      break;
      case 5000:
        detdatcorr=fabs(by/1.91221);
        sm2corr=1.0040/detdatcorr;
        sm1sm2corr=1.0040/detdatcorr;
      break;
    }
    if(faster_mode) if_first=0;
  }

    vector<PaTrack>& TrVec=  e.vTrack();
    vector<PaParticle>& PaVec=  e.vParticle();

    for(unsigned int i=0; i<TrVec.size();i++)
    {

      int TrackType=-1;
      if(TrVec[i].ZLast()<z_tar) TrackType=0;  //beam
      if(TrackType!=0 && TrVec[i].ZFirst()<z_sm1 && TrVec[i].ZLast()<z_sm2) TrackType=1;  //sm1
      if(TrackType!=0 && TrVec[i].ZFirst()>z_sm1 && TrVec[i].ZLast()>z_sm2) TrackType=2;  //sm2
      if(TrackType!=0 && TrVec[i].ZFirst()<z_sm1 && TrVec[i].ZLast()>z_sm2) TrackType=3;  //sm1+sm2


      if(!TrackType) continue; // do nothing with beam track

      vector<PaTPar>& PaTTrVec=TrVec[i].vTPar();
      for(unsigned int j=0; j<PaTTrVec.size();j++)
      {
        switch(TrackType)
        {
          case 1:
            PaTTrVec[j](5)/=sm1corr;
            break;
          case 2:
            PaTTrVec[j](5)/=sm2corr;
            break;
          case 3:
            PaTTrVec[j](5)/=sm1sm2corr;
            break;
          // I don't touch covariance matrix.
        }
      }

      if(TrVec[i].iParticle()<0) continue; //some tracks don't have corresponding particle
      vector<PaTPar>& PaTPaVec=PaVec[TrVec[i].iParticle()].vFitPar();
      for(unsigned int j=0; j<PaTPaVec.size();j++)
      {
        switch(TrackType)
        {
          case 1:
            PaTPaVec[j](5)/=sm1corr;
            break;
          case 2:
            PaTPaVec[j](5)/=sm2corr;
            break;
          case 3:
            PaTPaVec[j](5)/=sm1sm2corr;
            break;
        }
      }
    }

}

// Only for 2006 data
void RemoveBadMiddleTrigger7102(PaEvent& e, const PaTrack& mup_tr) {
	// Extract information about the Hodoscope planes once and store it in std::map for fast access
	// If run number changes refresh the information
	static bool first = 1;
	static int run_number = 0;
	static double minX = -9999;
	static double Z = 0;
	if(first || run_number != e.RunNum()) {
		run_number = e.RunNum();
		int idet;
		idet = PaSetup::Ref().iDetector("HM05X1_d");
		if(idet>0) {
			const PaDetect& HM05 =  PaSetup::Ref().Detector(idet);
			minX = HM05.Uorig() + 2.5*HM05.Pitch(); // Remove 3 strips
			Z = HM05.Z();
		} else {
			cerr<<"RemoveBadMiddleTrigger ERROR: Detector HM05Y1_d not found!"<<endl;
			exit(1);
		}
		first = 0;
	}

	if((e.TrigMask()&0x102) == 0) return;

	int Npars = mup_tr.NTPar();
	if(Npars == 0) {
		cerr<<"RemoveBadMiddleTrigger ERROR: track does not have associated helixes!"<<endl;
		return;
	}

	PaTPar partr = mup_tr.vTPar(Npars-1); //Track parameters in last measured point
	PaTPar result;
	bool res = partr.Extrapolate(Z, result, false);
	if(result(1) < minX || !res) {
		int new_mask = e.TrigMask() & (~0x102);
		// cout<<"Old: "<<hex<<e.TrigMask()<<" New:"<<new_mask<<dec<<endl;
		e.vHeader()[1] = new_mask;
	}
}

void UserEvent7102(PaEvent& e){
	// RescaleMom7102(e,true);
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

	static TTree* tree(NULL);
	static int Run, TriggerMask;
	static Long64_t Evt;
	static float vx, vy ,vz;
	static double Q2, xbj, y;
	static double Emiss;
	static TLorentzVector lv_beam, lv_scat, lv_pip, lv_pim, lv_rho;

	static double pp_x, pp_y, pp_mom, pp_theta, pm_x, pm_y, pm_mom, pm_theta;
	static int Nout;
	static double pp_lh[7],pm_lh[7];
	static double pp_pt, pm_pt, alpha;
	static double pi_thr, k_thr, p_thr;

	const float m_pi = G3partMass[8];
	const float m_k = G3partMass[11];
	const float m_mu = G3partMass[5];
	const float m_p = G3partMass[14];
	const float m_rho = G3partMass[63];

  TargetCell* fTcell;
	fTcell = new TargetCell();

	static bool first(true);
	if(first){
		tree = new TTree("tree","selected events");
		tree->Branch("Run",    			&Run,    			"Run/I");
		tree->Branch("Evt",    			&Evt,    			"Evt/L");
		tree->Branch("TriggerMask", 	&TriggerMask, 		"TriggerMask/I");

		tree->Branch("lv_beam",			"TLorentzVector", 	&lv_beam);
		tree->Branch("lv_scat",			"TLorentzVector", 	&lv_scat);

		tree->Branch("Q2",				&Q2,				"Q2/D");
		tree->Branch("xbj",				&xbj,				"xbj/D");
		tree->Branch("y",				&y,					"y/D");
		tree->Branch("Emiss",			&Emiss,				"Emiss/D");
		tree->Branch("vx",				&vx,				"vx/F");
		tree->Branch("vy",				&vy,				"vy/F");
		tree->Branch("vz",				&vz,				"vz/F");
		tree->Branch("Nout", 			&Nout,		 		"Nout/I");
		tree->Branch("k_thr",			&k_thr,				"k_thr/D");
		tree->Branch("pi_thr",			&pi_thr,			"pi_thr/D");
		tree->Branch("p_thr",			&p_thr,				"p_thr/D");

		tree->Branch("lv_pip",			"TLorentzVector", 	&lv_pip);
		tree->Branch("lv_pim",			"TLorentzVector", 	&lv_pim);
		tree->Branch("lv_rho",			"TLorentzVector", 	&lv_rho);

		tree->Branch("pp_pt",			&pp_pt,				"pp_pt/D");
		tree->Branch("pm_pt",			&pm_pt,				"pm_pt/D");
		tree->Branch("alpha",			&alpha,				"alpha/D");

		tree->Branch("pp_x",			&pp_x,				"pp_x/D");
		tree->Branch("pp_y",			&pp_y,				"pp_y/D");
		tree->Branch("pp_mom",			&pp_mom,			"pp_mom/D");
		tree->Branch("pp_theta",		&pp_theta,			"pp_theta/D");
		tree->Branch("pp_lh",			pp_lh,				"pp_lh[7]/D");

		tree->Branch("pm_x",			&pm_x,				"pm_x/D");
		tree->Branch("pm_y",			&pm_y,				"pm_y/D");
		tree->Branch("pm_mom",			&pm_mom,			"pm_mom/D");
		tree->Branch("pm_theta",		&pm_theta,			"pm_theta/D");
		tree->Branch("pm_lh",			pm_lh,				"pm_lh[7]/D");

		first=false;
	}

	double n_rich = 1 + PaMetaDB::Ref().NminusOne(Run)/1e6;
	if (PaMetaDB::Ref().NminusOne(Run) < 0) n_rich = 1 +  PaSetup::Ref().NminusOne()/1e6;
	k_thr = (m_k/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));
	pi_thr = (m_pi/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));
	p_thr = (m_p/n_rich) * 1/TMath::Sqrt(1-1/(n_rich*n_rich));

	Run        =  e.RunNum();
	Evt        =  e.UniqueEvNum();

	const PaSetup& setup = PaSetup::Ref();
	const double z_rich = setup.Rich().DetPos(0).Z();// - 1.65;

	PaPid pid;

	int Bvtx = e.iBestPrimaryVertex();
	for(int iv = 0; iv < e.NVertex(); iv++){
		const PaVertex& v = e.vVertex(iv);
		if(! v.IsPrimary()) continue;
		if(Bvtx !=-1 && Bvtx != iv) continue;

		int imu =v.iMuPrim(false,true,true,true);
		int ibeam = v.InParticle();

		if(ibeam ==-1) continue;
		if(imu == -1) continue;

		const PaParticle & pa_scat = e.vParticle(imu);
		const PaParticle & pa_beam = e.vParticle(ibeam);

		const PaTPar & par_scat = pa_scat.ParInVtx(iv);
		const PaTPar & par_beam = pa_beam.ParInVtx(iv);

		const PaTrack & tr_scat = e.vTrack(pa_scat.iTrack());
		// const PaTrack & tr_beam = e.vTrack(pa_beam.iTrack());

		// RemoveBadMiddleTrigger7102(e, tr_scat);
		TriggerMask = (e.TrigMask() & 0xffff);

    //if( !fTcell->TargetCell::InTarget(par_beam,Run) ) continue;  THIS FUNCTION IS FUCKED UP. TODO SEE WHY
    if( !fTcell->TargetCell::CrossCells(par_beam, Run)) continue;

		vx = v.X();
		vy = v.Y();
		vz = v.Z();

		lv_scat = par_scat.LzVec(m_mu);
		lv_beam = par_beam.LzVec(m_mu);

		Q2 = PaAlgo::Q2(lv_beam,lv_scat);
		xbj = PaAlgo::xbj(lv_beam,lv_scat);
		y = (lv_beam.E()-lv_scat.E())/lv_beam.E();

		if(y<0.1 || y>0.9) continue;

		Nout = v.NOutParticles();

		if(Nout < 3) continue; // at least mu' + pi+ + pi-
		// BackPropLH = pa_beam.BackPropLH();
		// if(BackPropLH<1 && BackPropLH>0.005) cut_backprob = 1;
		// cov14 = par_beam.Cov14();

		for(int l = 0; l< v.NOutParticles();l++){
			int i_part1 = v.iOutParticle(l);
			if(i_part1 == imu) continue;
			for(int m = l; m< v.NOutParticles();m++){
				int i_part2 = v.iOutParticle(m);
				if(i_part2 == imu) continue;
				if(i_part2 == i_part1) continue;

				PaParticle p1 = e.vParticle(i_part1);
				PaParticle p2 = e.vParticle(i_part2);
				if(p1.Q() == p2.Q()) continue;
				if(p1.Q() < p2.Q()) swap(p1,p2);   // now first particle is positive

				const PaTPar & par1 = p1.ParInVtx(iv);
				const PaTPar & par2 = p2.ParInVtx(iv);
				const PaTrack & tr1 = e.vTrack(p1.iTrack());
				const PaTrack & tr2 = e.vTrack(p2.iTrack());

				lv_pip = par1.LzVec(m_pi);
				TLorentzVector lv_kp = par1.LzVec(m_k);
				int zf1 = tr1.ZFirst();
				int zl1 = tr1.ZLast();
				int xx0_1 = tr1.XX0();
				float pip_chi = tr1.Chi2tot()/(float)tr1.Ndf();
				// if(tr_part1.CanBeMuon() ) can_be = 1;
				if(pip_chi > 10) continue;
				if(xx0_1 > 10) continue;
				if(zl1 <350) continue;
				if(zf1 >350) continue;
				if(zl1 >3300) continue;
				if(!tr1.HasMom()) continue;
				if(lv_pip.Vect().Mag()<1.) continue;

				lv_pim = par2.LzVec(m_pi);
				TLorentzVector lv_km = par2.LzVec(m_k);
				int zf2 = tr2.ZFirst();
				int zl2 = tr2.ZLast();
				int xx0_2 = tr2.XX0();
				float pim_chi = tr2.Chi2tot()/(float)tr2.Ndf();
				// if(tr_part2.CanBeMuon() ) can_be = 1;
				if(pim_chi > 10) continue;
				if(xx0_2 > 10) continue;
				if(zl2 <350) continue;
				if(zf2 >350) continue;
				if(zl2 >3300) continue;
				if(!tr2.HasMom()) continue;
				if(lv_pim.Vect().Mag()<1.) continue;

				lv_rho = lv_pip+lv_pim;

				if( TMath::Abs(lv_rho.Mag() - m_rho) >=0.25) continue;
				if( (lv_kp+lv_km).Mag() <1.04) continue;

				for(int i = 0; i<6; i++) pp_lh[i] = pid.GetLike(i, tr1);
				for(int i = 0; i<6; i++) pm_lh[i] = pid.GetLike(i, tr2);
				pm_lh[6] = pid.SecondLike(tr1);
				pp_lh[6] = pid.SecondLike(tr2);

				PaTPar tpar1, tpar2;
				if( !tr1.Extrapolate(z_rich,tpar1) ) continue;
				if( !tr2.Extrapolate(z_rich,tpar2) ) continue;;

				TVector3 p1_mom = par1.Mom3();
				TVector3 p2_mom = par2.Mom3();

				TVector3 p3_mom = p1_mom + p2_mom;
				if(p1_mom.Pt(p3_mom) <0.023) continue;
				if(p2_mom.Pt(p3_mom) <0.023) continue;

				pp_pt = p1_mom.Pt(p3_mom);
				pm_pt = p2_mom.Pt(p3_mom);
				alpha = (TMath::Sqrt(p1_mom.Mag2() - pp_pt*pp_pt) - TMath::Sqrt(p2_mom.Mag2() - pm_pt*pm_pt))/ (TMath::Sqrt(p1_mom.Mag2() - pp_pt*pp_pt) + TMath::Sqrt(p2_mom.Mag2() - pm_pt*pm_pt));

				pp_x     = tpar1.X();
				pp_y     = tpar1.Y();
				pp_mom   = tpar1.Mom();
				pp_theta = tpar1.Theta();
				pm_x     = tpar2.X();
				pm_y     = tpar2.Y();
				pm_mom   = tpar2.Mom();
				pm_theta = tpar2.Theta();

				TLorentzVector lv_ppp;
				lv_ppp.SetXYZM(0.,0.,0.,m_p);
				double MX2 = (lv_beam-lv_scat+lv_ppp-lv_rho).Mag2();
				Emiss = (MX2 - m_p*m_p)/(2*m_p);
				tree->Fill();
			}
		}
	}
}
