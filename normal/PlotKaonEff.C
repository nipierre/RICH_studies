void PlotKaonEff()
{
  TGraphErrors* Eff[2][3];
  Double_t EfficienciesKp1[10];
  Double_t EfficienciesKp2[10];
  Double_t EfficienciesKp3[10];
  Double_t EfficienciesKp1e[10];
  Double_t EfficienciesKp2e[10];
  Double_t EfficienciesKp3e[10];
  Double_t EfficienciesKm1[10];
  Double_t EfficienciesKm2[10];
  Double_t EfficienciesKm3[10];
  Double_t EfficienciesKm1e[10];
  Double_t EfficienciesKm2e[10];
  Double_t EfficienciesKm3e[10];
  Double_t ph[10] = {12.5,14,16,18,20.5,23.5,26,28.5,32.5,37.5};
  TLine l1(10,0,40,0); l1.SetLineStyle(2);
  TLine l2(10,1,40,1); l2.SetLineStyle(2);

  ifstream kp("KplusEff.txt");

  for(int i=0; i<10; i++)
  {
    kp >> EfficienciesKp1[i];
    kp >> EfficienciesKp1e[i];
    kp >> EfficienciesKp2[i];
    kp >> EfficienciesKp2e[i];
    kp >> EfficienciesKp3[i];
    kp >> EfficienciesKp3e[i];
  }

  kp.close();

  ifstream km("KminusEff.txt");

  for(int i=0; i<10; i++)
  {
    km >> EfficienciesKm1[i];
    km >> EfficienciesKm1e[i];
    km >> EfficienciesKm2[i];
    km >> EfficienciesKm2e[i];
    km >> EfficienciesKm3[i];
    km >> EfficienciesKm3e[i];
  }

  km.close();

  TCanvas* c1 = new TCanvas("Efficiencies K^+","Efficiencies K^+");
  TCanvas* c2 = new TCanvas("Efficiencies K^-","Efficiencies K^-");


  Eff[1][0] = new TGraphErrors(10,ph,EfficienciesKp1,0,EfficienciesKp1e);
  Eff[1][1] = new TGraphErrors(10,ph,EfficienciesKp2,0,EfficienciesKp2e);
  Eff[1][2] = new TGraphErrors(10,ph,EfficienciesKp3,0,EfficienciesKp3e);
  Eff[0][0] = new TGraphErrors(10,ph,EfficienciesKm1,0,EfficienciesKm1e);
  Eff[0][1] = new TGraphErrors(10,ph,EfficienciesKm2,0,EfficienciesKm2e);
  Eff[0][2] = new TGraphErrors(10,ph,EfficienciesKm3,0,EfficienciesKm3e);

  TLegend leg1(.55,.38,.85,.62);
  leg1.AddEntry(Eff[1][0], "#font[12]{P(#pi^{+} #rightarrow K^{+})}");
  leg1.AddEntry(Eff[1][1], "#font[12]{P(K^{+} #rightarrow K^{+})}");
  leg1.AddEntry(Eff[1][2], "#font[12]{P(p #rightarrow K^{+})}");

  TLegend leg2(.55,.38,.85,.62);
  leg2.AddEntry(Eff[0][0], "#font[12]{P(#pi^{-} #rightarrow K^{-})}");
  leg2.AddEntry(Eff[0][1], "#font[12]{P(K^{-} #rightarrow K^{-})}");
  leg2.AddEntry(Eff[0][2], "#font[12]{P(#bar{p} #rightarrow K^{-})}");

  c1->cd();
  Eff[1][0]->SetMarkerStyle(20);
  Eff[1][0]->SetMarkerSize(1.5);
  Eff[1][0]->SetMarkerColor(4);
  Eff[1][0]->SetMinimum(-0.05);
  Eff[1][0]->SetMaximum(1.05);
  Eff[1][0]->SetTitle("");
  Eff[1][0]->GetXaxis()->SetTitle("#font[12]{p_{h}} #font[12]{(}GeV/#font[12]{c)}");
  Eff[1][0]->GetXaxis()->SetTitleSize(0.05);
  Eff[1][0]->GetXaxis()->SetTitleOffset(0.8);
  Eff[1][0]->GetYaxis()->SetTitle("#font[12]{P(h^{+} #rightarrow K^{+})}");
  Eff[1][0]->GetYaxis()->SetTitleSize(0.05);
  Eff[1][0]->GetYaxis()->SetTitleOffset(0.8);
  Eff[1][0]->Draw("AP");
  Eff[1][1]->SetMarkerStyle(33);
  Eff[1][1]->SetMarkerSize(2);
  Eff[1][1]->SetMarkerColor(8);
  Eff[1][1]->Draw("PSAME");
  Eff[1][2]->SetMarkerStyle(21);
  Eff[1][2]->SetMarkerSize(1.5);
  Eff[1][2]->SetMarkerColor(2);
  Eff[1][2]->Draw("PSAME");
  c1->Update();
  l1.Draw("SAME"); l2.Draw("SAME"); leg1.Draw("SAME");
  c1->Update();

  c2->cd();
  Eff[0][0]->SetMarkerStyle(20);
  Eff[0][0]->SetMarkerSize(1.5);
  Eff[0][0]->SetMarkerColor(4);
  Eff[0][0]->SetMinimum(-0.05);
  Eff[0][0]->SetMaximum(1.05);
  Eff[0][0]->SetTitle("");
  Eff[0][0]->GetXaxis()->SetTitle("#font[12]{p_{h}} #font[12]{(}GeV/#font[12]{c)}");
  Eff[0][0]->GetXaxis()->SetTitleSize(0.05);
  Eff[0][0]->GetXaxis()->SetTitleOffset(0.8);
  Eff[0][0]->GetYaxis()->SetTitle("#font[12]{P(h^{-} #rightarrow K^{-})}");
  Eff[0][0]->GetYaxis()->SetTitleSize(0.05);
  Eff[0][0]->GetYaxis()->SetTitleOffset(0.8);
  Eff[0][0]->Draw("AP");
  Eff[0][1]->SetMarkerStyle(33);
  Eff[0][1]->SetMarkerSize(2);
  Eff[0][1]->SetMarkerColor(8);
  Eff[0][1]->Draw("PSAME");
  Eff[0][2]->SetMarkerStyle(21);
  Eff[0][2]->SetMarkerSize(1.5);
  Eff[0][2]->SetMarkerColor(2);
  Eff[0][2]->Draw("PSAME");
  c2->Update();
  l1.Draw("SAME"); l2.Draw("SAME"); leg2.Draw("SAME");
  c2->Update();

  c1->Print("KplusEff.png");
  c2->Print("KminusEff.png");

}
