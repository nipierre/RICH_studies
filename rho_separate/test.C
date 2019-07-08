void test()
{
//=========Macro generated from canvas: c/
//=========  (Thu Apr 25 15:54:21 2019) by ROOT version 6.14/06
   TCanvas *c = new TCanvas("c", "",0,0,800,600);
   c->SetHighLightColor(2);
   c->Range(0,0,1,1);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   
   TH2D *test_hist__1 = new TH2D("test_hist__1","",500,-1,1,500,0,0.3);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   test_hist__1->SetLineColor(ci);
   test_hist__1->GetXaxis()->SetLabelFont(42);
   test_hist__1->GetXaxis()->SetLabelSize(0.035);
   test_hist__1->GetXaxis()->SetTitleSize(0.035);
   test_hist__1->GetXaxis()->SetTitleFont(42);
   test_hist__1->GetYaxis()->SetLabelFont(42);
   test_hist__1->GetYaxis()->SetLabelSize(0.035);
   test_hist__1->GetYaxis()->SetTitleSize(0.035);
   test_hist__1->GetYaxis()->SetTitleOffset(0);
   test_hist__1->GetYaxis()->SetTitleFont(42);
   test_hist__1->GetZaxis()->SetLabelFont(42);
   test_hist__1->GetZaxis()->SetLabelSize(0.035);
   test_hist__1->GetZaxis()->SetTitleSize(0.035);
   test_hist__1->GetZaxis()->SetTitleFont(42);
   test_hist__1->Draw("col");
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
