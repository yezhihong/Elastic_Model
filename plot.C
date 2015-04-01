void plot(){

   const double Norm_E = 1.90e-7;
   const double Norm_Q = 2.09e-4;
   const double Norm_M = 1.90e-7;
   
   TCanvas *c1 = new TCanvas("c1","c1",800,600);
  
   c1->Divide(1,4);

   ////////////////////////////////////////////////////////////////////////////////
   c1->cd(1); 
   TFile *f1 = new TFile("XS_Deut_Elastic_1100.root","r");
   TH1D* h_theta = new TH1D("h_theta","Deut at E0=1100MeV",1000, 0.4,8.4);
   h_theta->SetXTitle("#theta (Degree)");
   h_theta->SetYTitle("Counts");

   T->SetLineColor(2);
   T->Draw("Theta>>h_theta",Form("XS*%e*(Ep>1100*0.9)", Norm_E));

   TFile *f2 = new TFile("XS_Deut_QE_1100.root","r");
   T->SetLineColor(4);
   T->Draw("Theta",Form("(XS_QE+XS_DIS)*%e*(Ep>1100*0.9)", Norm_Q),"same");

   TFile *f3 = new TFile("XS_Moller_1100_New.root","r");
   T->SetLineColor(6);
   T->Draw("Theta",Form("XS_M*%e*(Ep>1100*0.9)", Norm_M),"same");

   ////////////////////////////////////////////////////////////////////////////////
   c1->cd(2); 
   TFile *f11 = new TFile("XS_Deut_Elastic_2200.root","r");
   TH1D* h_theta2 = new TH1D("h_theta2","Deut at E0=2200MeV",1000, 0.0,8.40);
   h_theta2->SetXTitle("#theta (Degree)");
   h_theta2->SetYTitle("Counts");

   T->SetLineColor(2);
   T->Draw("Theta>>h_theta2",Form("XS*%e*(Ep>2200*0.9)", Norm_E));

   TFile *f21 = new TFile("XS_Deut_QE_2200.root","r");
   T->SetLineColor(4);
   T->Draw("Theta",Form("(XS_QE+XS_DIS)*%e*(Ep>2200*0.9)", Norm_Q*2200/1100),"same");

   TFile *f31 = new TFile("XS_Moller_2200_New.root","r");
   T->SetLineColor(6);
   T->Draw("Theta",Form("XS_M*%e*(Ep>2200*0.9)", Norm_M),"same");


   ////////////////////////////////////////////////////////////////////////////////
   c1->cd(3); 

   TFile *f12 = new TFile("XS_Deut_Elastic_3300.root","r");
   TH1D* h_theta3 = new TH1D("h_theta3","Deut at E0=3300MeV",1000, 0.00,8.40);
   h_theta3->SetXTitle("#theta (Degree)");
   h_theta3->SetYTitle("Counts");

   T->SetLineColor(2);
   T->Draw("Theta>>h_theta3",Form("XS*%e*(Ep>3300*0.9)", Norm_E));

   TFile *f22 = new TFile("XS_Deut_QE_3300.root","r");
   T->SetLineColor(4);
   T->Draw("Theta",Form("(XS_QE+XS_DIS)*%e*(Ep>3300*0.9)", Norm_Q*3300/1100),"same");

   TFile *f32 = new TFile("XS_Moller_3300_New.root","r");
   T->SetLineColor(6);
   T->Draw("Theta",Form("XS_M*%e*(Ep>3300*0.9)", Norm_M),"same");


   ////////////////////////////////////////////////////////////////////////////////
   c1->cd(4); 

   TFile *f13 = new TFile("XS_Deut_Elastic_4400.root","r");
   TH1D* h_theta4 = new TH1D("h_theta4","Deut at E0=4400MeV",1000, 0.00,8.40);
   h_theta4->SetXTitle("#theta (Degree)");
   h_theta4->SetYTitle("Counts");

   T->SetLineColor(2);
   T->Draw("Theta>>h_theta4",Form("XS*%e*(Ep>4400*0.9)", Norm_E));

   TFile *f23 = new TFile("XS_Deut_QE_4400.root","r");
   T->SetLineColor(4);
   T->Draw("Theta",Form("(XS_QE+XS_DIS)*%e*(Ep>4400*0.9)", Norm_Q*4400/1100),"same");

   TFile *f33 = new TFile("XS_Moller_4400_New.root","r");
   T->SetLineColor(6);
   T->Draw("Theta",Form("XS_M*%e*(Ep>4400*0.9)", Norm_M),"same");



}
