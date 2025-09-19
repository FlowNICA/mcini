#include <fstream>

R__LOAD_LIBRARY(libMcIniData.so)

using namespace std;

// Arguments: nev = -1 (all events), N>0 - specific event number; inFileName/outFileName - input/output; suffix = sym, asym, auau;
void getEcc(long nev, std::string inFileList, std::string outFileName, std::string suffix = "sym", double t_cut = -1.)
{
  const int nbins_b = 20;
  const int nbins_mult = 500;

  TChain *tree = new TChain("events");
  std::ifstream infile(inFileList.c_str());
  std::string line;
  if (!infile.is_open()){
    std::cerr << "Cannot open the file: " << inFileList.c_str() << "! Exit." << std::endl;
    return;
  }
  std::cout << "getEcc: Input file list: " << inFileList.c_str() << std::endl;
  while (std::getline(infile,line)){
    //std::cout << "getEcc: adding file " << line.c_str() << "." << std::endl;
    tree->Add(line.c_str());
  }
  EventInitialState *iniState{nullptr};
  UEvent *event{nullptr};
  tree->SetBranchAddress("iniState", &iniState);
  tree->SetBranchAddress("event", &event);
  
  // Output
  TFile *fo = new TFile(outFileName.c_str(), "recreate");

  TProfile *p1ecc1vsb = new TProfile(Form("p1_%s_ecc1vsb", suffix.c_str()),"#varepsilon_{1} vs b;b, fm;#varepsilon_{1}", nbins_b, 0., 20.); 
  TProfile *p1ecc2vsb = new TProfile(Form("p1_%s_ecc2vsb", suffix.c_str()),"#varepsilon_{2} vs b;b, fm;#varepsilon_{2}", nbins_b, 0., 20.); 
  TProfile *p1ecc3vsb = new TProfile(Form("p1_%s_ecc3vsb", suffix.c_str()),"#varepsilon_{3} vs b;b, fm;#varepsilon_{3}", nbins_b, 0., 20.); 
  TProfile *p1ecc4vsb = new TProfile(Form("p1_%s_ecc4vsb", suffix.c_str()),"#varepsilon_{4} vs b;b, fm;#varepsilon_{4}", nbins_b, 0., 20.); 
  TProfile *p1ecc5vsb = new TProfile(Form("p1_%s_ecc5vsb", suffix.c_str()),"#varepsilon_{5} vs b;b, fm;#varepsilon_{5}", nbins_b, 0., 20.);
  TProfile *p1ecc1p2vsb = new TProfile(Form("p1_%s_ecc1p2vsb", suffix.c_str()),"#varepsilon_{1}^{2} vs b;b, fm;#varepsilon_{1}^{2}", nbins_b, 0., 20.); 
  TProfile *p1ecc2p2vsb = new TProfile(Form("p1_%s_ecc2p2vsb", suffix.c_str()),"#varepsilon_{2}^{2} vs b;b, fm;#varepsilon_{2}^{2}", nbins_b, 0., 20.); 
  TProfile *p1ecc3p2vsb = new TProfile(Form("p1_%s_ecc3p2vsb", suffix.c_str()),"#varepsilon_{3}^{2} vs b;b, fm;#varepsilon_{3}^{2}", nbins_b, 0., 20.); 
  TProfile *p1ecc4p2vsb = new TProfile(Form("p1_%s_ecc4p2vsb", suffix.c_str()),"#varepsilon_{4}^{2} vs b;b, fm;#varepsilon_{4}^{2}", nbins_b, 0., 20.); 
  TProfile *p1ecc5p2vsb = new TProfile(Form("p1_%s_ecc5p2vsb", suffix.c_str()),"#varepsilon_{5}^{2} vs b;b, fm;#varepsilon_{5}^{2}", nbins_b, 0., 20.);
  TProfile *p1ecc1p4vsb = new TProfile(Form("p1_%s_ecc1p4vsb", suffix.c_str()),"#varepsilon_{1}^{4} vs b;b, fm;#varepsilon_{1}^{4}", nbins_b, 0., 20.); 
  TProfile *p1ecc2p4vsb = new TProfile(Form("p1_%s_ecc2p4vsb", suffix.c_str()),"#varepsilon_{2}^{4} vs b;b, fm;#varepsilon_{2}^{4}", nbins_b, 0., 20.); 
  TProfile *p1ecc3p4vsb = new TProfile(Form("p1_%s_ecc3p4vsb", suffix.c_str()),"#varepsilon_{3}^{4} vs b;b, fm;#varepsilon_{3}^{4}", nbins_b, 0., 20.); 
  TProfile *p1ecc4p4vsb = new TProfile(Form("p1_%s_ecc4p4vsb", suffix.c_str()),"#varepsilon_{4}^{4} vs b;b, fm;#varepsilon_{4}^{4}", nbins_b, 0., 20.); 
  TProfile *p1ecc5p4vsb = new TProfile(Form("p1_%s_ecc5p4vsb", suffix.c_str()),"#varepsilon_{5}^{4} vs b;b, fm;#varepsilon_{5}^{4}", nbins_b, 0., 20.);

  TH1D     *h1mult       = new TH1D(Form("h1_%s_mult", suffix.c_str()),"Multiplicity;dN/dN_{ch};n_{ch}", nbins_mult, 0., (double)nbins_mult);
  TProfile *p1ecc1vsmult = new TProfile(Form("p1_%s_ecc1vsmult", suffix.c_str()),"#varepsilon_{1} vs N_{ch};N_{ch};#varepsilon_{1}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc2vsmult = new TProfile(Form("p1_%s_ecc2vsmult", suffix.c_str()),"#varepsilon_{2} vs N_{ch};N_{ch};#varepsilon_{2}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc3vsmult = new TProfile(Form("p1_%s_ecc3vsmult", suffix.c_str()),"#varepsilon_{3} vs N_{ch};N_{ch};#varepsilon_{3}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc4vsmult = new TProfile(Form("p1_%s_ecc4vsmult", suffix.c_str()),"#varepsilon_{4} vs N_{ch};N_{ch};#varepsilon_{4}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc5vsmult = new TProfile(Form("p1_%s_ecc5vsmult", suffix.c_str()),"#varepsilon_{5} vs N_{ch};N_{ch};#varepsilon_{5}", nbins_mult, 0., (double)nbins_mult);
  TProfile *p1ecc1p2vsmult = new TProfile(Form("p1_%s_ecc1p2vsmult", suffix.c_str()),"#varepsilon_{1}^{2} vs N_{ch};N_{ch};#varepsilon_{1}^{2}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc2p2vsmult = new TProfile(Form("p1_%s_ecc2p2vsmult", suffix.c_str()),"#varepsilon_{2}^{2} vs N_{ch};N_{ch};#varepsilon_{2}^{2}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc3p2vsmult = new TProfile(Form("p1_%s_ecc3p2vsmult", suffix.c_str()),"#varepsilon_{3}^{2} vs N_{ch};N_{ch};#varepsilon_{3}^{2}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc4p2vsmult = new TProfile(Form("p1_%s_ecc4p2vsmult", suffix.c_str()),"#varepsilon_{4}^{2} vs N_{ch};N_{ch};#varepsilon_{4}^{2}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc5p2vsmult = new TProfile(Form("p1_%s_ecc5p2vsmult", suffix.c_str()),"#varepsilon_{5}^{2} vs N_{ch};N_{ch};#varepsilon_{5}^{2}", nbins_mult, 0., (double)nbins_mult);
  TProfile *p1ecc1p4vsmult = new TProfile(Form("p1_%s_ecc1p4vsmult", suffix.c_str()),"#varepsilon_{1}^{4} vs N_{ch};N_{ch};#varepsilon_{1}^{4}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc2p4vsmult = new TProfile(Form("p1_%s_ecc2p4vsmult", suffix.c_str()),"#varepsilon_{2}^{4} vs N_{ch};N_{ch};#varepsilon_{2}^{4}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc3p4vsmult = new TProfile(Form("p1_%s_ecc3p4vsmult", suffix.c_str()),"#varepsilon_{3}^{4} vs N_{ch};N_{ch};#varepsilon_{3}^{4}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc4p4vsmult = new TProfile(Form("p1_%s_ecc4p4vsmult", suffix.c_str()),"#varepsilon_{4}^{4} vs N_{ch};N_{ch};#varepsilon_{4}^{4}", nbins_mult, 0., (double)nbins_mult); 
  TProfile *p1ecc5p4vsmult = new TProfile(Form("p1_%s_ecc5p4vsmult", suffix.c_str()),"#varepsilon_{5}^{4} vs N_{ch};N_{ch};#varepsilon_{5}^{4}", nbins_mult, 0., (double)nbins_mult);

  long nEvents_tree = tree->GetEntries();
  long nEvents = (nev<0) ? nEvents_tree : std::min(nev, nEvents_tree);
  int Npart;
  double aversin1, avercos1, aver1;
  double aversin2, avercos2, aver2;
  double aversin3, avercos3, aver3;
  double aversin4, avercos4, aver4;
  double aversin5, avercos5, aver5;
  double averX, averY;
  double posX, posY;
  double ecc1, ecc2, ecc3, ecc4, ecc5;
  double ncount;
  double phi, rsqr, r;
  double pt, eta, p, charge;
  long refMult;

  for(int iEvent = 0; iEvent < nEvents; iEvent++) {
    if (iEvent%1000 == 0) std::cout << "Event [" << iEvent << "/" << nEvents << "]" << std::endl;
    tree->GetEntry(iEvent);
    Npart = 0;
    vector<Nucleon> nucleons = iniState->getNucleons();
    aversin1 = 0.; avercos1 = 0.; aver1 = 0.;
    aversin2 = 0.; avercos2 = 0.; aver2 = 0.;
    aversin3 = 0.; avercos3 = 0.; aver3 = 0.;
    aversin4 = 0.; avercos4 = 0.; aver4 = 0.;
    aversin5 = 0.; avercos5 = 0.; aver5 = 0.;
    averX = 0.; averY = 0.;
    ncount = 0;
    for(auto nucleon : nucleons)
    {
      TLorentzVector position = nucleon.getPosition();
      if (t_cut > 0 && position.T() >= t_cut) continue;
      if(nucleon.getCollisionType()==3 || nucleon.getCollisionType()==4)
      {
        ncount++;
        averX += position.X();
        averY += position.Y();
      }
    }
    if (ncount == 0) continue;
    averX /= ncount;
    averY /= ncount;

    ncount = 0;
    for(auto nucleon : nucleons)
    {
      TLorentzVector position = nucleon.getPosition();
      if (t_cut > 0 && position.T() >= t_cut) continue;
      if(nucleon.getCollisionType()==3 || nucleon.getCollisionType()==4)
      {
        posX = position.X() - averX;
        posY = position.Y() - averY;
        rsqr = posX*posX + posY*posY;
        r    = sqrt(rsqr);
        phi = TMath::ATan2(posY,posX);
        
        Npart++;

        aversin1 += TMath::Power(r,3.)*TMath::Sin(1.*phi);
        avercos1 += TMath::Power(r,3.)*TMath::Cos(1.*phi);
        aver1    += TMath::Power(r,3.);
        
        aversin2 += rsqr*TMath::Sin(2.*phi);
        avercos2 += rsqr*TMath::Cos(2.*phi);
        aver2    += rsqr;

        aversin3 += TMath::Power(r,3.)*TMath::Sin(3.*phi);
        avercos3 += TMath::Power(r,3.)*TMath::Cos(3.*phi);
        aver3    += TMath::Power(r,3.);

        aversin4 += TMath::Power(r,4.)*TMath::Sin(4.*phi);
        avercos4 += TMath::Power(r,4.)*TMath::Cos(4.*phi);
        aver4    += TMath::Power(r,4.);

        aversin5 += TMath::Power(r,5.)*TMath::Sin(5.*phi);
        avercos5 += TMath::Power(r,5.)*TMath::Cos(5.*phi);
        aver5    += TMath::Power(r,5.);

        ncount++;
      }
    }
    if (ncount == 0.) ncount = 1.;
    aversin1 /= ncount;
    avercos1 /= ncount;
    aver1    /= ncount;
    aversin2 /= ncount;
    avercos2 /= ncount;
    aver2    /= ncount;
    aversin3 /= ncount;
    avercos3 /= ncount;
    aver3    /= ncount;
    aversin4 /= ncount;
    avercos4 /= ncount;
    aver4    /= ncount;
    aversin5 /= ncount;
    avercos5 /= ncount;
    aver5    /= ncount;
    if (aver1 != 0) ecc1 = TMath::Sqrt(avercos1*avercos1 + aversin1*aversin1)/ aver1;
    if (aver2 != 0) ecc2 = TMath::Sqrt(avercos2*avercos2 + aversin2*aversin2)/ aver2;
    if (aver3 != 0) ecc3 = TMath::Sqrt(avercos3*avercos3 + aversin3*aversin3)/ aver3;
    if (aver4 != 0) ecc4 = TMath::Sqrt(avercos4*avercos4 + aversin4*aversin4)/ aver4;
    if (aver5 != 0) ecc5 = TMath::Sqrt(avercos5*avercos5 + aversin5*aversin5)/ aver5;
    if (iEvent%1000 == 0) std::cout << "\tNpart(std) = " << iniState->getNPart()
      << " Npart(custom) = " << Npart+1 << std::endl;

    p1ecc1vsb->Fill(event->GetB(), ecc1);
    p1ecc2vsb->Fill(event->GetB(), ecc2);
    p1ecc3vsb->Fill(event->GetB(), ecc3);
    p1ecc4vsb->Fill(event->GetB(), ecc4);
    p1ecc5vsb->Fill(event->GetB(), ecc5);
    p1ecc1p2vsb->Fill(event->GetB(), pow(ecc1, 2));
    p1ecc2p2vsb->Fill(event->GetB(), pow(ecc2, 2));
    p1ecc3p2vsb->Fill(event->GetB(), pow(ecc3, 2));
    p1ecc4p2vsb->Fill(event->GetB(), pow(ecc4, 2));
    p1ecc5p2vsb->Fill(event->GetB(), pow(ecc5, 2));
    p1ecc1p4vsb->Fill(event->GetB(), pow(ecc1, 4));
    p1ecc2p4vsb->Fill(event->GetB(), pow(ecc2, 4));
    p1ecc3p4vsb->Fill(event->GetB(), pow(ecc3, 4));
    p1ecc4p4vsb->Fill(event->GetB(), pow(ecc4, 4));
    p1ecc5p4vsb->Fill(event->GetB(), pow(ecc5, 4));

    refMult = 0;
    for (int iPart=0; iPart<event->GetNpa(); iPart++)
    {
      UParticle *particle = event->GetParticle(iPart);
      auto partPDG = (TParticlePDG*)TDatabasePDG::Instance()->GetParticle(particle->GetPdg());
      if (!partPDG) charge = 0.;
      else charge = 3.*partPDG->Charge();
      if (charge == 0) continue;
      pt = particle->GetMomentum().Pt();
      eta= particle->GetMomentum().Eta();
      if (pt < 0.15) continue;
      if (TMath::Abs(eta) > 0.5) continue;

      refMult++;
    }

    h1mult->Fill(refMult);
    p1ecc1vsmult->Fill(refMult, ecc1);
    p1ecc2vsmult->Fill(refMult, ecc2);
    p1ecc3vsmult->Fill(refMult, ecc3);
    p1ecc4vsmult->Fill(refMult, ecc4);
    p1ecc5vsmult->Fill(refMult, ecc5);
    p1ecc1p2vsmult->Fill(refMult, pow(ecc1, 2));
    p1ecc2p2vsmult->Fill(refMult, pow(ecc2, 2));
    p1ecc3p2vsmult->Fill(refMult, pow(ecc3, 2));
    p1ecc4p2vsmult->Fill(refMult, pow(ecc4, 2));
    p1ecc5p2vsmult->Fill(refMult, pow(ecc5, 2));
    p1ecc1p4vsmult->Fill(refMult, pow(ecc1, 4));
    p1ecc2p4vsmult->Fill(refMult, pow(ecc2, 4));
    p1ecc3p4vsmult->Fill(refMult, pow(ecc3, 4));
    p1ecc4p4vsmult->Fill(refMult, pow(ecc4, 4));
    p1ecc5p4vsmult->Fill(refMult, pow(ecc5, 4));

  }

  fo->cd();

  p1ecc1vsb->Write();
  p1ecc2vsb->Write();
  p1ecc3vsb->Write();
  p1ecc4vsb->Write();
  p1ecc5vsb->Write();
  p1ecc1p2vsb->Write();
  p1ecc2p2vsb->Write();
  p1ecc3p2vsb->Write();
  p1ecc4p2vsb->Write();
  p1ecc5p2vsb->Write();
  p1ecc1p4vsb->Write();
  p1ecc2p4vsb->Write();
  p1ecc3p4vsb->Write();
  p1ecc4p4vsb->Write();
  p1ecc5p4vsb->Write();

  h1mult->Write();
  p1ecc1vsmult->Write();
  p1ecc2vsmult->Write();
  p1ecc3vsmult->Write();
  p1ecc4vsmult->Write();
  p1ecc5vsmult->Write();
  p1ecc1p2vsmult->Write();
  p1ecc2p2vsmult->Write();
  p1ecc3p2vsmult->Write();
  p1ecc4p2vsmult->Write();
  p1ecc5p2vsmult->Write();
  p1ecc1p4vsmult->Write();
  p1ecc2p4vsmult->Write();
  p1ecc3p4vsmult->Write();
  p1ecc4p4vsmult->Write();
  p1ecc5p4vsmult->Write();
}
