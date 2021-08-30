#include "tglaubermc.hh"

using namespace std;


ClassImp(TGlauNucleon)
  //---------------------------------------------------------------------------------
void TGlauNucleon::RotateXYZ(Double_t phi, Double_t theta)
{
  TVector3 v(fX,fY,fZ);
  TVector3 vr;
  vr.SetMagThetaPhi(1,theta,phi);
  v.RotateUz(vr);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

void TGlauNucleon::RotateXYZ_3D(Double_t psiX, Double_t psiY, Double_t psiZ)
{
  TVector3 v(fX,fY,fZ);
  v.RotateX(psiX);
  v.RotateY(psiY);
  v.RotateZ(psiZ);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------
TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
  TNamed(iname,""),
  fN(iN),fR(iR),fA(ia),fW(iw),fR2(0),fA2(0),fW2(0),fBeta2(0),fBeta4(0),
  fMinDist(0.4),fNodeDist(0.0),fSmearing(0.0),fRecenter(1),fLattice(0),fSmax(99),
  fF(0),fTrials(0),fNonSmeared(0),fFunc1(ifunc),fFunc2(0),fFunc3(0),fNucleons(0),
  fPhiRot(0),fThetaRot(0),fNucCounter(-1),fIsUsed(0),fMaxR(14)
{
  if (fN==0) {
    cout << "Setting up nucleus " << iname << endl;
    Lookup(iname);
  }
}

TGlauNucleus::~TGlauNucleus()
{
  if (fIsUsed)
    delete fIsUsed;
  if (fNucleons)
    delete fNucleons;
  delete fFunc1;
  delete fFunc2;
  delete fFunc3;
}

void TGlauNucleus::Draw(Double_t xs, Int_t colp, Int_t cols)
{
  Double_t r = 0.5*TMath::Sqrt(xs/TMath::Pi()/10.);
  TEllipse en;
  en.SetLineStyle(1);
  en.SetLineWidth(1);
  en.SetFillStyle(1001);
  for (Int_t i = 0; i<fNucleons->GetEntries(); ++i) {
    TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
    if (!gn->IsSpectator()) {
      en.SetFillColor(colp);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    } else {
      en.SetFillColor(cols);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    }
  }
}

void TGlauNucleus::Lookup(const char* name)
{
  SetName(name);
  TString tmp(name);
  Double_t r0=0, r1=0, r2=0;

  if      (TString(name) == "p")       {fN = 1;   fR = 0.234;      fA = 0;      fW =  0;       fF = 0;  fZ=1;}
  else if (TString(name) == "pg")      {fN = 1;   fR = 0.514;      fA = 0;      fW =  0;       fF = 9;  fZ=1;} 
  else if (TString(name) == "pdg")     {fN = 1;   fR = 1;          fA = 0;      fW =  0;       fF = 10; fZ=1;} // from arXiv:1101.5953
  else if (TString(name) == "dpf")     {fN = 2;   fR = 0.01;       fA = 0.5882; fW =  0;       fF = 1;  fZ=1;} // deuteron 2pf (tuned to Hulthen)
  else if (TString(name) == "dh")      {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 3;  fZ=1;} // deuteron Hulthen free
  else if (TString(name) == "d")       {fN = 2;   fR = 0.2283;     fA = 1.1765; fW =  0;       fF = 4;  fZ=1;} // deuteron Hulthen constrained
  else if (TString(name) == "He3")     {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=1;} // read configurations from file
  else if (TString(name) == "H3")      {fN = 3;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "He4")     {fN = 4;   fR = 0.00;       fA = 0.0000; fW =  0;       fF = 6;  fZ=2;} // read configurations from file
  else if (TString(name) == "C")       {fN = 12;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 6;  fZ=6;} // read configurations from file  
  else if (TString(name) == "O")       {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 6;  fZ=8;} // read configurations from file
  else if (TString(name) == "Opar")    {fN = 16;  fR = 2.608;      fA = 0.513;  fW = -0.051;   fF = 1;  fZ=8;} // WS parameterization
  else if (TString(name) == "Oho")     {fN = 16;  fR = 1.833;      fA = 1.544;  fW =  0;       fF = 15; fZ=8;} // Harmonic oscillator parameterization
  else if (TString(name) == "Al")      {fN = 27;  fR = 3.34;       fA = 0.580;  fW = 0.0;      fF = 8;  fZ=13; fBeta2=-0.448; fBeta4=0.239;}
  else if (TString(name) == "Si")      {fN = 28;  fR = 3.34;       fA = 0.580;  fW = -0.233;   fF = 1;  fZ=14;}
  else if (TString(name) == "Si2")     {fN = 28;  fR = 3.34;       fA = 0.580;  fW =  0;       fF = 8;  fZ=14; fBeta2=-0.478; fBeta4=0.250;}
  else if (TString(name) == "S")       {fN = 32;  fR = 2.54;       fA = 2.191;  fW =  0.16;    fF = 2;  fZ=16;}
  else if (TString(name) == "Ar")      {fN = 40;  fR = 3.53;       fA = 0.542;  fW =  0;       fF = 1;  fZ=18;}
  else if (TString(name) == "Ca")      {fN = 40;  fR = 3.766;      fA = 0.586;  fW = -0.161;   fF = 1;  fZ=20;}
  else if (TString(name) == "Ni")      {fN = 58;  fR = 4.309;      fA = 0.517;  fW = -0.1308;  fF = 1;  fZ=28;}
  else if (TString(name) == "Cu")      {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 1;  fZ=29;}
  else if (TString(name) == "Curw ")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 12; fZ=29; r0=1.00898; r1=-0.000790403; r2=-0.000389897;} 
  else if (TString(name) == "Cu2")     {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 8;  fZ=29; fBeta2=0.162; fBeta4=-0.006;}  
  else if (TString(name) == "Cu2rw")   {fN = 63;  fR = 4.20;       fA = 0.596;  fW =  0;       fF = 14; fZ=29; fBeta2=0.162; fBeta4=-0.006; r0=1.01269; r1=-0.00298083; r2=-9.97222e-05;}  
  else if (TString(name) == "CuHN")    {fN = 63;  fR = 4.28;       fA = 0.5;    fW =  0;       fF = 1;  fZ=29;} // from arXiv:0904.4080v1
  else if (TString(name) == "Xe")      {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 1;  fZ=54;} // adapted from arXiv:1703.04278
  else if (TString(name) == "Xes")     {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 1;  fZ=54;} // scale from Sb (Antimony, A=122, r=5.32) by 1.019 = (129/122)**0.333
  else if (TString(name) == "Xe2")     {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.161; fBeta4=-0.003;} // adapted from arXiv:1703.04278 and Z. Physik (1974) 270: 113
  else if (TString(name) == "Xe2a")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 8;  fZ=54; fBeta2=0.18; fBeta4=0;} // ALICE parameters (see public note from 2018 at https://cds.cern.ch/collection/ALICE%20Public%20Notes?ln=en)
  else if (TString(name) == "Xerw")    {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 12; fZ=54; r0=1.00911; r1=-0.000722999; r2=-0.0002663;}
  else if (TString(name) == "Xesrw")   {fN = 129; fR = 5.42;       fA = 0.57;   fW =  0;       fF = 12; fZ=54; r0=1.0096; r1=-0.000874123; r2=-0.000256708;}
  else if (TString(name) == "Xe2arw")  {fN = 129; fR = 5.36;       fA = 0.59;   fW =  0;       fF = 14; fZ=54; fBeta2=0.18; fBeta4=0; r0=1.01246; r1=-0.0024851; r2=-5.72464e-05;} 
  else if (TString(name) == "W")       {fN = 186; fR = 6.58;       fA = 0.480;  fW =  0;       fF = 1;  fZ=74;}
  else if (TString(name) == "Au")      {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 1;  fZ=79;}
  else if (TString(name) == "Aurw")    {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 12; fZ=79; r0=1.00899; r1=-0.000590908; r2=-0.000210598;}
  else if (TString(name) == "Au2")     {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 8;  fZ=79; fBeta2=-0.131; fBeta4=-0.031; }
  else if (TString(name) == "Au2rw")   {fN = 197; fR = 6.38;       fA = 0.535;  fW =  0;       fF = 14; fZ=79; fBeta2=-0.131; fBeta4=-0.031; r0=1.01261; r1=-0.00225517; r2=-3.71513e-05;}
  else if (TString(name) == "AuHN")    {fN = 197; fR = 6.42;       fA = 0.44;   fW =  0;       fF = 1;  fZ=79;} // from arXiv:0904.4080v1
  else if (TString(name) == "Pb")      {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbrw")    {fN = 208; fR = 6.62;       fA = 0.546;  fW =  0;       fF = 12; fZ=82; r0=1.00863; r1=-0.00044808; r2=-0.000205872;} //only Pb 207 was tested but should be the same for 208
  else if (TString(name) == "Pb*")     {fN = 208; fR = 6.624;      fA = 0.549;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "PbHN")    {fN = 208; fR = 6.65;       fA = 0.460;  fW =  0;       fF = 1;  fZ=82;}
  else if (TString(name) == "Pbpn")    {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 11; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  else if (TString(name) == "Pbpnrw")  {fN = 208; fR = 6.68;       fA = 0.447;  fW =  0;       fF = 13; fZ=82; fR2=6.69; fA2=0.56; fW2=0;}
  // Uranium description taken from Heinz & Kuhlman, nucl-th/0411054.  In this code, fR is defined as 6.8*0.91, fW=6.8*0.26
  else if (TString(name) == "U")       {fN = 238; fR = 6.188;      fA = 0.54;   fW =  1.77;    fF = 5;  fZ=92;}  
  else if (TString(name) == "U2")      {fN = 238; fR = 6.67;       fA = 0.44;   fW =  0;       fF = 8;  fZ=92; fBeta2=0.280; fBeta4=0.093;}
  else {
    cout << "Could not find nucleus " << name << endl;
    return;
  }

  switch (fF) {
    case 0: // Proton exp
      fFunc1 = new TF1(name,"x*x*exp(-x/[0])",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 2: // 3pG
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      break;
    case 3: // Hulthen (see nucl-ex/0603010)
    case 4: // same but constrain the neutron opposite to the proton event-by-event
      fFunc1 = new TF1(name,"x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;
    case 5: // Ellipsoid (Uranium)
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,0); // same as 3pF but setting W to zero
      break;
    case 6: // He3/H3
      fFunc1 = 0; // read in file instead
      break;
    case 7: // Deformed nuclei, box method
#ifndef HAVE_MATHMORE
      cerr << "Need libMathMore.so for deformed nuclei" << endl;
      gSystem->Exit(123);
#endif
     fFunc1 = 0; // no func: only need beta parameters and use uniform box distribution
      break;
    case 8: // Deformed nuclei, TF2 method
      fFunc3 = new TF2(name,"x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4);
      break;
    case 9: // Proton gaus
      fFunc1 = new TF1(name,"x*x*exp(-x*x/[0]/[0]/2)",0,5);
      fFunc1->SetParameter(0,fR);
      break;
    case 10: // Proton dgaus
      fFunc1 = new TF1(name,"x*x*((1-[0])/[1]^3*exp(-x*x/[1]/[1])+[0]/(0.4*[1])^3*exp(-x*x/(0.4*[1])^2))",0,5);
      fFunc1->SetParameter(0,0.5);
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // 3pF for proton and neutrons
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW);
      fFunc2 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2);
      break;
    case 12: // reweighted
      fFunc1 = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,r0,r1,r2); 
      fRecenter=1;
      fSmax=0.1;
      break;
    case 13: // Pb for proton and neutrons reweighted
      fFunc1 = new TF1(Form("%s_prot",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA,fW,1.00866,-0.000461484,-0.000203571);
      fFunc2 = new TF1(Form("%s_neut",name),"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))/([3]+[4]*x+[5]*x^2)",0,fMaxR);
      fFunc2->SetParameters(fR2,fA2,fW2,1.00866,-0.000461484,-0.000203571);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 14: // Deformed nuclei, TF2 method, reweighted
      fFunc3 = new TF2(name,"x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))/([4]+[5]*x+[6]*x^2)",0,fMaxR,0.0,TMath::Pi());
      fFunc3->SetNpx(120);
      fFunc3->SetNpy(120);
      fFunc3->SetParameters(fR,fA,fBeta2,fBeta4,r0,r1,r2);
      fRecenter=1;
      fSmax=0.1;
      break;
    case 15: // harmonic oscillator model 
      fFunc1 = new TF1(name,"x^2*(1+[0]*(x^2/[1]^2))*exp(-x^2/[1]^2)",0,fMaxR);
      fFunc1->SetParameters(fR,fA);
      break;
    default:
      cerr << "Could not find function type " << fF << endl;
  }
  return;
}

void TGlauNucleus::SetA(Double_t ia, Double_t ia2)
{
  fA  = ia;
  fA2 = ia2;
  switch (fF) {
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(1,fA);
      break;
    case 8:
      fFunc3->SetParameter(1,fA);
      break;
    case 11: //p&n
      fFunc1->SetParameter(1,fA);//proton
      fFunc2->SetParameter(1,fA2);//neutron
      break;
    default:
      cout << "Error: fA not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetBeta(Double_t b2, Double_t b4) 
{
  fBeta2=b2; 
  fBeta4=b4;      
  if (fFunc3) {
    fFunc3->SetParameter(2,fBeta2);
    fFunc3->SetParameter(3,fBeta4);
  }
}

void TGlauNucleus::SetR(Double_t ir, Double_t ir2)
{
  fR  = ir;
  fR2 = ir2;
  switch (fF) {
    case 0:  // Proton exp
    case 9:  // Proton gaus
    case 1:  // 3pF
    case 12: // 3pF with pol2 normalization
    case 2:  // 3pG
    case 5:  // Ellipsoid (Uranium)
      fFunc1->SetParameter(0,fR);
      break;
    case 8:
      fFunc3->SetParameter(0,fR);
      break;
    case 10: // Proton
      fFunc1->SetParameter(1,fR);
      break;
    case 11: // p&n
      fFunc1->SetParameter(0,fR);//proton
      fFunc2->SetParameter(0,fR2);//neutron
      break;
    default:
      cout << "Error: fR not needed for function " << fF <<endl;
  }
}

void TGlauNucleus::SetW(Double_t iw)
{
  fW = iw;
  switch (fF) {
    case 1: // 3pF
    case 2: // 3pG
      fFunc1->SetParameter(2,fW);
      break;
    default:
      cout << "Error: fW not needed for function " << fF <<endl;
  }
}

Bool_t TGlauNucleus::TestMinDist(Int_t n, Double_t x, Double_t y, Double_t z) const
{
  if (fMinDist<=0)
    return kTRUE;
  const Double_t md2 = fMinDist*fMinDist; 
  for (Int_t j = 0; j<n; ++j) {
    TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
    Double_t xo=other->GetX();
    Double_t yo=other->GetY();
    Double_t zo=other->GetZ();
    Double_t dist2 = (x-xo)*(x-xo)+
		                 (y-yo)*(y-yo)+
		                 (z-zo)*(z-zo);
    if (dist2<md2) {
      return kFALSE;
    }
  }
  return kTRUE;
}

TVector3 &TGlauNucleus::ThrowNucleons(Double_t xshift)
{

  if (fNucleons==0) {
    fNucleons=new TObjArray(fN);
    fNucleons->SetOwner();
    for (Int_t i=0; i<fN; ++i) {
      TGlauNucleon *nucleon=new TGlauNucleon(); 
      nucleon->SetType(0);
      if (i<fZ) 
        nucleon->SetType(1);
      fNucleons->Add(nucleon); 
    }
  } 
  if (1) { //randomize p and n in nucleus
    for (Int_t i=0,iz=0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)fNucleons->At(i);
      Double_t frac=double(fZ-iz)/double(fN-i);
      Double_t rn=gRandom->Uniform(0,1);
      if (rn<frac) {
        nucleon->SetType(1);
        ++iz;
      } else {
        nucleon->SetType(0);
      }
    }
  }

 cmscheck: /* start over here in case shift was too large */

  fTrials = 0;
  fNonSmeared = 0;
  fPhiRot = gRandom->Rndm()*2*TMath::Pi();
  const Double_t cosThetaRot = 2*gRandom->Rndm()-1;
  fThetaRot = TMath::ACos(cosThetaRot);
  fXRot = gRandom->Rndm()*2*TMath::Pi();
  fYRot = gRandom->Rndm()*2*TMath::Pi();
  fZRot = gRandom->Rndm()*2*TMath::Pi();

  const Bool_t hulthen = (fF==3||fF==4);
  TString tmpname(GetName());
  Bool_t nucleonsfromfile = false;
  if ((tmpname=="He3") || (tmpname=="H3") ||
      (tmpname=="He4") || (tmpname=="C")   || 
      (tmpname=="O"))
    nucleonsfromfile = true;  
  
  if (fN==1) { //special treatment for proton
    Double_t r = fFunc1->GetRandom();
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(0));
    nucleon->Reset();
    nucleon->SetXYZ(r * stheta * TMath::Cos(phi),
		    r * stheta * TMath::Sin(phi),
		    r * ctheta);
    fTrials = 1;

  } else if (fN==2 && hulthen) { //special treatment for Hulten

    Double_t r = fFunc1->GetRandom()/2;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi();
    Double_t ctheta = 2*gRandom->Rndm() - 1;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);

    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    nucleon1->Reset();
    nucleon1->SetXYZ(r * stheta * TMath::Cos(phi),
                     r * stheta * TMath::Sin(phi),
                     r * ctheta);
    nucleon2->Reset();
    if (fF==4) { // place opposite of 1
      nucleon2->SetXYZ(-nucleon1->GetX(),
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());
    } else {
      r = fFunc1->GetRandom()/2;
      phi = gRandom->Rndm() * 2 * TMath::Pi();
      ctheta = 2*gRandom->Rndm() - 1;
      stheta = TMath::Sqrt(1-ctheta*ctheta);
      nucleon2->SetXYZ(r * stheta * TMath::Cos(phi),
		       r * stheta * TMath::Sin(phi),
		       r * ctheta);
    }
    fTrials = 1;

  } else if (fN > 2 && fN < 20 && nucleonsfromfile) { 

    // if the first call, then read in the file configurations
    if (fNucCounter == -1) {
      // read in the ascii file into the array and step through the counter
      char filename[100] = "foo.dat";
      if (tmpname=="He3") {
        sprintf(filename,"he3_plaintext.dat");
      } else if (tmpname=="H3") {
        sprintf(filename,"h3_plaintext.dat");
      } else if (tmpname=="He4") {
        sprintf(filename,"he4_plaintext.dat");
      } else if (tmpname=="C") {
        sprintf(filename,"carbon_plaintext.dat");
      } else if (tmpname=="O") {
        sprintf(filename,"oxygen_plaintext.dat");
      }
      cout << "Reading in " << filename << " for nucleon configurations with fN = " << fN << endl;
      ifstream myfile;
      myfile.open(filename);
      if (!myfile) {
        cout << "ERROR:  no file for nucleon configurations found with name = " << filename << endl;
        gSystem->Exit(123);
      }

      Int_t inputcounter = 0;
      while (myfile) {
        if (inputcounter > 5999) break;
        Double_t foo;
	if (fN == 3) {
	  myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
		 >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
		 >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
		 >> foo >> foo >> foo >> foo;
	} else if (fN == 4) {
	  // no extra data with isospin information at the end of the nucleon configurations
	  myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
		 >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
		 >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
		 >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2];
	} else if (fN == 12) {
	  // no extra data with isospin information at the end of the nucleon configurations
	  // two extra words at the beginning --> foo foo
	  myfile >> foo >> foo 
	         >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
		 >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
		 >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
		 >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2]
		 >> fNucArr[inputcounter][4][0] >> fNucArr[inputcounter][4][1] >> fNucArr[inputcounter][4][2]
		 >> fNucArr[inputcounter][5][0] >> fNucArr[inputcounter][5][1] >> fNucArr[inputcounter][5][2]
		 >> fNucArr[inputcounter][6][0] >> fNucArr[inputcounter][6][1] >> fNucArr[inputcounter][6][2]
		 >> fNucArr[inputcounter][7][0] >> fNucArr[inputcounter][7][1] >> fNucArr[inputcounter][7][2]
		 >> fNucArr[inputcounter][8][0] >> fNucArr[inputcounter][8][1] >> fNucArr[inputcounter][8][2]
		 >> fNucArr[inputcounter][9][0] >> fNucArr[inputcounter][9][1] >> fNucArr[inputcounter][9][2]
		 >> fNucArr[inputcounter][10][0] >> fNucArr[inputcounter][10][1] >> fNucArr[inputcounter][10][2]
		 >> fNucArr[inputcounter][11][0] >> fNucArr[inputcounter][11][1] >> fNucArr[inputcounter][11][2];
	} else if (fN == 16) {
	  // no extra data with isospin information at the end of the nucleon configurations
	  myfile >> fNucArr[inputcounter][0][0] >> fNucArr[inputcounter][0][1] >> fNucArr[inputcounter][0][2]
		 >> fNucArr[inputcounter][1][0] >> fNucArr[inputcounter][1][1] >> fNucArr[inputcounter][1][2]
		 >> fNucArr[inputcounter][2][0] >> fNucArr[inputcounter][2][1] >> fNucArr[inputcounter][2][2]
		 >> fNucArr[inputcounter][3][0] >> fNucArr[inputcounter][3][1] >> fNucArr[inputcounter][3][2]
		 >> fNucArr[inputcounter][4][0] >> fNucArr[inputcounter][4][1] >> fNucArr[inputcounter][4][2]
		 >> fNucArr[inputcounter][5][0] >> fNucArr[inputcounter][5][1] >> fNucArr[inputcounter][5][2]
		 >> fNucArr[inputcounter][6][0] >> fNucArr[inputcounter][6][1] >> fNucArr[inputcounter][6][2]
		 >> fNucArr[inputcounter][7][0] >> fNucArr[inputcounter][7][1] >> fNucArr[inputcounter][7][2]
		 >> fNucArr[inputcounter][8][0] >> fNucArr[inputcounter][8][1] >> fNucArr[inputcounter][8][2]
		 >> fNucArr[inputcounter][9][0] >> fNucArr[inputcounter][9][1] >> fNucArr[inputcounter][9][2]
		 >> fNucArr[inputcounter][10][0] >> fNucArr[inputcounter][10][1] >> fNucArr[inputcounter][10][2]
		 >> fNucArr[inputcounter][11][0] >> fNucArr[inputcounter][11][1] >> fNucArr[inputcounter][11][2]
		 >> fNucArr[inputcounter][12][0] >> fNucArr[inputcounter][12][1] >> fNucArr[inputcounter][12][2]
		 >> fNucArr[inputcounter][13][0] >> fNucArr[inputcounter][13][1] >> fNucArr[inputcounter][13][2]
		 >> fNucArr[inputcounter][14][0] >> fNucArr[inputcounter][14][1] >> fNucArr[inputcounter][14][2]
		 >> fNucArr[inputcounter][15][0] >> fNucArr[inputcounter][15][1] >> fNucArr[inputcounter][15][2];
	}
        ++inputcounter;
      }
      myfile.close();
      fNucCounter=0;
    } // done reading in the file the first time

    if (fNucCounter > 5999) 
      fNucCounter = 0;

    // change to loop over fN nucleons!
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      nucleon->SetXYZ(fNucArr[fNucCounter][i][0],
		      fNucArr[fNucCounter][i][1],
		      fNucArr[fNucCounter][i][2]);
      nucleon->RotateXYZ(fPhiRot,fThetaRot);
    }

    ++fNucCounter;
    fTrials = 1;
  } else { // all other nuclei 

    const Double_t startingEdge  = 20; // throw nucleons within a cube of this size (fm)
    const Double_t startingEdgeX = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeY = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Double_t startingEdgeZ = startingEdge + fNodeDist*gRandom->Rndm() - 0.5*fNodeDist;
    const Int_t nslots = 2*startingEdge/fNodeDist+1;
    if (fNodeDist>0) {
      if (fMinDist>fNodeDist) {
        cout << "Minimum distance (nucleon hard core diameter) [" 
          << fMinDist << "] cannot be larger than the nodal spacing of the grid [" 
          << fNodeDist << "]." << endl;
        cout << "Quitting...." << endl;
        gSystem->Exit(123);
      }
      if (!fIsUsed)
        fIsUsed = new TBits(nslots*nslots*nslots);
      else
        fIsUsed->ResetAllBits();
    }
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while (1) {
        ++fTrials;
        Bool_t nucleon_inside = 0;
        Double_t x=999, xsmeared=999;
        Double_t y=999, ysmeared=999;
        Double_t z=999, zsmeared=999;
        if (fF==5||fF==7) { // the extended way, throw in a box and test the weight
          while (!nucleon_inside) {
            x = (fR*2)*(gRandom->Rndm() * 2 - 1);
            y = (fR*2)*(gRandom->Rndm() * 2 - 1);
            z = (fR*2)*(gRandom->Rndm() * 2 - 1);
            Double_t r = TMath::Sqrt(x*x+y*y);
            Double_t theta = TMath::ATan2(r,z);
            Double_t R = TMath::Sqrt(x*x+y*y+z*z);
            Double_t Rtheta = fR;
            if (fF==5)
              Rtheta= fR + fW*TMath::Cos(theta)*TMath::Cos(theta);
            if (fF==7)
#ifdef HAVE_MATHMORE
              Rtheta = fR*(1+fBeta2*ROOT::Math::sph_legendre(2,0,theta)+fBeta4*ROOT::Math::sph_legendre(4,0,theta));
#else
            cerr << "Should not end here because you do not have libMathMore" << endl;
#endif
            Double_t prob = 1/(1+TMath::Exp((R-Rtheta)/fA));
            if (gRandom->Rndm()<prob) 
              nucleon_inside=1;
          }
        } else if ((fF==8) || (fF==14)) { // use TF2
          Double_t r;
          Double_t theta;
          fFunc3->GetRandom2(r,theta);
          Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
          x = r * TMath::Sin(phi) * TMath::Sin(theta);
          y = r * TMath::Cos(phi) * TMath::Sin(theta);
          z = r *                   TMath::Cos(theta);
        } else { // all other types
          TF1 *ff = fFunc1;
          if ((fFunc2) && (nucleon->GetType()==0))
            ff = fFunc2;
          if (fNodeDist<=0) { // "continuous" mode
            Double_t r = ff->GetRandom();
            Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
            Double_t ctheta = 2*gRandom->Rndm() - 1 ;
            Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
            x = r * stheta * TMath::Cos(phi);
            y = r * stheta * TMath::Sin(phi);
            z = r * ctheta;
          } else { // "grid/lattice" mode
            Int_t iNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t jNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t kNode = Int_t((2*startingEdge/fNodeDist)*gRandom->Rndm());
            Int_t index=iNode*nslots*nslots+jNode*nslots+kNode;
            if (fIsUsed->TestBitNumber(index))
              continue;
            if (fLattice==1) {       // Primitive cubic system (PCS) -> https://en.wikipedia.org/wiki/Cubic_crystal_system
              x = fNodeDist*(iNode) - startingEdgeX;
              y = fNodeDist*(jNode) - startingEdgeY;
              z = fNodeDist*(kNode) - startingEdgeZ;
            } else if (fLattice==2) { //Body centered cubic (BCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(-iNode+jNode+kNode) - 0.5*startingEdgeX;
              y = 0.5*fNodeDist*(+iNode-jNode+kNode) - 0.5*startingEdgeY;
              z = 0.5*fNodeDist*(+iNode+jNode-kNode) - 0.5*startingEdgeZ;
            } else if (fLattice==3) { //Face Centered Cubic (FCC) -> http://mathworld.wolfram.com/CubicClosePacking.html
              x = 0.5*fNodeDist*(jNode+kNode) - startingEdgeX;
              y = 0.5*fNodeDist*(iNode+kNode) - startingEdgeY;
              z = 0.5*fNodeDist*(iNode+jNode) - startingEdgeZ;
            } else {                  //Hexagonal close packing (HCP) -> https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
              x = 0.5*fNodeDist*(2*iNode+((jNode+kNode)%2))          - startingEdgeX;
              y = 0.5*fNodeDist*(TMath::Sqrt(3)*(jNode+(kNode%2)/3)) - startingEdgeY;
              z = 0.5*fNodeDist*(kNode*2*TMath::Sqrt(6)/3)           - startingEdgeZ;
            }
            const Double_t r2 = x*x + y*y + z*z;
            const Double_t r  = TMath::Sqrt(r2);
	    if ((r>fMaxR)||(r2*gRandom->Rndm()>ff->Eval(r)))
	      continue;
            if (fSmearing>0.0) {
              Int_t nAttemptsToSmear = 0;
              while (1) {
                xsmeared = x*gRandom->Gaus(1.0,fSmearing);
                ysmeared = y*gRandom->Gaus(1.0,fSmearing);
                zsmeared = z*gRandom->Gaus(1.0,fSmearing);
                nAttemptsToSmear++;
                if (TestMinDist(i,xsmeared,ysmeared,zsmeared)) {
                  x = xsmeared;
                  y = ysmeared;
                  z = zsmeared;
                  break;
                }
                if (nAttemptsToSmear>=99) {
                  cerr << "Could not place on this node :: [" << x <<","<< y <<","<< z <<"] r = " << TMath::Sqrt(x*x+y*y+z*z) << " fm; "
                    << "Node (" << iNode << "," << jNode << "," << kNode << ") not smeared !!!" << endl;
                  ++fNonSmeared;
                  break;
                }
              }
            }
            fIsUsed->SetBitNumber(index);
          } /* end "grid/lattice mode" */
	}
	nucleon->SetXYZ(x,y,z);
	if (fF==5||fF==7||fF==8||fF==14) 
	  nucleon->RotateXYZ(fPhiRot,fThetaRot); // Uranium etc.
	if (fNodeDist>0) {
	  nucleon->RotateXYZ_3D(fXRot,fYRot,fZRot);
	  break;
	}
	if (TestMinDist(i,x,y,z))
	  break;
      }
    }
  }    

  // calculate center of mass
  Double_t sumx=0;       
  Double_t sumy=0;       
  Double_t sumz=0;       
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    sumx += nucleon->GetX();
    sumy += nucleon->GetY();
    sumz += nucleon->GetZ();
  }
  sumx = sumx/fN;
  sumy = sumy/fN;
  sumz = sumz/fN;

  static TVector3 finalShift;
  finalShift.SetXYZ(sumx,sumy,sumz);
  if (finalShift.Mag()>fSmax)
    goto cmscheck;
  Double_t fsumx = 0;
  Double_t fsumy = 0;
  Double_t fsumz = 0;
  if (fRecenter==1) {
    fsumx = sumx;
    fsumy = sumy;
    fsumz = sumz;
  } else if (fRecenter==2) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(fN-1));
    Double_t x = nucleon->GetX() - fN*sumx;
    Double_t y = nucleon->GetY() - fN*sumy;
    Double_t z = nucleon->GetZ() - fN*sumz;
    nucleon->SetXYZ(x,y,z);
  } else if ((fRecenter==3)||(fRecenter==4)) {
    TVector3 zVec;
    zVec.SetXYZ(0,0,1);
    TVector3 shiftVec;
    shiftVec.SetXYZ(sumx,sumy,sumz);
    TVector3 orthVec;
    orthVec = shiftVec.Cross(zVec);
    TRotation myRot;
    myRot.Rotate(shiftVec.Angle(zVec),orthVec);
    TVector3 myNuc;
    for (Int_t i = 0; i<fN; ++i) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      myNuc.SetXYZ(nucleon->GetX(),nucleon->GetY(),nucleon->GetZ());
      myNuc.Transform(myRot);
      nucleon->SetXYZ(myNuc.X(), myNuc.Y(), myNuc.Z());
    }
    if (fRecenter==3)
      fsumz = shiftVec.Mag();
  }

  // recenter and shift
  for (Int_t i = 0; i<fN; ++i) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    nucleon->SetXYZ(nucleon->GetX()-fsumx + xshift,
		    nucleon->GetY()-fsumy,
		    nucleon->GetZ()-fsumz);
  }

  return finalShift;

}

//---------------------------------------------------------------------------------
  ClassImp(TGlauberMC)
  ClassImp(TGlauberMC::Event)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect, Double_t xsectsigma) :
  fANucleus(NA),fBNucleus(NB),
  fXSect(xsect),fXSectOmega(0),fXSectLambda(0),fXSectEvent(0),
  fNucleonsA(0),fNucleonsB(0),fNucleons(0),
  fAN(0),fBN(0),fNt(0),
  fEvents(0),fTotalEvents(0),fBmin(0),fBmax(20),fHardFrac(0.65),
  fDetail(99),fCalcArea(0),fCalcLength(0), fDoCore(0), fDoAAGG(1),
  fMaxNpartFound(0),f2Cx(0),fPTot(0),fNNProf(0),
  fEv()
{
  if (xsectsigma>0) {
    fXSectOmega = xsectsigma;
    fXSectLambda = 1;
    fPTot = new TF1("fPTot","((x/[2])/(x/[2]+[0]))*exp(-(((x/[2])/[0]-1 )**2)/([1]*[1]))/[2]",0,300);
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    fPTot->SetNpx(1000);
    fXSectLambda = fXSect/fPTot->GetHistogram()->GetMean();
    cout << "final lambda=" << fXSectLambda << endl;
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    cout << "final <sigma>=" << fPTot->GetHistogram()->GetMean() << endl;
  }

  TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
  SetName(name);
  SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
  // calc next event
  if (!fNucleonsA) {
    fNucleonsA = fANucleus.GetNucleons();
    fAN = fANucleus.GetN();
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
    }
  }

  if (!fNucleonsB) {
    fNucleonsB = fBNucleus.GetNucleons();
    fBN = fBNucleus.GetN();
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
    }
  }

  Double_t xsecA[999] = {0};
  Double_t xsecB[999] = {0};
  if (fPTot) {
    fXSectEvent = fPTot->GetRandom();
    if (fDoAAGG) {
      for (Int_t i = 0; i<fAN; ++i)
	xsecA[i] = fPTot->GetRandom();
      for (Int_t i = 0; i<fBN; ++i)
	xsecB[i] = fPTot->GetRandom();
    }
  } else 
    fXSectEvent = fXSect;

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
  Double_t bh = TMath::Sqrt(d2*fHardFrac);
  if (fNNProf) {
    Double_t xmin=0,xmax=0;
    fNNProf->GetRange(xmin,xmax);
    d2 = xmax*xmax;
  }

  fEv.Reset();
  memset(fBC,0,sizeof(Bool_t)*999*999);
  Int_t nc=0,nh=0;
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Bool_t tB=nucleonB->GetType();
    for (Int_t j = 0; j<fAN; ++j) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      if (fDoAAGG && fPTot)
	d2 = 0.5*(xsecA[j]+xsecB[i])/(TMath::Pi()*10);
      if (dij>d2) 
        continue;
      Double_t bij = TMath::Sqrt(dij);
      if (fNNProf) {
        Double_t val = fNNProf->Eval(bij);
        Double_t ran = gRandom->Uniform();
        if (ran>val)
          continue;
      }
      nucleonB->Collide();
      nucleonA->Collide();
      fBC[i][j] = 1;
      fEv.BNN  += bij;
      ++nc;
      if (bij<bh)
        ++nh;
      Bool_t tA=nucleonA->GetType();
      if (tA!=tB)
        ++fEv.Ncollpn;
      else if (tA==1)
        ++fEv.Ncollpp;
      else
        ++fEv.Ncollnn;
      if (nc==1) {
        fEv.X0 = (nucleonA->GetX()+nucleonB->GetX())/2;
        fEv.Y0 = (nucleonA->GetY()+nucleonB->GetY())/2;
      }
    }
  }
  fEv.B = bgen;
  ++fTotalEvents;
  if (nc>0) {
    ++fEvents;
    fEv.Ncoll     = nc;
    fEv.Nhard     = nh;
    fEv.BNN      /= nc;
    return CalcResults(bgen);
  }
  return kFALSE;
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  Double_t sumW=0;
  Double_t sumWA=0;
  Double_t sumWB=0;

  Double_t sinphi[10] = {0};
  Double_t cosphi[10] = {0};
  Double_t rn[10]     = {0};

  const Int_t kNc = fDoCore; // used later for core/corona

  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fEv.MeanXSystem  += xA;
    fEv.MeanYSystem  += yA;
    fEv.MeanXA  += xA;
    fEv.MeanYA  += yA;
    if (nucleonA->IsWounded()) {
      Double_t w = nucleonA->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonA->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartA;
      sumW   += w;
      sumWA  += w;
      fEv.MeanX  += xA * w;
      fEv.MeanY  += yA * w;
      fEv.MeanX2 += xA * xA * w;
      fEv.MeanY2 += yA * yA * w;
      fEv.MeanXY += xA * yA * w;
    }
  }

  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fEv.MeanXSystem  += xB;
    fEv.MeanYSystem  += yB;
    fEv.MeanXB  += xB;
    fEv.MeanYB  += yB;
    if (nucleonB->IsWounded()) {
      Double_t w = nucleonB->Get2CWeight(f2Cx);
      ++fEv.Npart;
      if (nucleonB->GetNColl()==1)
	++fEv.Npart0;
      ++fEv.NpartB;
      sumW   += w;
      sumWB  += w;
      fEv.MeanX  += xB * w;
      fEv.MeanY  += yB * w;
      fEv.MeanX2 += xB * xB * w;
      fEv.MeanY2 += yB * yB * w;
      fEv.MeanXY += xB * yB * w;
    }
  }
  if (fEv.Npart>0) {
    fEv.MeanX  /= sumW;
    fEv.MeanY  /= sumW;
    fEv.MeanX2 /= sumW;
    fEv.MeanY2 /= sumW;
    fEv.MeanXY /= sumW;
  } else {
    fEv.MeanX = 0;
    fEv.MeanY  = 0;
    fEv.MeanX2 = 0;
    fEv.MeanY2 = 0;
    fEv.MeanXY = 0;
  }

  if (fAN+fBN>0) {
    fEv.MeanXSystem /= (fAN + fBN);
    fEv.MeanYSystem /= (fAN + fBN);
  } else {
    fEv.MeanXSystem = 0;
    fEv.MeanYSystem = 0;
  }
  if (fAN>0) {
    fEv.MeanXA /= fAN;
    fEv.MeanYA /= fAN;
  } else {
    fEv.MeanXA = 0;
    fEv.MeanYA = 0;
  }
  if (fBN>0) {
    fEv.MeanXB /= fBN;
    fEv.MeanYB /= fBN;
  } else {
    fEv.MeanXB = 0;
    fEv.MeanYB = 0;
  }

  fEv.VarX  = fEv.MeanX2-(fEv.MeanX*fEv.MeanX);
  fEv.VarY  = fEv.MeanY2-(fEv.MeanY*fEv.MeanY);
  fEv.VarXY = fEv.MeanXY-fEv.MeanX*fEv.MeanY;
  Double_t tmpa = fEv.VarX*fEv.VarY-fEv.VarXY*fEv.VarXY;
  if (tmpa<0) 
    fEv.AreaW = -1;
  else 
    fEv.AreaW = TMath::Sqrt(tmpa);

  if (fEv.Npart>0) {
    // do full moments relative to meanX and meanY
    for (Int_t n = 1; n<10; ++n) {
      for (Int_t ia = 0; ia<fAN; ++ia) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(ia));
	if (nucleonA->GetNColl()<=kNc) 
	  continue;
	Double_t xA=nucleonA->GetX() - fEv.MeanX;
	Double_t yA=nucleonA->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xA*xA+yA*yA);
	Double_t phi = TMath::ATan2(yA,xA);
	Double_t w = n;
	if (n==1) 
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      for (Int_t ib = 0; ib<fBN; ++ib) {
        TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(ib));
	if (nucleonB->GetNColl()<=kNc) 
	  continue;
	Double_t xB=nucleonB->GetX() - fEv.MeanX;
	Double_t yB=nucleonB->GetY() - fEv.MeanY;
	Double_t r = TMath::Sqrt(xB*xB+yB*yB);
	Double_t phi = TMath::ATan2(yB,xB);
	Double_t w = n;
	if (n==1)
	  w = 3; // use r^3 weighting for Ecc1/Psi1
	Double_t rw = TMath::Power(r,w);
	cosphi[n] += rw*TMath::Cos(n*phi);
	sinphi[n] += rw*TMath::Sin(n*phi);
	rn[n] += rw;
      }
      cosphi[n] /= fEv.Npart;
      sinphi[n] /= fEv.Npart;
      rn[n] /= fEv.Npart;
      if (rn[n]>0) {
	fPsiN[n] = (TMath::ATan2(sinphi[n],cosphi[n]) + TMath::Pi())/n;
	fEccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
      } else {
	fPsiN[n] = -1;
	fEccN[n] = -1;
      }
    }
    if (!kNc) { //silly test but useful to catch errors 
      Double_t t=TMath::Sqrt(TMath::Power(fEv.VarY-fEv.VarX,2)+4.*fEv.VarXY*fEv.VarXY)/(fEv.VarY+fEv.VarX)/fEccN[2];
      if (t<0.99||t>1.01)
        cout << "Error: Expected t=1 but found t=" << t << endl;
    }
  }

  fEv.B      = bgen;
  fEv.PhiA   = fANucleus.GetPhiRot();
  fEv.ThetaA = fANucleus.GetThetaRot();
  fEv.PhiB   = fBNucleus.GetPhiRot();
  fEv.ThetaB = fBNucleus.GetThetaRot();
  fEv.Psi1   = fPsiN[1];
  fEv.Ecc1   = fEccN[1];
  fEv.Psi2   = fPsiN[2];
  fEv.Ecc2   = fEccN[2];
  fEv.Psi3   = fPsiN[3];
  fEv.Ecc3   = fEccN[3];
  fEv.Psi4   = fPsiN[4];
  fEv.Ecc4   = fEccN[4];
  fEv.Psi5   = fPsiN[5];
  fEv.Ecc5   = fEccN[5];

  if (fCalcArea) {
    const Int_t nbins=200;
    const Double_t ell=10;
    const Double_t da=2*ell*2*ell/nbins/nbins;
    const Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2
    const Double_t r2 = d2/4.;
    const Double_t mx = fEv.MeanX;
    const Double_t my = fEv.MeanY;
    TH2D areaA("hAreaA",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    TH2D areaB("hAreaB",";x (fm);y (fm)",nbins,-ell,ell,nbins,-ell,ell);
    for (Int_t i = 0; i<fAN; ++i) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      if (!nucleonA->IsWounded())
        continue;
      if (nucleonA->GetNColl()==kNc)
        continue;
      Double_t x = nucleonA->GetX()-mx;
      Double_t y = nucleonA->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaA.GetBin(xi,yi);
          Double_t val=areaA.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaA.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaA.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaA.SetBinContent(bin,1);
        }
      }
    }
    for (Int_t i = 0; i<fBN; ++i) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      if (!nucleonB->IsWounded())
        continue;
      if (nucleonB->GetNColl()==kNc)
        continue;
      Double_t x = nucleonB->GetX()-mx;
      Double_t y = nucleonB->GetY()-my;
      for (Int_t xi=1; xi<=nbins; ++xi) {
        for (Int_t yi=1; yi<=nbins; ++yi) {
          Int_t bin = areaB.GetBin(xi,yi);
          Double_t val=areaB.GetBinContent(bin);
          if (val>0)
            continue;
          Double_t dx=x-areaB.GetXaxis()->GetBinCenter(xi);
          Double_t dy=y-areaB.GetYaxis()->GetBinCenter(yi);
          if (dx*dx+dy*dy<r2)
            areaB.SetBinContent(bin,1);
        }
      }
    }
    Double_t overlap1=0;
    Double_t overlap2=0;
    for (Int_t xi=1; xi<=nbins; ++xi) {
      for (Int_t yi=1; yi<=nbins; ++yi) {
        Int_t bin = areaA.GetBin(xi,yi);
        Double_t vA=areaA.GetBinContent(bin);
        Double_t vB=areaB.GetBinContent(bin);
        if (vA>0&&vB>0)
          ++overlap1;
        if (vA>0||vB>0)
          ++overlap2;
      }
    }
    fEv.AreaO = overlap1*da;
    fEv.AreaA = overlap2*da;
  }

  if (fCalcLength) {
    const Double_t krhs = TMath::Sqrt(fXSectEvent/40./TMath::Pi());
    const Double_t ksg  = krhs/TMath::Sqrt(5);
    const Double_t kDL  = 0.1;
    TF1 rad("rad","2*pi/[0]/[0]*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,5*ksg); 
    rad.SetParameter(0,ksg);
    const Double_t minval = rad.Eval(5*ksg);
    fEv.Phi0         = gRandom->Uniform(0,TMath::TwoPi());
    Double_t kcphi0  = TMath::Cos(fEv.Phi0);
    Double_t ksphi0  = TMath::Sin(fEv.Phi0);
    Double_t x       = fEv.X0;
    Double_t y       = fEv.Y0;
    Double_t i0a     = 0;
    Double_t i1a     = 0;
    Double_t l       = 0;
    Double_t val     = CalcDens(rad,x,y);
    while (val>minval) {
      x     += kDL * kcphi0;
      y     += kDL * ksphi0;
      i0a   += val;
      i1a   += l*val;
      l+=kDL;
      val    = CalcDens(rad,x,y);
    }
    fEv.Length = 2*i1a/i0a;
  }

  if (fEv.Npart > fMaxNpartFound) 
    fMaxNpartFound = fEv.Npart;

  return kTRUE;
}

Double_t TGlauberMC::CalcDens(TF1 &prof, Double_t xval, Double_t yval) const
{
  Double_t rmin=0,rmax=0;
  prof.GetRange(rmin,rmax);
  Double_t r2max = rmax*rmax;
  Double_t ret = 0;
  for (Int_t i = 0; i<fAN; ++i) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    if (!nucleonA->IsWounded())
      continue;
    Double_t x = nucleonA->GetX();
    Double_t y = nucleonA->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    if (!nucleonB->IsWounded())
      continue;
    Double_t x = nucleonB->GetX();
    Double_t y = nucleonB->GetY();
    Double_t r2=(xval-x)*(xval-x)+(yval-y)*(yval-y);
    if (r2>r2max)
      continue;
    ret += prof.Eval(TMath::Sqrt(r2));
  }
  return ret;
}

void TGlauberMC::Draw(Option_t* option)
{
  static TH2F *h2f = new TH2F("hGlauberMC",";x (fm);y(fm)",1,-18,18,1,-12,12);

  h2f->Reset();
  h2f->SetStats(0);
  h2f->Draw();

  TEllipse e;
  e.SetFillColor(0);
  e.SetFillStyle(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
  fANucleus.Draw(fXSect, kMagenta, kYellow);
  fBNucleus.Draw(fXSect, kMagenta, kOrange);

  TString opt(option);
  if (opt.IsNull())
    return;

  Double_t sy2 = GetSy2();
  Double_t sx2 = GetSx2();
  Double_t phase = 0;
  if (sy2<sx2) {
    Double_t d = sx2;
    sx2 = sy2;
    sy2 = d;
    phase = TMath::Pi()/2.;
  }
  Double_t x1 = (0.5*(sy2-sx2)+TMath::Sqrt(TMath::Power(sy2-sx2,2.)-4*TMath::Power(GetSxy(),2)));
  Double_t ang = TMath::ATan2(-GetSxy(),x1)+phase;
  TLine l;
  l.SetLineWidth(3);
  l.DrawLine(-10*TMath::Cos(ang),-10*TMath::Sin(ang),10*TMath::Cos(ang),10*TMath::Sin(ang));
}

Double_t TGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBmax*fBmax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
  if (!fNucleonsA || !fNucleonsB) return 0;
  if (fNucleons) return fNucleons;

  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; ++i) {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; ++i) {
    allnucleons->Add(fNucleonsB->At(i));
  }
  fNucleons = allnucleons;
  return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  if (bgen<0) 
    bgen = TMath::Sqrt((fBmax*fBmax-fBmin*fBmin)*gRandom->Rndm()+fBmin*fBmin);

  fANucleus.ThrowNucleons(-bgen/2.);
  fBNucleus.ThrowNucleons(bgen/2.);
  return CalcEvent(bgen);
}

Bool_t TGlauberMC::ReadNextEvent(Bool_t calc, const char *fname)
{
  static TFile *inf = 0;
  static Int_t iev  = 0;
  if (fname) {
    cout << "ReadNextEvent: Setting up file " << fname << endl;
    delete inf;
    inf = TFile::Open(fname);
    if (!inf) 
      return 0;
    if (!fNucleonsA) {
      fANucleus.ThrowNucleons();
      fNucleonsA = fANucleus.GetNucleons();
      fAN = fANucleus.GetN();
      for (Int_t i = 0; i<fAN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
	nucleonA->SetInNucleusA();
      }
    }
    if (!fNucleonsB) {
      fBNucleus.ThrowNucleons();
      fNucleonsB = fBNucleus.GetNucleons();
      fBN = fBNucleus.GetN();
      for (Int_t i = 0; i<fBN; ++i) {
	TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
	nucleonB->SetInNucleusB();
      }
    }
    if (calc)
      return 1;
    fNt = dynamic_cast<TNtuple*>(inf->Get(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName())));
    if (!fNt) {
      cerr << "ReadNextEvent: Could not find ntuple!" << endl;
      inf->ls();
      return 0;
    }
    fNt->SetBranchAddress("Npart",&fEv.Npart);
    fNt->SetBranchAddress("Ncoll",&fEv.Ncoll);
    fNt->SetBranchAddress("B",&fEv.B);
    fNt->SetBranchAddress("BNN",&fEv.BNN);
    fNt->SetBranchAddress("VarX",&fEv.VarX);
    fNt->SetBranchAddress("VarY",&fEv.VarY);
    fNt->SetBranchAddress("VarXY",&fEv.VarXY);
    fNt->SetBranchAddress("NpartA",&fEv.NpartA);
    fNt->SetBranchAddress("NpartB",&fEv.NpartB);
    fNt->SetBranchAddress("Npart0",&fEv.Npart0);
    fNt->SetBranchAddress("Psi1",&fEv.Psi1);
    fNt->SetBranchAddress("Ecc1",&fEv.Ecc1);
    fNt->SetBranchAddress("Psi2",&fEv.Psi2);
    fNt->SetBranchAddress("Ecc2",&fEv.Ecc2);
    fNt->SetBranchAddress("Psi3",&fEv.Psi3);
    fNt->SetBranchAddress("Ecc3",&fEv.Ecc3);
    fNt->SetBranchAddress("Psi4",&fEv.Psi4);
    fNt->SetBranchAddress("Ecc4",&fEv.Ecc4);
    fNt->SetBranchAddress("Psi5",&fEv.Psi5);
    fNt->SetBranchAddress("Ecc5",&fEv.Ecc5);
    return 1;
  }
  if ((!inf)||(!fNt&&!calc)) {
    cerr << "ReadNextEvent was not initialized" <<endl;
    return 0;
  }
  TObjArray *arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  if (!arr) {
    if (iev==0) {
      cerr << "ReadNextEvent could not read nucleon array for event " << iev << endl;
      return 0;
    }
    iev = 0;
    cerr << "ReadNextEvent resetting to first event" << endl;
    arr = dynamic_cast<TObjArray*>(inf->Get(Form("nucleonarray%d",iev)));
  }

  Double_t bgenA=0, bgenB=0;
  Int_t inA=0, inB=0;
  const Int_t nNucls = arr->GetEntries();
  for (Int_t iNucl=0; iNucl<nNucls; ++iNucl) {
    TGlauNucleon *nuclinp = static_cast<TGlauNucleon*>(arr->At(iNucl));
    TGlauNucleon *nuclout = 0;
    if (nuclinp->IsInNucleusB()) { 
      nuclout = static_cast<TGlauNucleon*>(fNucleonsB->At(inB));
      bgenB += nuclinp->GetX();
      ++inB;
    } else {
      nuclout = static_cast<TGlauNucleon*>(fNucleonsA->At(inA));
      bgenA += nuclinp->GetX();
      ++inA;
    }
    nuclout->Reset();
    nuclout->SetXYZ(nuclinp->GetX(),nuclinp->GetY(),nuclinp->GetZ());
    nuclout->SetType(nuclinp->GetType());
    nuclout->SetEnergy(nuclinp->GetEnergy());
    if (!calc)
      nuclout->SetNColl(nuclinp->GetNColl());
  }
  delete arr;
  Double_t bgen = bgenB/inB-bgenA/inA;
  if (calc) {
    Bool_t ret = CalcEvent(bgen);
    if (0) 
      cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
    ++iev;
    return ret;
  }
  Int_t ret = fNt->GetEntry(iev);
  if (ret<=0) 
    return 0;
  fEccN[1]=fEv.Ecc1;
  fEccN[2]=fEv.Ecc2;
  fEccN[3]=fEv.Ecc3;
  fEccN[4]=fEv.Ecc4;
  fEccN[5]=fEv.Ecc5;
  if (0) 
    cout << iev << ": " << fEv.B << " " << fEv.Npart << " " << fEv.Ncoll << " " << fEv.Npart0 << endl;
  if (0) { // test ntuple values vs re-calculated values
    Double_t npart = fEv.Npart;
    Double_t ncoll = fEv.Ncoll;
    Double_t ecc2  = fEv.Ecc2;
    CalcEvent(bgen);
    if (npart!=fEv.Npart) 
      cout << iev << " differ in npart " << npart << " " << fEv.Npart << endl;
    if (ncoll!=fEv.Ncoll) 
      cout << iev << " differ in ncoll " << ncoll << " " << fEv.Ncoll << endl;
    if (TMath::Abs(ecc2-fEv.Ecc2)>0.001) 
      cout << iev << " differ in ecc2 " << ecc2 << " " << fEv.Ecc2 << endl;
  }
  ++iev;
  return 1;
}

void TGlauberMC::Run(Int_t nevents, Double_t b)
{
  if (fNt == 0) {
    TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
    TString title(Form("%s + %s (x-sect = %.1f mb) str %s",fANucleus.GetName(),fBNucleus.GetName(),fXSect,Str()));
    TString vars("Npart:Ncoll:Nhard:B:BNN:Ncollpp:Ncollpn:Ncollnn:VarX:VarY:VarXY:NpartA:NpartB:Npart0:AreaW");
    if (fDetail>1)
      vars+=":Psi1:Ecc1:Psi2:Ecc2:Psi3:Ecc3:Psi4:Ecc4:Psi5:Ecc5";
    if (fDetail>2)
      vars+=":AreaO:AreaA:X0:Y0:Phi0:Length";
    if (fDetail>3)
      vars+=":MeanX:MeanY:MeanX2:MeanY2:MeanXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB";
    if (fDetail>4)
      vars+=":PhiA:ThetaA:PhiB:ThetaB";
    fNt = new TNtuple(name,title,vars);
    fNt->SetDirectory(0);
    TObjArray *l = fNt->GetListOfBranches();
    for (Int_t i=0; i<l->GetEntries(); ++i) {
      TBranch *br = dynamic_cast<TBranch*>(l->At(i));
      if (br)
        br->SetCompressionLevel(9);
    }
  }

  for (Int_t i = 0; i<nevents; ++i) {
    while (!NextEvent(b)) {}
    fNt->Fill((Float_t*)(&fEv.Npart));
    if ((i>0)&&(i%100)==0) 
      cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;
  }
  if (nevents>99)
    cout << endl << "Done!" << endl;
}



