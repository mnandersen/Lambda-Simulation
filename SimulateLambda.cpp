//To be implimented inside the ROOT framework (available for free from CERN)

#define V0MASS 1115.683
#define PIMASS 139.57018
#define PRMASS 938.27231

//GLOBAL DECLARATIONS
double trueRight = 0;
double trueLeft  = 0;

TLorentzVector pionArray[10];
TLorentzVector protonArray[10];

TRandom3 r;

//FUNCTION DEFINTIONS
//Generate Lambda 4-momentum in lab frame
TLorentzVector CreateLambda(TRandom3 random) {
    double y = random.Uniform(-1, 1);
    double Pt = random.Exp(150);
    double Phi = random.Uniform(0, 2*3.14159);
    
    double Px = TMath::Sqrt((Pt*Pt)/(1+tan(Phi)*tan(Phi)));
    double Py = TMath::Sqrt(Pt*Pt - Px*Px);
    double Pz = TMath::Sqrt(V0MASS*V0MASS + Pt*Pt)*TMath::SinH(y);
    double E = TMath::Sqrt(V0MASS*V0MASS + Pt*Pt)*TMath::CosH(y);
    
    TLorentzVector Lambda (Px, Py, Pz, E);
    return Lambda;
}

//Generate Pion 4-Momenta in lambda frame
TLorentzVector CreatePion(TRandom3 random) {
    TF1 ICos = TF1("ICos", "1 - (.64*x)", -1, 1);
    
    double P = TMath::Sqrt((V0MASS*V0MASS - (PIMASS + PRMASS)*(PIMASS + PRMASS))*
                           (V0MASS*V0MASS - (PIMASS - PRMASS)*(PIMASS - PRMASS)))/(2*V0MASS);
    double Phi = random.Uniform(0, 2*3.14159);
    double Theta = (TMath::Pi() - TMath::ACos(ICos.GetRandom()));
    double E = TMath::Sqrt((P)*(P) + (PIMASS)*(PIMASS));
    
    double Pz = P*cos(Theta);
    double Pt = P*sin(Theta);
    double Py = Pt*sin(Phi);
    double Px = Pt*cos(Phi);
    
    TLorentzVector Pion(Px, Py, Pz, E);
    return Pion;
}

//Generate Proton 4-Momenta in lambda frame
TLorentzVector CreateProton(TLorentzVector Pion) {
    double Px = -Pion.Px();
    double Py = -Pion.Py();
    double Pz = -Pion.Pz();
    double E = TMath::Sqrt(Pion.P()*Pion.P() + PRMASS*PRMASS);
    
    TLorentzVector Proton (Px, Py, Pz, E);
    return Proton;
}

//Rotate Lambda Frame to align with Lab Frame; return daughter momenta in new frame
TLorentzVector RotationTransform(TLorentzVector Lambda, TLorentzVector Daughter) {
    double V0Theta = Lambda.Theta();
    double V0Phi = Lambda.Phi();
    
    TRotation rotMatrix;
    TRotation r1;
    TRotation r2;
    r2.RotateZ(V0Phi);
    r1.RotateY(V0Theta);
    rotMatrix = r2 * r1;
    
    TVector3 rotatedDaughter3;
    rotatedDaughter3 = rotMatrix * (Daughter.Vect());
    TLorentzVector rotatedDaughter (rotatedDaughter3, Daughter.E());
    return rotatedDaughter;
}

//Transform daughter particles into lab frame
TLorentzVector LorentzTransform(TLorentzVector Lambda, TLorentzVector Daughter) {
    TVector3 LambdaBoostVector = Lambda.BoostVector();
    Daughter.Boost(LambdaBoostVector);
    return Daughter;
}

//Simulate the effect of detector resolution on measurments of daughter momenta
TLorentzVector SimiulateResolution (string includeRes, TRandom3 random, TLorentzVector DaughterInLab) {
    if (includeRes == "y") {
        DaughterInLab.SetPx(random.Gaus(DaughterInLab.Px(), .03*DaughterInLab.Px()));
        DaughterInLab.SetPy(random.Gaus(DaughterInLab.Py(), .03*DaughterInLab.Py()));
    }
    return DaughterInLab;
}

//Simulate a measurement of a daughter track by storing it in an array (if it passes the cuts)
void SimulateMeasurement(int index, string includeCuts, int momentumCuts, TLorentzVector DaughterInLab, TLorentzVector daughterArray[]) {
    //Cuts simulating the detector's geometry and particle kinematics
    if ((includeCuts == "y" && TMath::Abs(DaughterInLab.Eta()) > 0.5)) return;
    if ((includeCuts == "y" && DaughterInLab.Pt() < momentumCuts)) return;
    daughterArray[index] = DaughterInLab;
}

//Simple Helicity calculator
double GetHelicity(TLorentzVector Lambda, TLorentzVector Proton) {
    double Helicity = Lambda.Vect()*Proton.Vect();
    return Helicity;
}

void SimulateLambda() {
    //Initialize histograms
    TFile *hOutput = new TFile("simulationHistos.root", "recreate");
    TH1F *hV0M   = new TH1F("hV0M", "Invariant Mass;Mass [MeV/c2];Counts", 500, 1000, 1600);
    
    string includeRes;
    cout << "Include Resolution (y/n)? ";
    cin >> includeRes;
    
    string includeCuts;
    int momentumCuts;
    cout << "Include Cuts (y/n)? ";
    cin >> includeCuts;
    if (includeCuts == "y") {
        cout << "Enter where you want to cut on transverse momentum [MeV/c]: ";
        cin >> momentumCuts;
    }
    
    //Simulate Events
    for (int event=0; event<10000; event++) {
        //Generate 10 Lambdas/Pions/Protons for each event
        for (int V0Num=0; V0Num<10; V0Num++) {
            r.SetSeed();
            TLorentzVector LambdaInLab = CreateLambda(r);
            TLorentzVector PionInLambda = CreatePion(r);
            TLorentzVector ProtonInLambda = CreateProton(PionInLambda);
            
            TLorentzVector PionInRot = RotationTransform(LambdaInLab, PionInLambda);
            TLorentzVector ProtonInRot = RotationTransform(LambdaInLab, ProtonInLambda);
            
            TLorentzVector PionInLab = LorentzTransform(LambdaInLab, PionInRot);
            TLorentzVector ProtonInLab = LorentzTransform(LambdaInLab, ProtonInRot);
            
            //Determine True Helicity
            double trueHelicity = GetHelicity(LambdaInLab, ProtonInRot);
            if (trueHelicity > 0) trueRight++;
            if (trueHelicity < 0) trueLeft++;
            
            //Simulate detector conditions if desired (i.e. Resolution and Cuts on daughter particles) and store Pion/Proton momenta
            SimiulateResolution(includeRes, r, PionInLab);
            SimiulateResolution(includeRes, r, ProtonInLab);
            
            SimulateMeasurement(V0Num, includeCuts, momentumCuts, PionInLab, pionArray);
            SimulateMeasurement(V0Num, includeCuts, momentumCuts, ProtonInLab, protonArray);
        }//End Lambda/Pion/Proton Creation Loop
        
        //Loop over Proton-Pion pairs in one event
        for (int i=0; i<10; i++) {
            for (int j=0; j<10; j++) {
                //Calculate invariant mass
                double prPx = protonArray[i].Px();
                double prPy = protonArray[i].Py();
                double prPz = protonArray[i].Pz();;
                double prE = protonArray[i].E();;
                
                double piPx = pionArray[j].Px();
                double piPy = pionArray[j].Py();
                double piPz = pionArray[j].Pz();
                double piE = pionArray[j].E();
                
                double V0E = TMath::Sqrt((prE + piE)*(prE + piE));
                
                TVector3 V03 ((prPx + piPx), (prPy + piPy), (prPz + piPz));
                double V0M = TMath::Sqrt(V0E*V0E - V03.Mag2());
                hV0M->Fill(V0M);
                
            }//End pion track loop
         }//End proton track loop
    }//End Event Loop
    
    hV0M->Write();
    hOutput->Close();
    
    //Calculate Helicity results
    double HR = trueLeft/trueRight;
    double uHR = TMath::Sqrt((1/trueRight)*(TMath::Sqrt(trueLeft))*(1/trueRight)*(TMath::Sqrt(trueLeft)) + (-trueLeft/(trueRight*trueRight))*(TMath::Sqrt(trueRight))*(-trueLeft/(trueRight*trueRight))*(TMath::Sqrt(trueRight)));
    
    //Write results to file
    ofstream simulationResults;
    simulationResults.open("Simulation\ Results.txt");
    simulationResults << "True #Left-Handed: "  << trueLeft  << endl;
    simulationResults << "True #Right-Handed: " << trueRight << endl;
    simulationResults << std::setprecision(3) << "True Helicity Ratio (Left/Right): ";
    simulationResults << HR;
    simulationResults << " ±";
    simulationResults << std::setprecision(2) << uHR << endl;
    simulationResults.close();
    
}//End SimulateLambda
