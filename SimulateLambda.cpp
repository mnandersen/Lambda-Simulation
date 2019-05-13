//To be implimented inside the ROOT framework (available for free from CERN)
#define V0MASS 1115.683
#define PIMASS 139.57018
#define PRMASS 938.27231

double trueRight = 0;
double trueLeft  = 0;

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
}//End CreateLambda Function

//Generate Pion 4-Momenta in lambda frame
TLorentzVector CreatePion(TRandom3 random) {
    double P = TMath::Sqrt((V0MASS*V0MASS - (PIMASS + PRMASS)*(PIMASS + PRMASS))*
                           (V0MASS*V0MASS - (PIMASS - PRMASS)*(PIMASS - PRMASS)))/(2*V0MASS);
    double Phi = random.Uniform(0, 2*3.14159);
    double Theta = 0;//(TMath::Pi() - TMath::ACos(ICos.GetRandom()));
    double E = TMath::Sqrt((P)*(P) + (PIMASS)*(PIMASS));
    
    double Pz = P*cos(Theta);
    double Pt = P*sin(Theta);
    double Py = Pt*sin(Phi);
    double Px = Pt*cos(Phi);
    
    TLorentzVector Pion(Px, Py, Pz, E);
    return Pion;
}//End CreatePion Function

//Generate Proton 4-Momenta in lambda frame
TLorentzVector CreateProton(TLorentzVector Pion) {
    double Px = -Pion.Px();
    double Py = -Pion.Py();
    double Pz = -Pion.Pz();
    double E = TMath::Sqrt(Pion.P()*Pion.P() + PRMASS*PRMASS);
    
    TLorentzVector Proton (Px, Py, Pz, E);
    return Proton;
}//End CreateProton Function

//Rotate Lambda Frame to align with Lab Frame; return daughter momenta in new frame
TLorentzVector RotationTransform(TLorentzVector Lambda, TLorentzVector Daughter) {
    double V0Theta = Lambda.Theta();
    double V0Phi = Lambda.Phi();
    
    TRotation rotMatrix;
    TRotation r1;
    TRotation r2;
    r1.RotateZ(-V0Phi);
    r2.RotateY(V0Theta);
    rotMatrix = r2 * r1;
    
    TVector3 rotatedDaughter3;
    rotatedDaughter3 = rotMatrix * (Daughter.Vect());
    TLorentzVector rotatedDaughter (rotatedDaughter3, Daughter.E());
    return rotatedDaughter;
}//End RotateCoordinates Function

//Transform daughter particles into lab frame
TLorentzVector LorentzTransform(TLorentzVector Lambda, TLorentzVector rotatedDaughter) {
    TVector3 LambdaBoostVector = Lambda.BoostVector();
    rotatedDaughter.Boost(LambdaBoostVector);
    TLorentzVector boostedDaughter = rotatedDaughter;
    return boostedDaughter;
}//End LorentzTransform Function

//Simulate the effect of detector resolution on measurments of daughter momenta
TLorentzVector SimiulateResolution (string includeRes, TRandom3 random, TLorentzVector DaughterInLab) {
    if (includeRes == "y") {
        DaughterInLab.SetPx(random.Gaus(DaughterInLab.Px(), .03*DaughterInLab.Px()));
        DaughterInLab.SetPy(random.Gaus(DaughterInLab.Py(), .03*DaughterInLab.Py()));
    }
    return DaughterInLab;
}//End SimulateResolution Function

//Simulate a measurement of a daughter track by storing it in an array (if it passes the cuts)
void SimulateMeasurement(string includeCuts, int momentumCuts, TLorentzVector DaughterInLab, TObjArray& daughterTracks) {
    //Cuts simulating the detector's geometry and particle kinematics
    if ((includeCuts == "y" && TMath::Abs(DaughterInLab.Eta()) > 0.5)) return;
    if ((includeCuts == "y" && DaughterInLab.Pt() < momentumCuts)) return;
    daughterTracks.Add(&DaughterInLab);
}//End SimulateMeasurement Function

double GetHelicity(TLorentzVector Lambda, TLorentzVector Proton) {
    double Helicity = Lambda.Vect()*Proton.Vect();
    return Helicity;
}

void SimulateLambda() {
    //Initialize histograms
    TFile *hOutput = new TFile("simulationHistos.root", "recreate");
    TH1F *hV0M   = new TH1F("hV0M", "Invariant Mass;Mass [MeV/c2];Counts", 500, 1000, 1600);

    //Random number generator
    TRandom3 r;
    
    TF1 ICos = TF1("ICos", "1 - (.64*x)", -1, 1);
    
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
        TObjArray piTracks;
        TObjArray prTracks;
        //Generate 10 Lambdas/Pions/Protons for each event
        for (int V0Num=0; V0Num<10; V0Num++) {
            r.SetSeed();
            TLorentzVector LambdaInLab = CreateLambda(r);
            TLorentzVector PionInLambda = CreatePion(r);
            TLorentzVector ProtonInLambda = CreateProton(PionInLambda);
            
            TLorentzVector PionInRot = RotationTransform(LambdaInLab, PionInLambda);
            TLorentzVector ProtonInRot = RotationTransform(LambdaInLab, ProtonInLambda);
            
            cout<< PionInRot.Phi() << ", " << PionInRot.Theta() << endl;
            cout<< LambdaInLab.Phi() << ", " << LambdaInLab.Theta() << endl;
            
            TLorentzVector PionInLab = LorentzTransform(LambdaInLab, PionInRot);
            TLorentzVector ProtonInLab = LorentzTransform(LambdaInLab, ProtonInRot);
            
            //Determine True Helicity
            double trueHelicity = GetHelicity(LambdaInLab, ProtonInRot);
            if (trueHelicity > 0) trueRight++;
            if (trueHelicity < 0) trueLeft++;
            
            //Simulate detector conditions if desired (i.e. Resolution and Cuts on daughter particles) and store Pion/Proton momenta
            SimiulateResolution(includeRes, r, PionInLab);
            SimiulateResolution(includeRes, r, ProtonInLab);
            
            SimulateMeasurement(includeCuts, momentumCuts, PionInLab, piTracks);
            SimulateMeasurement(includeCuts, momentumCuts, ProtonInLab, prTracks);
            
        }//End Lambda/Pion/Proton Creation Loop
        
        //Loop over Proton-Pion pairs in one event
        int prNum = prTracks.GetEntries();
        int piNum = piTracks.GetEntries();
        for (int i=0; i<prNum; i++) {
            for (int j=0; j<piNum; j++) {

                //Calculate invariant mass
                double prPx = ((TLorentzVector*) prTracks.At(i))->Px();
                double prPy = ((TLorentzVector*) prTracks.At(i))->Py();
                double prPz = ((TLorentzVector*) prTracks.At(i))->Pz();
                double prE = ((TLorentzVector*) prTracks.At(i))->E();
                
                double piPx = ((TLorentzVector*) piTracks.At(j))->Px();
                double piPy = ((TLorentzVector*) piTracks.At(j))->Py();
                double piPz = ((TLorentzVector*) piTracks.At(j))->Pz();
                double piE = ((TLorentzVector*) piTracks.At(j))->E();
                
                double V0E = TMath::Sqrt((prE + piE)*(prE + piE));
                
                TVector3 V03 ((prPx + prPx), (prPy + prPy), (prPz + prPz));
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
