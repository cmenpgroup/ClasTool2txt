#include "Riostream.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TClasTool.h"
#include "TIdentificator.h"
#include "TMath.h"
#include "TString.h"
#include "TMath.h"
#include "massConst.h"
using namespace std;

#define MAX_ELECTRONS 1
#define MAX_PIPLUS 0
#define MAX_PIMINUS 0
#define MAX_PHOTONS 3

#define PDG_PHOTON 22
#define PDG_POSITRON -11
#define PDG_ELECTRON 11
#define PDG_PIPLUS 211
#define PDG_PIMINUS -211
#define PDG_KPLUS 321
#define PDG_KMINUS -321
#define PDG_NEUTRON 2112
#define PDG_PROTON 2212

#define GEANT3_PHOTON 1
#define GEANT3_POSITRON 2
#define GEANT3_ELECTRON 3
#define GEANT3_PIPLUS 8
#define GEANT3_PIMINUS 9
#define GEANT3_KPLUS 11
#define GEANT3_KMINUS 12
#define GEANT3_NEUTRON 13
#define GEANT3_PROTON 14

#define CUT_Q2 1.0
#define CUT_W 2.0
#define CUT_NU 0.85

#define EBEAM 5.015  // e- beam energy in GeV

//declarations of functions
void PrintAnalysisTime(float tStart, float tStop);
void PrintUsage(char *processName);
int GetPID(string partName, int kind);

int main(int argc, char **argv)
{
    gROOT->SetBatch(true);
    extern char *optarg;
    int c;
    extern int optind;

    int i, j, k;
    int nRows, kind;
    int simTypes = 1; // number of simulation banks (1: data, 2: reconstructed + generated)
    int dEvents = 1000; // increment of events for processing print statement
    int MaxEvents = 0; // max. number of events to process
    int nfiles = 0; // number of processed files
    int genCtr = 0; // counter for generated (kind=1) events
    int recCtr = 0; // counter for reconstructed (kind=0) events
    int minRows = MAX_ELECTRONS + MAX_PHOTONS; // number of particles in an event to start filtering

    TString catPid;

    bool bBatchMode = false;    // quiet mode
    bool simul_key = false;  // simulation flag (true = simulation, false = data)
    int tgt_key = 1;  // intitialize target flag 1 = Carbon (default), 2 = Iron, 3 = Lead
    string target; // solid target name

    char *inFile;
    string outFile = "ClasTool2txt.txt";

    bool topology = false;
    vector<int> elecIndex;
    vector<int> gamIndex;

    TVector3 *myVertex = new TVector3(0., 0., 0.); // constructor is necessary

    float timeStart = clock(); // start time

    TClasTool *input = new TClasTool();
    input->InitDSTReader("ROOTDSTR");

    TIdentificator *t = new TIdentificator(input);

    for (i = 0; i < argc; ++i) cerr << argv[i] << " "; cerr << endl;
    while ((c = getopt(argc,argv, "o:M:D:t:Sih")) != -1 ) {
        switch (c) {
            case 'o': outFile = optarg; break;
            case 'M': MaxEvents = atoi(optarg); break;
            case 'D': dEvents = atoi(optarg); break;
            case 't': tgt_key = atoi(optarg); break;
            case 'S': simul_key = true; break;
            case 'i': bBatchMode = true; break;
            case 'h':
                PrintUsage(argv[0]);
                exit(0);
                break;

            default:
                cerr << "Unrecognized argument: " << optarg << endl;
                PrintUsage(argv[0]);
                exit(0);
                break;
        }
    }

    ofstream txtOut;
    txtOut.open(outFile.c_str());

    // check target selection
    switch(tgt_key){
        case 1: target = "C"; break;
        case 2: target = "Fe"; break;
        case 3: target = "Pb"; break;
        default: cout<<"Unknown target "<<target<<endl; exit(0); break;
    }
    cout<<"Analyzing " << target << " target data"<<endl;

    if(simul_key) simTypes = 2;

    for (i = optind; i < argc; ++i) {
        inFile = argv[i]; // process all arguments on command line.
        if (*inFile != '-') { // we have a file to process
            cout << "Analyzing file " << inFile << endl; // let user know which file is being processed
            input->Add(inFile); // read file into ClasTool object
            nfiles++; // increment file counter
        }
    }

    Long_t nEntries = (Long_t) input->GetEntries(); // get total number of events

    cout<<"Analyzing "<<nEntries<<" from "<<nfiles<< " files."<<endl; // print out stats

    input->Next(); // load the first event

    if(MaxEvents == 0) MaxEvents = nEntries; // if user does not set max. number of events, set to nEntries

    for(k=0; k < MaxEvents; k++) {
    	if (!bBatchMode && ((k % dEvents) == 0)) cerr << k << "\r"; // print the event number
      for(kind=0; kind<simTypes; kind++){
        elecIndex.clear(); // clear out the electron list
        gamIndex.clear(); // clear out the photon list
        topology = false; // init. the event topology cut

        if(kind == 0) nRows = input->GetNRows("EVNT");
        if(kind == 1) nRows = input->GetNRows("GSIM");
        if(nRows >= minRows){ // check that the minimum number of particles in event
          for (j = 0; j < nRows; j++) { // count particles in topology
            catPid = t -> GetCategorizationParticle(j,kind);
            if(catPid.EqualTo("electron")) elecIndex.push_back(j);
            if(catPid.EqualTo("gamma")) gamIndex.push_back(j);
          } // for loop for searching for particle in topology
          // check event topology
          topology = (elecIndex.size()>=MAX_ELECTRONS && gamIndex.size()>=MAX_PHOTONS);
          if(topology && t->Q2(kind) > CUT_Q2 && t->W(kind) > CUT_W) {
            if(kind==1){
              genCtr++;
              myVertex->SetXYZ(t->X(0, kind), t->Y(0, kind), t->Z(0, kind));
            }else{
              recCtr++;
              myVertex = t->GetCorrectedVert();
            }
            txtOut << (Int_t) t->NEvent() << "\t" << kind << endl;
            for (j = 0; j < nRows; j++) {
              txtOut << setw(10) << (Int_t) t->Id(j, kind)
		            << fixed << setprecision(5)
                << setw(10) << t->Betta(j, kind)
                << setprecision(3)
                << setw(10) << t->Px(j, kind)
                << setw(10) << t->Py(j, kind)
                << setw(10) << t->Pz(j, kind)
                << setw(10) << myVertex->X()
                << setw(10) << myVertex->Y()
                << setw(10) << myVertex->Z() << endl;
            } // for loop for printing
          } // if to check topology and cuts
        } // if to check nRows > minRows
      } // for loop for kind counter
      input->Next(); // load the next event
    }
    txtOut.close();
    cout<<"Reconstructed events found = " << recCtr << endl;
    cout<<"Generated events found = " << genCtr << endl;
    float timeStop = clock();
    PrintAnalysisTime(timeStart,timeStop);
    return 0;
}

void PrintUsage(char *processName)
{
    cerr << processName << " <options> <filename>\n";
    cerr << "\toptions are:\n";
    cerr << "\t-o<filename>\tROOT output file (def. = f1.root).\n";
    cerr << "\t-M#\t\tprocess maximum # of events.\n";
    cerr << "\t-D#\t\tinform user when # of events have been processed (def. = 1000).\n";
    cerr << "\t-t#\t\tTarget # of 1=C, 2=Fe, or 3=Pb (def. = 1).\n";
    cerr << "\t-S\t\tAnalyze simulation.\n";
    cerr << "\t-i\t\tquiet mode (no counter).\n";
    cerr << "\t-h\t\tprint the above" << endl;
}


void PrintAnalysisTime(float tStart, float tStop){
    //time to complete function
    float minutes = 0;
    float seconds = 0;
    minutes = (tStop - tStart)/1000000;
    minutes = (minutes)/60;
    seconds = fmod(minutes,1);
    minutes = minutes-seconds;
    seconds = seconds*60;

    if (minutes==0){
        cout<<endl<<"Completed in "<<seconds<<" seconds."<<endl<<endl;
    }
    else{
        cout<<endl<<"Completed in "<<minutes<<" minutes and "<<seconds<<" seconds."<<endl<<endl;
    }
}

int GetPID(string partName, int kind){

    int ret = 0;

    if(kind==0){
        if(partName.compare("Electron")==0){
            ret = PDG_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = PDG_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = PDG_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = PDG_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = PDG_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = PDG_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = PDG_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = PDG_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = PDG_PROTON;
        }else{
            cerr<<"GetPid(): Unknown PDG particle "<<partName.c_str()<<endl; exit(0);
        }
    }else if(kind==1){
        if(partName.compare("Electron")==0){
            ret = GEANT3_ELECTRON;
        }else if(partName.compare("Positron")==0){
            ret = GEANT3_POSITRON;
        }else if(partName.compare("Photon")==0){
            ret = GEANT3_PHOTON;
        }else if(partName.compare("PiPlus")==0){
            ret = GEANT3_PIPLUS;
        }else if(partName.compare("PiMinus")==0){
            ret = GEANT3_PIMINUS;
        }else if(partName.compare("KPlus")==0){
            ret = GEANT3_KPLUS;
        }else if(partName.compare("KMinus")==0){
            ret = GEANT3_KMINUS;
        }else if(partName.compare("Neutron")==0){
            ret = GEANT3_NEUTRON;
        }else if(partName.compare("Proton")==0){
            ret = GEANT3_PROTON;
        }else{
            cerr<<"GetPid(): Unknown GEANT3 particle "<<partName.c_str()<<endl; exit(0);
        }
    }else{
        cerr<<"GetPID: Unknown analysis channel "<<kind<<endl;
    }
    return ret;
}
