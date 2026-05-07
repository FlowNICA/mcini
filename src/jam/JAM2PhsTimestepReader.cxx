#include <iostream>
#include <sstream>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <fstream>

#include <jam/JAMReader.h>

#include <URun.h>
#include <UEvent.h>
#include <UParticle.h>
#include <EventInitialState.h>

namespace {
using std::ifstream;
using std::stringstream;
using std::string;
using std::cout;
using std::endl;
}

enum class ELineType {
  kRunHeader,
  kEventHeader,
  kParticle,
  kNotInitializedYet,
  kUnknown
};

TFile *gOutputFile{nullptr};
TTree *gEventTree{nullptr};

unsigned int countWordsInTheString(const std::string &str) {
  std::stringstream stream(str);
  std::string word;
  unsigned int words = 0;
  while (stream >> word) {
    ++words;
  }
  return words;
}

// ------------------------------------------------------------------
// JAM2 does not provide JAMRUN.DAT.
// We still create URun in the output file, but fill it with zeros for now.
// ------------------------------------------------------------------
void initURunWithZeroes() {
  getEntity<URun>()->SetNEvents(0);
  // For now filled with zeroes
}

// ------------------------------------------------------------------
// JAM2 event header format (from fill()):
//   # iev  nv  ncollG  npartG  b  npart  ncoll  ncollBB
// Stored into:
//   UEvent: id, b
//   EventInitialState: id, NPart = npart, NColl = ncoll
// ------------------------------------------------------------------
void parseEventHeader(const string &line) {
  stringstream stream(line);

  string hashStr;
  stream >> hashStr;
  if (hashStr != "#") {
    throw std::runtime_error("Malformed event header: must start with hash");
  }

  unsigned int eventId{0};
  unsigned int nv{0};

  unsigned int ncollG{0};
  unsigned int npartG{0};

  double impactParameter{0.0};

  unsigned int npart{0};
  unsigned int ncoll{0};
  unsigned int ncollBB{0};

  stream >> eventId >> nv >> ncollG >> npartG
         >> impactParameter >> npart >> ncoll >> ncollBB;

  std::cout << "working on event#" << eventId << std::endl;

  getEntity<UEvent>()->Clear();
  getEntity<UEvent>()->SetParameters(eventId, (float)impactParameter, 0., 0, 0, 0.);

  getEntity<EventInitialState>()->clear();
  getEntity<EventInitialState>()->setId(eventId);
  getEntity<EventInitialState>()->setNPart(npartG);
  getEntity<EventInitialState>()->setNColl(ncoll);

}

void parseParticle(const string &line) {
  stringstream stream(line);

  int par0{-1};      // status
  int pdgCode{0};    // id
  int par2{-1};      // nColl (can be negative in JAM2)

  float mass{0};
  float px{0};
  float py{0};
  float pz{0};
  float energy{0};
  float x{0};
  float y{0};
  float z{0};
  float time{0};

  // NOTE:
  // JAM2 can append extra trailing columns (e.g. tform).
  // We intentionally read only the first 12 columns and ignore the rest.
  stream
      >> par0 >> pdgCode >> par2
      >> mass >> px >> py >> pz >> energy
      >> x >> y >> z >> time;
  //if (std::fabs(par2) == 1) return;
  getEntity<UParticle>()->SetMate(par2);
  getEntity<UParticle>()->SetStatus(par0);
  getEntity<UParticle>()->SetPdg(pdgCode);
  getEntity<UParticle>()->SetPosition(x, y, z, time);
  getEntity<UParticle>()->SetMomentum(px, py, pz, energy);
}

void particlePostProcess() {
  getEntity<UEvent>()->AddParticle(*getEntity<UParticle>());
}

void previousEventPostProcess() {
  // getEntity<UEvent>()->Print();

  // keep the same selection logic as in the legacy reader:
  // store only events with non-zero NColl
  if (getEntity<EventInitialState>()->getNColl() != 0) {
    gEventTree->Fill();
  }
}

int main(int argc, char **argv) {
  std::cout << "go!" << std::endl;

  string inputFileName{"phase.dat"};
  if (argc > 1) {
    inputFileName = string(argv[1]);
  }

  string jamOutputFileName{"jamOutput.root"};
  if (argc > 2) {
    jamOutputFileName = string(argv[2]);
  }

  ifstream inputFile(inputFileName);
  if (!inputFile.is_open()) {
    Error(__func__, "Cannot open input file: %s", inputFileName.c_str());
    return 2;
  }

  gOutputFile = TFile::Open(jamOutputFileName.c_str(), "recreate");

  gEventTree = new TTree("events", "JAM events");
  gEventTree->Branch("event", getEntity<UEvent>());
  gEventTree->Branch("iniState", getEntity<EventInitialState>());

  // Create URun and fill with zeroes for now (JAM2 has no JAMRUN.DAT)
  initURunWithZeroes();

  bool hasPreviousEvent{false};
  unsigned int nEvents{0};

  string line;
  while (std::getline(inputFile, line)) {
    unsigned int nWords = countWordsInTheString(line);

    ELineType currentLineType{ELineType::kNotInitializedYet};

    switch (nWords) {
      // JAM2 has many comment/config lines starting with '#'.
      // We ignore them unless they match the event header word count.
      case 9:
        currentLineType = ELineType::kEventHeader;

        if (hasPreviousEvent) previousEventPostProcess();
        parseEventHeader(line);

        hasPreviousEvent = true;
        nEvents++;
        break;

      case 13:
        currentLineType = ELineType::kParticle;
        parseParticle(line);
        particlePostProcess();
        break;

      default:
        currentLineType = ELineType::kUnknown;
    }
  }

  if (hasPreviousEvent) {
    previousEventPostProcess();
  }

  // Write outputs
  gEventTree->Write();

  // Write URun with zeroes for now (placeholder)
  getEntity<URun>()->Write();

  gOutputFile->Close();

  return 0;
}
