#include <TFile.h>
#include <TObject.h>
#include <TNamed.h>
#include <TClass.h>
#include <TDataMember.h>
#include <TString.h>
#include <TSystem.h>

#include <URun.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <stdexcept>

namespace
{
  constexpr double kNucleonMass = 0.9382720813; // GeV/c^2

  struct NucleusInfo {
    int A{0};
    int Z{0};
    std::string symbol;
  };

  struct JAM2RunConfig {
    int nEvents{0};

    NucleusInfo projectile;
    NucleusInfo target;

    double sqrtSNN{0.0};
    double bMin{0.0};
    double bMax{0.0};
  };

  std::string Trim(const std::string& s)
  {
    const auto begin = std::find_if_not(s.begin(), s.end(),
                                        [](unsigned char c) { return std::isspace(c); });

    const auto end = std::find_if_not(s.rbegin(), s.rend(),
                                      [](unsigned char c) { return std::isspace(c); }).base();

    if (begin >= end) return "";
    return std::string(begin, end);
  }

  std::string RemoveComment(const std::string& line)
  {
    const size_t hashPos = line.find('#');
    const size_t exclPos = line.find('!');

    size_t cutPos = std::string::npos;

    if (hashPos != std::string::npos) cutPos = hashPos;
    if (exclPos != std::string::npos) {
      if (cutPos == std::string::npos || exclPos < cutPos) {
        cutPos = exclPos;
      }
    }

    if (cutPos == std::string::npos) return line;
    return line.substr(0, cutPos);
  }

  bool StartsWith(const std::string& s, const std::string& prefix)
  {
    return s.rfind(prefix, 0) == 0;
  }

  std::map<std::string, int> BuildZMap()
  {
    return {
      {"H", 1},   {"He", 2},
      {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},   {"N", 7},   {"O", 8},
      {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14},
      {"P", 15},  {"S", 16},  {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20},
      {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24}, {"Mn", 25}, {"Fe", 26},
      {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32},
      {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38},
      {"Y", 39},  {"Zr", 40}, {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44},
      {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
      {"Sb", 51}, {"Te", 52}, {"I", 53},  {"Xe", 54}, {"Cs", 55}, {"Ba", 56},
      {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60}, {"Pm", 61}, {"Sm", 62},
      {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66}, {"Ho", 67}, {"Er", 68},
      {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74},
      {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
      {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85}, {"Rn", 86},
      {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90}, {"Pa", 91}, {"U", 92}
    };
  }

  NucleusInfo ParseNucleus(const std::string& input)
  {
    const std::string s = Trim(input);

    std::string aStr;
    std::string symbol;

    for (char c : s) {
      if (std::isdigit(static_cast<unsigned char>(c))) {
        aStr += c;
      } else if (std::isalpha(static_cast<unsigned char>(c))) {
        symbol += c;
      }
    }

    if (aStr.empty() || symbol.empty()) {
      throw std::runtime_error("Cannot parse nucleus name: \"" + input + "\"");
    }

    // Normalize symbol: xe -> Xe, XE -> Xe
    symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
    for (size_t i = 1; i < symbol.size(); ++i) {
      symbol[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[i])));
    }

    static const std::map<std::string, int> zMap = BuildZMap();

    const auto it = zMap.find(symbol);
    if (it == zMap.end()) {
      throw std::runtime_error("Unknown chemical symbol: \"" + symbol + "\" from \"" + input + "\"");
    }

    NucleusInfo n;
    n.A = std::stoi(aStr);
    n.Z = it->second;
    n.symbol = symbol;

    return n;
  }

  bool ParseKeyValueLine(const std::string& rawLine,
                         std::string& key,
                         std::string& value)
  {
    const std::string noComment = RemoveComment(rawLine);
    const size_t eqPos = noComment.find('=');

    if (eqPos == std::string::npos) {
      return false;
    }

    key = Trim(noComment.substr(0, eqPos));
    value = Trim(noComment.substr(eqPos + 1));

    return !key.empty() && !value.empty();
  }

  JAM2RunConfig ParseJAM2InputFile(const std::string& inputFileName)
  {
    std::ifstream input(inputFileName);

    if (!input.is_open()) {
      throw std::runtime_error("Cannot open JAM2 input file: " + inputFileName);
    }

    JAM2RunConfig cfg;

    bool foundNEvents = false;
    bool foundBeamA   = false;
    bool foundBeamB   = false;
    bool foundECM     = false;
    bool foundBMin    = false;
    bool foundBMax    = false;

    std::string line;

    while (std::getline(input, line)) {
      std::string key;
      std::string value;

      if (!ParseKeyValueLine(line, key, value)) {
        continue;
      }

      if (key == "Main:numberOfEvents") {
        cfg.nEvents = std::stoi(value);
        foundNEvents = true;
      } else if (key == "Beams:beamA") {
        cfg.projectile = ParseNucleus(value);
        foundBeamA = true;
      } else if (key == "Beams:beamB") {
        cfg.target = ParseNucleus(value);
        foundBeamB = true;
      } else if (key == "Beams:eCM") {
        cfg.sqrtSNN = std::stod(value);
        foundECM = true;
      } else if (key == "Beams:bmin") {
        cfg.bMin = std::stod(value);
        foundBMin = true;
      } else if (key == "Beams:bmax") {
        cfg.bMax = std::stod(value);
        foundBMax = true;
      }
    }

    if (!foundNEvents) throw std::runtime_error("Missing key: Main:numberOfEvents");
    if (!foundBeamA)   throw std::runtime_error("Missing key: Beams:beamA");
    if (!foundBeamB)   throw std::runtime_error("Missing key: Beams:beamB");
    if (!foundECM)     throw std::runtime_error("Missing key: Beams:eCM");
    if (!foundBMin)    throw std::runtime_error("Missing key: Beams:bmin");
    if (!foundBMax)    throw std::runtime_error("Missing key: Beams:bmax");

    return cfg;
  }

  double BeamMomentumCMPerNucleon(double sqrtSNN)
  {
    const double arg = 0.25 * sqrtSNN * sqrtSNN - kNucleonMass * kNucleonMass;
    return arg > 0.0 ? std::sqrt(arg) : 0.0;
  }

  TDataMember* GetMember(TObject* obj, const char* memberName)
  {
    if (!obj) return nullptr;

    TClass* cl = obj->IsA();
    if (!cl) return nullptr;

    TDataMember* dm = cl->GetDataMember(memberName);

    if (!dm) {
      std::cerr << "    Missing data member: " << memberName
                << " in class " << cl->GetName() << std::endl;
    }

    return dm;
  }

  bool SetIntMember(TObject* obj, const char* memberName, int value)
  {
    TDataMember* dm = GetMember(obj, memberName);
    if (!dm) return false;

    char* base = reinterpret_cast<char*>(obj);
    char* addr = base + dm->GetOffset();

    TString type = dm->GetTypeName();

    if (type == "Int_t" || type == "int") {
      *reinterpret_cast<Int_t*>(addr) = static_cast<Int_t>(value);
      return true;
    }

    if (type == "UInt_t" || type == "unsigned int") {
      *reinterpret_cast<UInt_t*>(addr) = static_cast<UInt_t>(value);
      return true;
    }

    if (type == "Long64_t" || type == "long long") {
      *reinterpret_cast<Long64_t*>(addr) = static_cast<Long64_t>(value);
      return true;
    }

    std::cerr << "    Unsupported integer member type for " << memberName
              << ": " << type << std::endl;

    return false;
  }

  bool SetDoubleMember(TObject* obj, const char* memberName, double value)
  {
    TDataMember* dm = GetMember(obj, memberName);
    if (!dm) return false;

    char* base = reinterpret_cast<char*>(obj);
    char* addr = base + dm->GetOffset();

    TString type = dm->GetTypeName();

    // Double32_t is stored as double in memory.
    if (type == "Double_t" || type == "double" || type == "Double32_t") {
      *reinterpret_cast<Double_t*>(addr) = static_cast<Double_t>(value);
      return true;
    }

    // Float16_t is stored as float in memory.
    if (type == "Float_t" || type == "float" || type == "Float16_t") {
      *reinterpret_cast<Float_t*>(addr) = static_cast<Float_t>(value);
      return true;
    }

    if (type == "Int_t" || type == "int") {
      *reinterpret_cast<Int_t*>(addr) = static_cast<Int_t>(value);
      return true;
    }

    if (type == "UInt_t" || type == "unsigned int") {
      *reinterpret_cast<UInt_t*>(addr) = static_cast<UInt_t>(value);
      return true;
    }

    std::cerr << "    Unsupported floating member type for " << memberName
              << ": " << type << std::endl;

    return false;
  }

  bool SetStringMember(TObject* obj, const char* memberName, const TString& value)
  {
    TDataMember* dm = GetMember(obj, memberName);
    if (!dm) return false;

    char* base = reinterpret_cast<char*>(obj);
    char* addr = base + dm->GetOffset();

    TString type = dm->GetTypeName();

    if (type == "TString") {
      *reinterpret_cast<TString*>(addr) = value;
      return true;
    }

    std::cerr << "    Unsupported string member type for " << memberName
              << ": " << type << std::endl;

    return false;
  }

  bool FillURunObject(TObject* run, const JAM2RunConfig& cfg)
  {
    if (!run) return false;

    std::cout << "  Read object class: " << run->ClassName() << std::endl;

    if (TString(run->ClassName()) != "URun") {
      std::cerr << "  ERROR: object named \"run\" was not read as URun.\n"
                << "  It was read as: " << run->ClassName() << "\n"
                << "  Check that libMcIniData.so / UniGen dictionary is loaded or linked."
                << std::endl;
      return false;
    }

    TNamed* runNamed = dynamic_cast<TNamed*>(run);

    if (!runNamed) {
      std::cerr << "  ERROR: URun object does not inherit from TNamed." << std::endl;
      return false;
    }

    runNamed->SetNameTitle("run", "Run Header");

    const double pCM = BeamMomentumCMPerNucleon(cfg.sqrtSNN);

    const TString generator = "JAM2";
    const TString comment   = "";
    const TString decayer   = "";

    // Following your current convention:
    //   fPProj = +p_cm
    //   fPTarg = -p_cm
    const double pProj =  pCM;
    const double pTarg = -pCM;

    // Following your current zeroed convention:
    const double bWeight = 0.0;
    const double phiMin  = 0.0;
    const double phiMax  = 0.0;
    const double sigma   = 0.0;

    bool ok = true;

    ok &= SetStringMember(run, "fGenerator", generator);
    ok &= SetStringMember(run, "fComment",   comment);
    ok &= SetStringMember(run, "fDecayer",   decayer);

    ok &= SetIntMember(run, "fAProj", cfg.projectile.A);
    ok &= SetIntMember(run, "fZProj", cfg.projectile.Z);

    ok &= SetIntMember(run, "fATarg", cfg.target.A);
    ok &= SetIntMember(run, "fZTarg", cfg.target.Z);

    ok &= SetDoubleMember(run, "fPProj", pProj);
    ok &= SetDoubleMember(run, "fPTarg", pTarg);

    ok &= SetDoubleMember(run, "fBMin",    cfg.bMin);
    ok &= SetDoubleMember(run, "fBMax",    cfg.bMax);
    ok &= SetDoubleMember(run, "fBWeight", bWeight);

    ok &= SetDoubleMember(run, "fPhiMin", phiMin);
    ok &= SetDoubleMember(run, "fPhiMax", phiMax);

    ok &= SetDoubleMember(run, "fSigma", sigma);

    ok &= SetIntMember(run, "fNEvents", cfg.nEvents);

    if (!ok) {
      std::cerr << "  ERROR: some URun fields were not updated correctly." << std::endl;
      return false;
    }

    std::cout << "  Filled URun values:" << std::endl;
    std::cout << "    AProj/ZProj = "
              << cfg.projectile.A << " / " << cfg.projectile.Z
              << " (" << cfg.projectile.symbol << ")" << std::endl;

    std::cout << "    ATarg/ZTarg = "
              << cfg.target.A << " / " << cfg.target.Z
              << " (" << cfg.target.symbol << ")" << std::endl;

    std::cout << "    sqrt(sNN)   = " << cfg.sqrtSNN << " GeV" << std::endl;
    std::cout << "    fPProj      = " << std::setprecision(9) << pProj << " GeV/c per nucleon" << std::endl;
    std::cout << "    fPTarg      = " << std::setprecision(9) << pTarg << " GeV/c per nucleon" << std::endl;
    std::cout << "    fBMin/fBMax = " << cfg.bMin << " / " << cfg.bMax << " fm" << std::endl;
    std::cout << "    fNEvents    = " << cfg.nEvents << std::endl;

    return true;
  }
}

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "\nUsage:\n"
              << "  " << argv[0] << " jam2_input.inp unigen.root\n\n"
              << "Example:\n"
              << "  " << argv[0] << " input.inp unigen.root\n"
              << std::endl;
    return 1;
  }

  const std::string jam2InputFileName = argv[1];
  const std::string unigenFileName    = argv[2];

  // Force dictionary reference when executable is linked with libMcIniData.
  // This helps ensure URun dictionary is registered before reading the file.
  if (!TClass::GetClass("URun")) {
    URun::Class();
  }

  if (!TClass::GetClass("URun")) {
    std::cerr << "ERROR: URun dictionary is not available.\n"
              << "Compile/link this executable with libMcIniData.so, or make sure the library is loaded."
              << std::endl;
    return 2;
  }

  JAM2RunConfig cfg;

  try {
    cfg = ParseJAM2InputFile(jam2InputFileName);
  } catch (const std::exception& e) {
    std::cerr << "ERROR while parsing JAM2 input file:\n  "
              << e.what() << std::endl;
    return 3;
  }

  std::cout << "\nParsed JAM2 input file:" << std::endl;
  std::cout << "  File        = " << jam2InputFileName << std::endl;
  std::cout << "  beamA       = " << cfg.projectile.A << cfg.projectile.symbol
            << "  Z=" << cfg.projectile.Z << std::endl;
  std::cout << "  beamB       = " << cfg.target.A << cfg.target.symbol
            << "  Z=" << cfg.target.Z << std::endl;
  std::cout << "  sqrt(sNN)   = " << cfg.sqrtSNN << " GeV" << std::endl;
  std::cout << "  bmin/bmax   = " << cfg.bMin << " / " << cfg.bMax << " fm" << std::endl;
  std::cout << "  nEvents     = " << cfg.nEvents << std::endl;

  TFile file(unigenFileName.c_str(), "UPDATE");

  if (file.IsZombie() || !file.IsOpen()) {
    std::cerr << "ERROR: cannot open UniGen file in UPDATE mode:\n  "
              << unigenFileName << std::endl;
    return 4;
  }

  TObject* run = file.Get("run");

  if (!run) {
    std::cerr << "ERROR: object named \"run\" was not found in file:\n  "
              << unigenFileName << std::endl;
    file.Close();
    return 5;
  }

  if (!FillURunObject(run, cfg)) {
    file.Close();
    return 6;
  }

  file.cd();

  const Int_t bytes = run->Write("run", TObject::kOverwrite);

  if (bytes <= 0) {
    std::cerr << "ERROR: failed to write updated URun object." << std::endl;
    file.Close();
    return 7;
  }

  file.Close();

  std::cout << "\nSuccessfully updated URun in file:\n  "
            << unigenFileName << std::endl;

  return 0;
}