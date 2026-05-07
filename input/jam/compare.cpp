    // CompareUnigenWithDat_JAM2.C
    //
    // Usage (ROOT):
    //   .L CompareUnigenWithDat_JAM2.C+
    //   CompareUnigenWithDat_JAM2("unigen.root","phase.dat", true, 1e-6, 1e-9, -1, false, true);
    //
    // Arguments:
    //   rootFile         - produced unigen.root (tree "events")
    //   datFile          - original JAM2 .dat file
    //   skipZeroCollDat  - skip .dat events with ncoll==0 (set true if your reader only stored ncoll!=0)
    //   relTol, absTol   - tolerances for floating comparisons
    //   maxEvents        - limit number of ROOT entries to check (-1 = all)
    //   stopOnFirst      - stop immediately on first mismatch
    //   checkMassFromEP  - also compare .dat mass to sqrt(E^2 - p^2) from ROOT values

    #include <iostream>
    #include <fstream>
    #include <sstream>
    #include <vector>
    #include <string>
    #include <cmath>
    #include <algorithm>
    #include <TSystem.h>

    #include <TFile.h>
    #include <TTree.h>
    #include <TTreeReader.h>
    #include <TTreeReaderValue.h>
    #include <TTreeReaderArray.h>

    namespace {
    using std::string;
    using std::ifstream;
    using std::stringstream;
    using std::vector;

    struct DatParticle {
    int    status{-999};
    int    pdg{0};
    int    nColl{-999};     // present in .dat, NOT available in root tree you showed
    double mass{0};

    double px{0}, py{0}, pz{0}, E{0};
    double x{0}, y{0}, z{0}, t{0};

    bool   hasTform{false};
    double tform{0};

    long   lineNo{-1};
    string lineText;
    };

    struct DatEvent {
    unsigned int eventId{0};
    unsigned int nv{0};

    unsigned int ncollG{0};
    unsigned int npartG{0};

    double b{0};

    unsigned int npart{0};
    unsigned int ncoll{0};
    unsigned int ncollBB{0};

    long   headerLineNo{-1};
    string headerLineText;

    vector<DatParticle> particles;
    };

    unsigned int countWordsInTheString(const std::string &str) {
    std::stringstream stream(str);
    std::string word;
    unsigned int words = 0;
    while (stream >> word) ++words;
    return words;
    }

    bool startsWithHash(const string &line) {
    auto p = line.find_first_not_of(" \t\r\n");
    return (p != string::npos && line[p] == '#');
    }

    bool approxEqual(double a, double b, double relTol, double absTol) {
    const double diff = std::fabs(a - b);
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    return diff <= std::max(absTol, relTol * scale);
    }

    // Reads next JAM2 event from .dat (including its nv particle lines).
    // Skips unrelated '#' lines (config/comments) automatically.
    // Returns false on EOF.
    bool readNextJam2Event(ifstream &fin, DatEvent &ev, long &lineCounter) {
    ev = DatEvent{};

    string line;
    while (std::getline(fin, line)) {
        ++lineCounter;
        if (countWordsInTheString(line) == 0) continue;

        // We only accept JAM2 event headers of the form:
        // # iev nv ncollG npartG b npart ncoll ncollBB
        if (!startsWithHash(line)) continue;

        // Many JAM2 files have lots of '#' comment/config lines; try parsing, if fail -> continue.
        stringstream ss(line);
        string hashStr;
        ss >> hashStr;
        if (hashStr != "#") continue;

        unsigned int iev=0, nv=0, ncollG=0, npartG=0, npart=0, ncoll=0, ncollBB=0;
        double b=0;

        if (!(ss >> iev >> nv >> ncollG >> npartG >> b >> npart >> ncoll >> ncollBB)) {
        continue; // not an event header
        }

        ev.eventId = iev;
        ev.nv = nv;
        ev.ncollG = ncollG;
        ev.npartG = npartG;
        ev.b = b;
        ev.npart = npart;
        ev.ncoll = ncoll;
        ev.ncollBB = ncollBB;

        ev.headerLineNo = lineCounter;
        ev.headerLineText = line;

        // Now read nv particle lines
        ev.particles.clear();
        ev.particles.reserve(ev.nv);

        for (unsigned int i = 0; i < ev.nv; ++i) {
        string pline;
        if (!std::getline(fin, pline)) {
            std::cerr << "[ERROR] EOF while reading particles for event " << ev.eventId
                    << " (expected nv=" << ev.nv << ").\n";
            return false;
        }
        ++lineCounter;

        // Allow blank lines (rare), keep reading until we get a non-empty line
        while (countWordsInTheString(pline) == 0) {
            if (!std::getline(fin, pline)) {
            std::cerr << "[ERROR] EOF while skipping empty lines in particles for event "
                        << ev.eventId << ".\n";
            return false;
            }
            ++lineCounter;
        }

        // Particle lines should NOT start with '#'
        if (startsWithHash(pline)) {
            std::cerr << "[ERROR] Unexpected header/comment line inside particle block at line "
                    << lineCounter << " for event " << ev.eventId << ":\n"
                    << "  " << pline << "\n";
            return false;
        }

        DatParticle p;
        p.lineNo = lineCounter;
        p.lineText = pline;

        // JAM2 particle output (from your fill()):
        // status  id  nColl  mass  px py pz pe  x y z t  tform
        // tform may be present (we treat as optional)
        stringstream ps(pline);
        if (!(ps >> p.status >> p.pdg >> p.nColl
                    >> p.mass
                    >> p.px >> p.py >> p.pz >> p.E
                    >> p.x >> p.y >> p.z >> p.t)) {
            std::cerr << "[ERROR] Failed to parse particle line " << lineCounter
                    << " (event " << ev.eventId << "):\n"
                    << "  " << pline << "\n";
            return false;
        }
        double extra = 0;
        if (ps >> extra) {
            p.hasTform = true;
            p.tform = extra;
        }

        ev.particles.push_back(p);
        }

        return true;
    }

    return false; // EOF
    }

    // Reads .dat until it finds the next event that should correspond to a ROOT entry.
    // If skipZeroCollDat==true, it will skip events with ncoll==0 (consuming their particle blocks).
    bool readNextComparableDatEvent(ifstream &fin, DatEvent &ev, long &lineCounter, bool skipZeroCollDat) {
    while (true) {
        const bool ok = readNextJam2Event(fin, ev, lineCounter);
        if (!ok) return false;
        if (skipZeroCollDat && ev.ncoll == 0) {
        // skip this event and continue
        continue;
        }
        return true;
    }
    }

    void printEventContext(const char* tag,
                        Long64_t iEntry,
                        int rootEventNr,
                        const DatEvent &dev) {
    std::cout << tag << " ROOT entry=" << iEntry
                << " rootEventNr=" << rootEventNr
                << " datEventId=" << dev.eventId
                << " datHeaderLine=" << dev.headerLineNo << "\n"
                << "  DAT header: " << dev.headerLineText << "\n";
    }

    } // namespace


    void CompareUnigenWithDat_JAM2(const char *rootFile="unigen_JAM2.root",
                                const char *datFile ="phase_10events.dat",
                                bool  skipZeroCollDat=true,
                                double relTol=1e-2,
                                double absTol=1e-2,
                                int   maxEvents=-1,
                                bool  stopOnFirst=false,
                                bool  checkMassFromEP=true)
    {
    gSystem->Load("/lhep/users/bikmenik/mcini/build/libMcIniData.so");
    // -------------------- open ROOT --------------------
    TFile *f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Cannot open ROOT file: " << rootFile << "\n";
        return;
    }

    TTree *t = (TTree*)f->Get("events");
    if (!t) {
        std::cerr << "[ERROR] Tree 'events' not found in: " << rootFile << "\n";
        f->Close();
        return;
    }

    // -------------------- ROOT reader --------------------
    TTreeReader r(t);

    // Event header-like fields
    TTreeReaderValue<Int_t>    r_evtNr(r, "event.fEventNr");
    TTreeReaderValue<Double_t> r_b    (r, "event.fB");
    //   TTreeReaderValue<Int_t>    r_nPartArr(r, "event.fParticles"); // number of particles in array

    // iniState fields
    TTreeReaderValue<UShort_t> r_iniId (r, "iniState.id");
    TTreeReaderValue<UInt_t>   r_iniNColl(r, "iniState.nColl");
    TTreeReaderValue<UInt_t>   r_iniNPart(r, "iniState.nPart");

    // Particle arrays
    TTreeReaderArray<Int_t>      r_pdg   (r, "event.fParticles.fPdg");
    TTreeReaderArray<Int_t>      r_stat  (r, "event.fParticles.fStatus");
    TTreeReaderArray<Double_t>   r_px    (r, "event.fParticles.fPx");
    TTreeReaderArray<Double_t>   r_py    (r, "event.fParticles.fPy");
    TTreeReaderArray<Double_t>   r_pz    (r, "event.fParticles.fPz");
    TTreeReaderArray<Double_t>   r_E     (r, "event.fParticles.fE");
    TTreeReaderArray<Double_t>   r_x     (r, "event.fParticles.fX");
    TTreeReaderArray<Double_t>   r_y     (r, "event.fParticles.fY");
    TTreeReaderArray<Double_t>   r_z     (r, "event.fParticles.fZ");
    TTreeReaderArray<Double_t>   r_t     (r, "event.fParticles.fT");


    // -------------------- open DAT --------------------
    ifstream fin(datFile);
    if (!fin.is_open()) {
        std::cerr << "[ERROR] Cannot open DAT file: " << datFile << "\n";
        f->Close();
        return;
    }
    long datLineCounter = 0;

    // -------------------- loop and compare --------------------
    Long64_t iEntry = 0;
    int nMismatchEvents = 0;
    int nMismatchParticles = 0;

    DatEvent dev;

    while (r.Next()) {
        if (maxEvents >= 0 && (int)iEntry >= maxEvents) break;

        // read corresponding dat event
        if (!readNextComparableDatEvent(fin, dev, datLineCounter, skipZeroCollDat)) {
        std::cerr << "[ERROR] Reached EOF in DAT before consuming all ROOT entries.\n";
        std::cerr << "        ROOT entry=" << iEntry << " rootEventNr=" << *r_evtNr << "\n";
        break;
        }

        bool eventMismatch = false;

        // -------- compare event id --------
        if ((unsigned int)(*r_evtNr) != dev.eventId) {
        printEventContext("[MISMATCH eventId]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT event.fEventNr=" << *r_evtNr
                    << "  DAT eventId=" << dev.eventId << "\n";
        eventMismatch = true;
        }

        // -------- compare impact parameter b --------
        if (!approxEqual(*r_b, dev.b, relTol, absTol)) {
        printEventContext("[MISMATCH b]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT b=" << *r_b << "  DAT b=" << dev.b
                    << "  (relTol=" << relTol << ", absTol=" << absTol << ")\n";
        eventMismatch = true;
        }

        // -------- compare iniState --------
        if ((unsigned int)(*r_iniId) != dev.eventId) {
        printEventContext("[MISMATCH iniState.id]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT iniState.id=" << (unsigned int)(*r_iniId)
                    << "  DAT eventId=" << dev.eventId << "\n";
        eventMismatch = true;
        }

        if ((unsigned int)(*r_iniNPart) != dev.npart) {
        printEventContext("[MISMATCH nPart]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT iniState.nPart=" << (unsigned int)(*r_iniNPart)
                    << "  DAT npart=" << dev.npart << "\n";
        eventMismatch = true;
        }

        if ((unsigned int)(*r_iniNColl) != dev.ncoll) {
        printEventContext("[MISMATCH nColl]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT iniState.nColl=" << (unsigned int)(*r_iniNColl)
                    << "  DAT ncoll=" << dev.ncoll << "\n";
        eventMismatch = true;
        }

        // -------- compare particle count --------
        const int rootNv = (int)r_pdg.GetSize();
        if ((unsigned int)rootNv != dev.nv) {
        printEventContext("[MISMATCH nv]", iEntry, *r_evtNr, dev);
        std::cout << "  ROOT nv(from arrays)=" << rootNv
                    << "  DAT nv=" << dev.nv << "\n";
        eventMismatch = true;
        }

        // -------- compare particles --------
        const int nCompare = std::min<int>(rootNv, (int)dev.particles.size());

        for (int i = 0; i < nCompare; ++i) {
        const DatParticle &dp = dev.particles[i];

        bool pm = false;

        if (r_pdg[i] != dp.pdg) pm = true;
        if (!approxEqual(r_px[i], dp.px, relTol, absTol)) pm = true;
        if (!approxEqual(r_py[i], dp.py, relTol, absTol)) pm = true;
        if (!approxEqual(r_pz[i], dp.pz, relTol, absTol)) pm = true;
        if (!approxEqual(r_E [i], dp.E,  relTol, absTol)) pm = true;
        if (!approxEqual(r_x [i], dp.x,  relTol, absTol)) pm = true;
        if (!approxEqual(r_y [i], dp.y,  relTol, absTol)) pm = true;
        if (!approxEqual(r_z [i], dp.z,  relTol, absTol)) pm = true;
        if (!approxEqual(r_t [i], dp.t,  relTol, absTol)) pm = true;

        // Optional: check mass via E^2 - p^2 (helpful because mass is not stored in the ROOT branches)
        if (checkMassFromEP) {
            const double p2 = r_px[i]*r_px[i] + r_py[i]*r_py[i] + r_pz[i]*r_pz[i];
            const double m2 = r_E[i]*r_E[i] - p2;
            const double m  = (m2 > 0 ? std::sqrt(m2) : 0.0);
            if (!approxEqual(m, dp.mass, relTol*5.0, absTol*5.0)) {
            pm = true;
            }
        }

        if (pm) {
            ++nMismatchParticles;
            if (!eventMismatch) {
            // Print event header context once
            printEventContext("[MISMATCH particles]", iEntry, *r_evtNr, dev);
            eventMismatch = true;
            }

            std::cout << "  Particle i=" << i
                    << "  DAT line=" << dp.lineNo << "\n"
                    << "    DAT: " << dp.lineText << "\n"
                    << "    ROOT: pdg=" << r_pdg[i]
                    << " px=" << r_px[i] << " py=" << r_py[i] << " pz=" << r_pz[i] << " E=" << r_E[i]
                    << " x=" << r_x[i]  << " y=" << r_y[i]  << " z=" << r_z[i]  << " t=" << r_t[i]
                    << "\n";

            // If you want to see the specific numeric diffs:
            std::cout << "    DIFF: dpx=" << (r_px[i]-dp.px)
                    << " dpy=" << (r_py[i]-dp.py)
                    << " dpz=" << (r_pz[i]-dp.pz)
                    << " dE="  << (r_E[i]-dp.E)
                    << " dx="  << (r_x[i]-dp.x)
                    << " dy="  << (r_y[i]-dp.y)
                    << " dz="  << (r_z[i]-dp.z)
                    << " dt="  << (r_t[i]-dp.t)
                    << "\n";

            if (stopOnFirst) break;
        }
        }

        if (eventMismatch) {
        ++nMismatchEvents;
        if (stopOnFirst) break;
        }

        ++iEntry;
    }

    // If DAT has more comparable events than ROOT entries, warn (optional)
    // (We don't consume to EOF to avoid long scans on big files.)

    std::cout << "\n=== Comparison summary ===\n";
    std::cout << "Checked ROOT entries: " << iEntry << "\n";
    std::cout << "Events with mismatches: " << nMismatchEvents << "\n";
    std::cout << "Particle mismatches: " << nMismatchParticles << "\n";
    std::cout << "skipZeroCollDat=" << (skipZeroCollDat ? "true" : "false")
                << "  relTol=" << relTol << "  absTol=" << absTol
                << "  massCheck=" << (checkMassFromEP ? "on" : "off") << "\n";

    fin.close();
    f->Close();
    }
