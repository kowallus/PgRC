#include <cstdlib>
#include <unistd.h>

#include "matcher.h"

using namespace std;

static const string PSEUDOGENOME_HEADER = "PGEN";
static const string OFFSETS_SUFFIX = "_matched_offsets.txt";
static const string MISSED_READS_SUFFIX = "_missed.txt";

static const int PSEUDOGENOME_PRECEDING_LINES = 13;

string getPgFromFile(ifstream &pgSrc) {
    string line;
    pgSrc >> line;
    if (line == PSEUDOGENOME_HEADER) {
        for (int i = 0; i < PSEUDOGENOME_PRECEDING_LINES; i++)
            pgSrc >> line;
    }
    return line;
}

void matchReadsInPgFile(const string &pgFile, const string &readsFile, const string &outPrefix) {
    std::ifstream pgSrc(pgFile, std::ios::in | std::ios::binary);
    if (pgSrc.fail()) {
        fprintf(stderr, "cannot open pseudogenome file %s\n", pgFile.c_str());
        exit(EXIT_FAILURE);
    }
    std::ifstream readsSrc(readsFile, std::ios::in | std::ios::binary);
    if (readsSrc.fail()) {
        fprintf(stderr, "cannot open readsfile %s\n", readsFile.c_str());
        exit(EXIT_FAILURE);
    }
    string pg = getPgFromFile(pgSrc);
    pgSrc.close();
    string offsetsFile = outPrefix + OFFSETS_SUFFIX;
    std::ofstream offsetsDest(offsetsFile, std::ios::out | std::ios::binary);
    if (offsetsDest.fail()) {
        fprintf(stderr, "cannot write to offsets file %s\n", offsetsFile.c_str());
        exit(EXIT_FAILURE);
    }
    string missedReadsFile = outPrefix + MISSED_READS_SUFFIX;
    std::ofstream missedReadsDest(missedReadsFile, std::ios::out | std::ios::binary);
    if (missedReadsDest.fail()) {
        fprintf(stderr, "cannot write to missed reads file %s\n", missedReadsFile.c_str());
        exit(EXIT_FAILURE);
    }

    exactMatchConstantLengthPatterns(pg, readsSrc, offsetsDest, missedReadsDest);

    offsetsDest.close();
    missedReadsDest.close();
    readsSrc.close();
}

int main(int argc, char *argv[])
{

    int opt; // current option

    while ((opt = getopt(argc, argv, "?")) != -1) {
        switch (opt) {
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s pgfile readsfile outputprefix\n\n",
                    argv[0]);
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 3) || optind < (argc - 4)) {
        fprintf(stderr, "%s: Expected 3 arguments (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string pgFile(argv[optind++]);
    string readsFile(argv[optind++]);
    string outPrefix(argv[optind++]);

    matchReadsInPgFile(pgFile, readsFile, outPrefix);
   
    exit(EXIT_SUCCESS);
}