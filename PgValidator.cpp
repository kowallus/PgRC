#include <cstdlib>
#include <unistd.h>

#include "pseudogenome/persistence/SeparatedPseudoGenomePersistence.h"
#include "pseudogenome/readslist/SeparatedExtendedReadsList.h"

using namespace std;
using namespace PgTools;

string
getRead(string pg, const ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &entry,
        uint_read_len_max readLength) {
    string read = pg.substr(entry.pos, readLength);
    for(uint8_t i = 0; i < entry.mismatchesCount; i++) {
        read[entry.mismatchOffset[i]] = code2mismatch(read[entry.mismatchOffset[i]], entry.mismatchCode[i]);
    }
    if (entry.revComp)
        read = reverseComplement(read);
    return read;
}

void
validatePg(const string &srcFastqFile, const string &pairFastqFile, const string &pgFilePrefix,
        uint_pg_len_max startPos) {
    clock_checkpoint();
    DefaultSeparatedExtendedReadsListIterator* rlIt = DefaultSeparatedExtendedReadsListIterator::getIterator(pgFilePrefix);
    PseudoGenomeHeader* pgh;
    bool plainTextReadMode;
    PgTools::SeparatedPseudoGenomePersistence::getPseudoGenomeProperties(pgFilePrefix, pgh, plainTextReadMode);
    uint_read_len_max readLength = pgh->getMaxReadLength();
    string pg = PgTools::SeparatedPseudoGenomePersistence::getPseudoGenome(pgFilePrefix);
    ReadsSourceIteratorTemplate<uint_read_len_max> *allReadsIterator = ReadsSetPersistence::createManagedReadsIterator(
            srcFastqFile, pairFastqFile);
    cout << "Reading reads set\n";
    PackedConstantLengthReadsSet *readsSet = new PackedConstantLengthReadsSet(allReadsIterator);
    readsSet->printout();
    delete (allReadsIterator);
    uint_reads_cnt_max pgReadsCounter = 0;
    uint_reads_cnt_max errorsCounter = 0;
    while(rlIt->moveNext()) {
        if (!(++pgReadsCounter % 100)) {
            cout << ".";
            if (!(pgReadsCounter % 10000)) cout << endl;
            cout.flush();
        }

        const ReadsListEntry<255, uint_read_len_max, uint_reads_cnt_max, uint_pg_len_max> &pgEntry = rlIt->peekReadEntry();
        if (pgEntry.pos < startPos)
            continue;
        string pgRead = getRead(pg, pgEntry, readLength);
        string rsRead = readsSet->getReadVirtual(pgEntry.idx);
        if (pgRead != rsRead) {
            errorsCounter++;
            if (errorsCounter < 10) {
                cout << endl << getRead(pg, pgEntry, readLength) << " (" << pgReadsCounter << ";" << pgEntry.pos << ")" << endl;
                cout << rsRead << " (" << pgEntry.idx << ";" << pgEntry.pos << ")" << endl << endl;
            }
        }
    }
    delete(readsSet);
    delete(pgh);
    delete(rlIt);
    cout << "Validated " << pgReadsCounter << " Pg entries, found " << errorsCounter << " errors." << endl;
    cout << "... validation completed in " << clock_millis() << " msec. " << endl;
}

int main(int argc, char *argv[])
{
    int opt; // current option
    uint_pg_len_max startPos = 0;

    while ((opt = getopt(argc, argv, "p:?")) != -1) {
        switch (opt) {
        case 'p':
            startPos = atoi(optarg);
            break;
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-p start_position] readssrcfile [pairsrcfile] pgfileprefix\n\n",
                    argv[0]);
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 or 3 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }

    string srcFastqFile(argv[optind++]);
    string pairFastqFile = "";
    if (optind == argc - 2)
        pairFastqFile = argv[optind++];
    string pgFilePrefix(argv[optind++]);

    validatePg(srcFastqFile, pairFastqFile, pgFilePrefix, startPos);

    exit(EXIT_SUCCESS);
}