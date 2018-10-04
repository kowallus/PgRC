#include "ReadsSetPersistence.h"

#include "../iterator/DivisionReadsSetDecorators.h"

namespace PgSAReadsSet {

    ReadsSourceIteratorTemplate<uint_read_len_max> *ReadsSetPersistence::createReadsIterator(string filename,
                                                                                        string pairfile) {
        istream* streamSource = new ifstream(filename, ios_base::in | ios_base::binary);
        istream* pairSource = 0;
        if (pairfile != "")
            pairSource = new ifstream(pairfile, ios_base::in | ios_base::binary);

        ReadsSourceIteratorTemplate<uint_read_len_max>* readsSource;
        if (filename.substr(filename.length() - 6) == ".fasta")
            readsSource = new FASTAReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
         else if (filename.substr(filename.length() - 6) == ".fastq")
            readsSource = new FASTQReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);
         else
            readsSource = new ConcatenatedReadsSourceIterator<uint_read_len_max>(streamSource);

        return readsSource;
    }

    FASTQReadsSourceIterator<uint_read_len_max> *ReadsSetPersistence::createFastQReadsIterator(string filename,
                                                                                             string pairfile) {
        istream* streamSource = new ifstream(filename, ios_base::in | ios_base::binary);
        istream* pairSource = 0;
        if (pairfile != "")
            pairSource = new ifstream(pairfile, ios_base::in | ios_base::binary);

        FASTQReadsSourceIterator<uint_read_len_max>* readsSource = 0;
        if (filename.substr(filename.length() - 6) == ".fastq")
            readsSource = new FASTQReadsSourceIterator<uint_read_len_max>(streamSource, pairSource);

        return readsSource;
    }

}