#include <string>
#include <cstdlib>
#include <unistd.h>


using namespace std;

void matchReadsInPg(const string &pgFile, const string &outPrefix) {

}

int main(int argc, char *argv[])
{

    int opt; // current option

    while ((opt = getopt(argc, argv, "?")) != -1) {
        switch (opt) {
        case '?':
        default: /* '?' */
            fprintf(stderr, "Usage: %s pgfile outputprefix\n\n",
                    argv[0]);
            fprintf(stderr, "\n\n");
            exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2) || optind < (argc - 3)) {
        fprintf(stderr, "%s: Expected 2 arguments (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        
        exit(EXIT_FAILURE);
    }
        
    string pgFile(argv[optind++]);
    string outPrefix(argv[optind++]);

    matchReadsInPg(pgFile, outPrefix);
   
    exit(EXIT_SUCCESS);
}