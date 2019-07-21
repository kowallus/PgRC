# PgRC: Pseudogenome based Read Compressor
 
Pseudogenome-based Read Compressor (PgRC) is an in-memory algorithm for compressing the DNA stream of FASTQ datasets, based on the idea of building an approximation of the shortest common superstring over high-quality reads.

### Installation on Linux
The following steps create PgRC executable. 
On Linux PgRC build requires installed cmake version >= 3.4 (check using ```cmake --version```):
```bash
git clone https://github.com/kowallus/PgRC.git
cd PgRC
mkdir build
cd build
cmake ..
make PgRC
```

## Related projects
PgSA - <i>Pseudogenome Suffix Array</i> (https://github.com/kowallus/PgSA)
