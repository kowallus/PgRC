# PgRC: Pseudogenome based Read Compressor
 
Pseudogenome-based Read Compressor (PgRC) is an in-memory algorithm 
for compressing the DNA stream of FASTQ datasets, based on the idea 
of building an approximation of the shortest common superstring over 
high-quality reads.

The current implementation supports constant-length reads limited
to 255 bases.

### Installation on Linux
The following steps create an PgRC executable. 
On Linux PgRC build requires installed cmake version >= 3.4 (check using ```cmake --version```):
```bash
git clone https://github.com/kowallus/PgRC.git
cd PgRC
mkdir build
cd build
cmake ..
make PgRC
```

### Basic usage

```
PgRC [-c compressionLevel] [-i seqSrcFile [pairSrcFile]] [-t noOfThreads] [-o] [-d] archiveName
   
   -c compression levels: 1 - fast; 2 - default; 3 - max
   -t number of threads used (8 - default)
   -d decompression mode
   -o preserve original read order information
```

compression of DNA stream in order non-preserving regime (SE mode):
```
./PgRC -i in.fastq comp.pgrc
```
compression of DNA stream in order preserving regime (SE_ORD mode):
```
./PgRC -o -i in.fastq comp.pgrc
```
compression of paired-end DNA stream in order non-preserving regime (PE mode):
```
./PgRC -i in1.fastq in2.fastq comp.pgrc
```
compression of paired-end DNA stream in order preserving regime (PE mode):
```
./PgRC -o -i in.fastq comp.pgrc
```

## Publications

[Tomasz M. Kowalski, Szymon Grabowski: PgRC: pseudogenome-based read compressor. Bioinformatics, Volume 36, Issue 7, pp. 2082–2089 (2020).](https://academic.oup.com/bioinformatics/article/36/7/2082/5670526)

[supplementary data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/36/7/10.1093_bioinformatics_btz919/1/btz919_supplementary_data.pdf?Expires=1589991362&Signature=AHL3zC8GJMsQxG4FYiWFS8cRPi5x~6iByqNvgr0rBkqs3KuVtr42-GASV0fBdzY5SQGGIvT4tB5QPm6cGjV9plKKBfhC5QMMQlQpTJjcYUfEnELSM9IKhWh0qw6Px4gsRuArzJYJ0zxQBhiHi8yw~vKQ68czbO7VxKl5jwC2TCjszX~0FrOI1WFKJpMHOAF0kHZb9O45i2WwQHkx6ZAgedGWYLk6DOi0KRYvNcRjgOH-q94TcEpHWdERburrrLt0mCpda~E6jW7xWVew7ymwZAM5W7wtPK5UrUotwEc9h1jY2DuYdcaxF4Wd4nxadKD2tGh4Nc9rdKNjtMQQZRIHXQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

[bioRxiv](https://www.biorxiv.org/content/10.1101/710822v1)

## Related projects
[PgSA - <i>Pseudogenome Suffix Array</i>](https://github.com/kowallus/PgSA)
