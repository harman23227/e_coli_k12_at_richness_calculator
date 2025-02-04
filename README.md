This Python script analyzes the AT richness of the E. coli K-12 genome and compares it to the AT content of the DnaA box motifs. The script fetches the genome sequence from NCBI using the accession number U00096.3, locates the DnaA box motifs, and computes both the overall AT richness of the genome and the AT richness of the motifs.

>Downloads the E. coli K-12 genome from NCBI in GenBank format.
>Computes the overall AT richness (percentage of adenine and thymine bases) of the genome.
>Locates all instances of the DnaA box motif (TTATCCACA).
>Compares the AT richness of the DnaA box motifs to the rest of the genome.


Requirements:
Biopython: For handling genomic data (downloading and parsing the genome).
requests and Entrez: To fetch genome data from NCBI.
