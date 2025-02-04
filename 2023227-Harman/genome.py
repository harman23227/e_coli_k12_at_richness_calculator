from Bio import Entrez, SeqIO


Entrez.email = "harman23227@iiitd.ac.in"  


def download_genome():
    accession_number = "U00096.3"  
    print(f"Fetching E. coli K-12 genome with accession number: {accession_number}")


    handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
    genome_record = SeqIO.read(handle, "genbank")
    handle.close()


    with open("e_coli_k12_genome.fasta", "w") as f:
        SeqIO.write(genome_record, f, "fasta")
    print("Genome saved to e_coli_k12_genome.fasta")

    return genome_record


def calculate_at_richness(sequence):
    at_count = sequence.count("A") + sequence.count("T")
    total_length = len(sequence)
    at_richness = (at_count / total_length) * 100
    return at_richness


def find_dnaa_motifs(sequence):
    dnaa_box = "TTATCCACA"
    motif_positions = [i for i in range(len(sequence) - len(dnaa_box) + 1) if sequence[i:i+len(dnaa_box)] == dnaa_box]
    return motif_positions


def compare_at_richness(sequence, motif_positions):

    dnaa_box = "TTATCCACA"
    motif_at_count = dnaa_box.count("A") + dnaa_box.count("T")
    motif_at_richness = (motif_at_count / len(dnaa_box)) * 100


    total_at_count = sequence.count("A") + sequence.count("T")
    total_length = len(sequence)
    genome_at_richness = (total_at_count / total_length) * 100

    return motif_at_richness, genome_at_richness


def main():

    genome_record = download_genome()
    genome_sequence = genome_record.seq


    at_richness = calculate_at_richness(genome_sequence)
    print(f"AT Richness of the genome: {at_richness:.2f}%")


    motif_positions = find_dnaa_motifs(genome_sequence)
    print(f"Number of DnaA box motifs found: {len(motif_positions)}")
    print(f"Positions of DnaA box motifs: {motif_positions}")


    motif_at_richness, genome_at_richness = compare_at_richness(genome_sequence, motif_positions)
    difference = motif_at_richness - genome_at_richness
    print(f"AT Richness of DnaA box motifs: {motif_at_richness:.2f}%")
    print(f"AT Richness of the rest of the genome: {genome_at_richness:.2f}%")
    print(f"Difference in AT Richness between genome , DnaA box motif :{difference:.2f}%")


    print("\nHypothesis about the DnaA box motif:")
    print("The hypothesis proposes that the DnaA protein preferentially binds to regions of the genome with higher AT content, particularly at the DnaA box motifs. This selective binding is essential for the initiation of DNA replication at the oriC (origin of replication) region in E. coli. By analyzing the AT content, we predict that the DnaA box motifs should have a higher AT content than the surrounding genomic regions, which would support the theory that these motifs are specialized for DnaA protein binding and initiation of replication.")
   

if __name__ == "__main__":
    main()
