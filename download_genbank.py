from Bio import Entrez, SeqIO
import time
import os

# Configuration
Entrez.email = "zeliz22@freeuni.edu.ge"
input_file = "data/LQ/fetch_targets.txt"  # Path to your LQ results
output_directory = "data/genbank_files"

def load_species_from_file(filepath):
    """Reads species list from a text file, one per line."""
    if not os.path.exists(filepath):
        print(f"Error: Target file {filepath} not found!")
        return []
    
    with open(filepath, 'r') as f:
        # Read lines, strip whitespace, and replace underscores with spaces for NCBI
        species = [line.strip().replace('_', ' ') for line in f if line.strip()]
    return species

def find_mitochondrial_accession(species_name):
    """Search NCBI for mitochondrial genome RefSeq accession number."""
    print(f"  Searching for accession number for {species_name}...")
    try:
        # search_term targets the RefSeq 'NC_' entries which are curated complete genomes
        search_term = f'"{species_name}"[Organism] AND mitochondrion[Title] AND complete genome AND RefSeq[Filter]'
        
        handle = Entrez.esearch(db="nuccore", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            print(f"  No RefSeq found, trying broader search...")
            search_term = f'"{species_name}"[Organism] AND mitochondrion[Title] AND complete genome'
            handle = Entrez.esearch(db="nuccore", term=search_term, retmax=5)
            record = Entrez.read(handle)
            handle.close()

        if record['IdList']:
            uid = record['IdList'][0]
            handle = Entrez.efetch(db="nuccore", id=uid, rettype="acc", retmode="text")
            accession = handle.read().strip()
            handle.close()
            print(f"  ✓ Found accession: {accession}")
            return accession
        else:
            print(f"  ✗ No mitochondrial genome found")
            return None
    except Exception as e:
        print(f"  ✗ Search error: {e}")
        return None

def download_genbank(species_name, accession, output_dir):
    """Download GenBank file and save to local directory."""
    try:
        print(f"  Downloading GenBank file...")
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Save with underscores for filename consistency
        species_filename = species_name.replace(' ', '_')
        output_file = os.path.join(output_dir, f"{species_filename}.gb")

        SeqIO.write(record, output_file, "genbank")
        print(f"  ✓ Saved to {output_file} ({len(record.seq):,} bp)")
        
        cds_count = sum(1 for f in record.features if f.type == "CDS")
        print(f"  CDS features: {cds_count}")
        return True
    except Exception as e:
        print(f"  ✗ Download error: {e}")
        return False

def main():
    # 1. Load the targets from your file
    species_list = load_species_from_file(input_file)
    
    if not species_list:
        print("No species found to process. Check your fetch_targets.txt file.")
        return

    os.makedirs(output_directory, exist_ok=True)

    print("=" * 70)
    print("AUTOMATED MITOCHONDRIAL GENOME DOWNLOADER")
    print("=" * 70)
    print(f"Reading from: {input_file}")
    print(f"Number of species: {len(species_list)}")
    print("=" * 70)

    successful, failed = [], []

    for i, species_name in enumerate(species_list, 1):
        print(f"\n[{i}/{len(species_list)}] Processing: {species_name}")
        
        accession = find_mitochondrial_accession(species_name)
        if accession:
            if download_genbank(species_name, accession, output_directory):
                successful.append(species_name)
            else:
                failed.append(species_name)
        else:
            failed.append(species_name)

        # Respect NCBI: Do not remove this sleep timer!
        time.sleep(0.5)

    print("\n" + "=" * 70)
    print(f"✓ DONE: {len(successful)} downloaded, {len(failed)} failed.")
    print(f"Files are in: {output_directory}")
    print("=" * 70)

if __name__ == "__main__":
    main()