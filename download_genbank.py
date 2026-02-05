from Bio import Entrez, SeqIO
import time
import os


Entrez.email = "zeliz22@freeuni.edu.ge"  # CHANGE THIS to your email!

species_list = [
    "Mus musculus",
    "Myotis lucifugus",
    "Loxodonta africana",
    "Mesocricetus auratus",
    "Orcinus orca",
    "Didelphis virginiana",
    "Balaena mysticetus",
    "Rattus norvegicus",
    # Add more species here as needed
    # "Homo sapiens",
    # "Heterocephalus glaber",
    # etc...
]

output_directory = "data/genbank_files"



def find_mitochondrial_accession(species_name):
    """
    Search NCBI for mitochondrial genome RefSeq accession number

    Args:
        species_name: Scientific name of species (e.g., "Homo sapiens")

    Returns:
        Accession number (e.g., "NC_012920") or None if not found
    """
    print(f"  Searching for accession number...")

    try:
        # Search for mitochondrial genome
        search_term = f'"{species_name}"[Organism] AND mitochondrion[Title] AND complete genome AND RefSeq[Filter]'

        handle = Entrez.esearch(db="nuccore", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            # Try alternative search without RefSeq filter
            print(f"  No RefSeq found, trying broader search...")
            search_term = f'"{species_name}"[Organism] AND mitochondrion[Title] AND complete genome'

            handle = Entrez.esearch(db="nuccore", term=search_term, retmax=5)
            record = Entrez.read(handle)
            handle.close()

        if record['IdList']:
            # Get the first result ID
            uid = record['IdList'][0]

            # Fetch summary to get accession number
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
    """
    Download GenBank file for a species

    Args:
        species_name: Species name
        accession: RefSeq accession number
        output_dir: Directory to save file

    Returns:
        True if successful, False otherwise
    """
    try:
        # Download GenBank format
        print(f"  Downloading GenBank file...")
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Create filename from species name
        species_filename = species_name.replace(' ', '_').replace('.', '')
        output_file = os.path.join(output_dir, f"{species_filename}.gb")

        # Save to file
        SeqIO.write(record, output_file, "genbank")

        print(f"  ✓ Saved to {output_file}")
        print(f"  Length: {len(record.seq):,} bp")

        # Count CDS features (protein-coding genes)
        cds_count = sum(1 for f in record.features if f.type == "CDS")
        print(f"  CDS features: {cds_count}")

        return True

    except Exception as e:
        print(f"  ✗ Download error: {e}")
        return False


def main():
    """Main function to process all species"""

    # Create output directory
    os.makedirs(output_directory, exist_ok=True)

    print("=" * 70)
    print("AUTOMATED MITOCHONDRIAL GENOME DOWNLOADER")
    print("=" * 70)
    print(f"Number of species to process: {len(species_list)}")
    print(f"Output directory: {output_directory}")
    print("=" * 70)

    # Track results
    successful = []
    failed = []
    accession_dict = {}

    # Process each species
    for i, species_name in enumerate(species_list, 1):
        print(f"\n[{i}/{len(species_list)}] Processing: {species_name}")
        print("-" * 70)

        # Step 1: Find accession number
        accession = find_mitochondrial_accession(species_name)

        if accession is None:
            failed.append(species_name)
            print(f"  ⚠️  Skipping {species_name} - no accession found")
            time.sleep(0.5)  # Be nice to NCBI
            continue

        # Store accession for later reference
        species_key = species_name.replace(' ', '_')
        accession_dict[species_key] = accession

        # Step 2: Download GenBank file
        success = download_genbank(species_name, accession, output_directory)

        if success:
            successful.append(species_name)
        else:
            failed.append(species_name)

        # Be nice to NCBI servers - wait between requests
        time.sleep(0.5)

    # Print summary
    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)
    print(f"✓ Successful: {len(successful)}/{len(species_list)}")
    print(f"✗ Failed: {len(failed)}/{len(species_list)}")

    if successful:
        print("\nSuccessfully downloaded:")
        for species in successful:
            print(f"  ✓ {species}")

    if failed:
        print("\nFailed to download:")
        for species in failed:
            print(f"  ✗ {species}")

    # Print accession dictionary for reference
    if accession_dict:
        print("\n" + "=" * 70)
        print("ACCESSION NUMBERS (for your reference)")
        print("=" * 70)
        print("species_data = {")
        for species, accession in accession_dict.items():
            print(f"    '{species}': '{accession}',")
        print("}")

    print("\n" + "=" * 70)
    print(f"✓ Complete! Files saved in: {output_directory}/")
    print("=" * 70)


if __name__ == "__main__":
    main()