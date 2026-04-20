import os
import pandas as pd
import time
from Bio import Entrez
import xml.etree.ElementTree as ET

from low_abundance_phyla import compute_low_abundance_phyla, normalize_taxon_name

# --- CONFIGURATION ---
Entrez.email = os.getenv("NCBI_EMAIL", "your.email@example.com")  # NCBI requires your email
API_KEY = os.getenv("NCBI_API_KEY")  # Optional: improves rate limits
if API_KEY:
    Entrez.api_key = API_KEY
if Entrez.email == "your.email@example.com":
    print("Warning: set NCBI_EMAIL for compliant NCBI usage.")

# Load your data
articles_df = pd.read_csv('unique_articles_ranked.csv')
taxa_df = pd.read_csv('taxa_ids_filtered.csv')

# Skip ranks that should not enter full-text scan / literature summarization.
EXCLUDED_RANKS = {"kingdom", "superkingdom", "domain"}
LOW_ABUNDANCE_PHYLA = {
    normalize_taxon_name(name).lower()
    for name in compute_low_abundance_phyla()
}
if "taxon_rank" in taxa_df.columns:
    before_n = len(taxa_df)
    ranks = taxa_df["taxon_rank"].astype(str).str.strip().str.lower()
    current_names = (
        taxa_df["current_scientific_name"]
        if "current_scientific_name" in taxa_df.columns
        else pd.Series([""] * len(taxa_df), index=taxa_df.index)
    )
    non_phylum_keep = ~ranks.isin(EXCLUDED_RANKS | {"phylum"})
    phylum_keep = ranks.eq("phylum") & (
        taxa_df["taxon"].astype(str).map(lambda x: normalize_taxon_name(x).lower()).isin(LOW_ABUNDANCE_PHYLA)
        | current_names.astype(str).map(lambda x: normalize_taxon_name(x).lower()).isin(LOW_ABUNDANCE_PHYLA)
    )
    taxa_df = taxa_df.loc[non_phylum_keep | phylum_keep].copy()
    print(
        f"Excluded {before_n - len(taxa_df)} taxa by rank filter: "
        f"{sorted(EXCLUDED_RANKS)} plus phyla not in low-abundance set"
    )

# --- HELPER FUNCTIONS ---

def get_article_text(article_id):
    """
    Fetches text from NCBI. 
    - For PMC IDs: Fetches full text XML and extracts all paragraphs.
    - For PMIDs: Fetches abstract.
    """
    try:
        if article_id.startswith('PMC'):
            # Fetch Full Text for PMC
            # Note: We strip 'PMC' from the ID for the request if needed, 
            # but Entrez often handles "PMC" prefix in db='pmc'.
            # It's safer to pass the number only for id, but let's try direct.
            clean_id = article_id.replace('PMC', '')
            handle = Entrez.efetch(db="pmc", id=clean_id, rettype="full", retmode="xml")
            record = handle.read()
            handle.close()
            
            # Parse XML to get all text content
            root = ET.fromstring(record)
            # Recursively find all text in the XML tree
            all_text = "".join(root.itertext())
            return all_text.lower()
            
        elif article_id.startswith('PMID') or article_id.isdigit():
            # Fetch Abstract for PubMed
            clean_id = article_id.replace('PMID', '')
            handle = Entrez.efetch(db="pubmed", id=clean_id, rettype="abstract", retmode="text")
            abstract = handle.read()
            handle.close()
            return abstract.lower()
            
    except Exception as e:
        print(f"Error fetching {article_id}: {e}")
        return "" # Return empty string on failure

def check_taxon_presence(text, taxon_name, scientific_name):
    """
    Checks if taxon name or scientific name appears in the text.
    """
    if not text:
        return 0
    
    # Simple substring search (case-insensitive done by lower() on text)
    # Add spaces around terms to avoid partial word matches (e.g. "actino" inside "actinobacteria")
    # This is a basic check; regex could be more robust for punctuation.
    
    t_name = taxon_name.lower()
    s_name = scientific_name.lower() if isinstance(scientific_name, str) else ""
    
    if t_name in text:
        return 1
    if s_name and s_name in text:
        return 1
        
    return 0

# --- MAIN EXECUTION ---

# Initialize the matrix
# Rows: Taxa, Columns: Articles
matrix = pd.DataFrame(0, index=taxa_df['taxon'], columns=articles_df['id'])

# Cache downloaded texts to avoid re-fetching if code is re-run or modified
article_texts = {}

print(f"Processing {len(articles_df)} articles...")

for art_id in articles_df['id']:
    print(f"Fetching {art_id}...")
    text = get_article_text(art_id)
    article_texts[art_id] = text
    
    # Be polite to the server (3 requests per second limit without API key)
    if not API_KEY:
        time.sleep(0.34) 

print("Scanning texts for taxa...")

# Iterate through every taxon and every article
for taxon_idx, taxon_row in taxa_df.iterrows():
    taxon = taxon_row['taxon']
    sci_name = taxon_row['current_scientific_name']
    
    for art_id in articles_df['id']:
        text = article_texts.get(art_id, "")
        
        # Check if present
        is_present = check_taxon_presence(text, taxon, sci_name)
        
        if is_present:
            matrix.loc[taxon, art_id] = 1

# Save the result
matrix.to_csv('taxon_article_fulltext_matrix.csv')
print("Done! Saved to taxon_article_fulltext_matrix.csv")
