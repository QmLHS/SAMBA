def parse_taxonomy_files(nodes_file, names_file):
    """Parse files
    From file: nodes.dmp; names.dmp 
    To Dict: nodes_info; scientific_names
    """
    nodes_info = {}
    with open(nodes_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t|\t')
            tax_id = fields[0]
            parent_tax_id = fields[1]
            rank = fields[2]
            nodes_info[tax_id] = {
                'parent_tax_id': parent_tax_id,
                'rank': rank
            }
    scientific_names = {}
    with open(names_file, 'r') as f:
        for line in f:
            fields = [field.strip().strip('|').strip() for field in line.strip().split('\t|\t')]
            tax_id = fields[0]
            name = fields[1]
            name_class = fields[3]
            if name_class == "scientific name":
                scientific_names[tax_id] = name
    
    return nodes_info, scientific_names


def get_lineage_info(tax_id, nodes_info, scientific_names):
    """Get all lineage information 
    """
    lineage = []
    current_tax_id = tax_id
    visited = set()
    while current_tax_id in nodes_info and current_tax_id not in visited:
        visited.add(current_tax_id)
        node_info = nodes_info[current_tax_id]
        if current_tax_id in scientific_names:
            rank = node_info['rank'].lower()
            name = scientific_names[current_tax_id]
            lineage.append({
                'rank': rank,
                'name': name
            })
        # go to parent node until no more parent available or root reached
        current_tax_id = node_info['parent_tax_id']
        if current_tax_id == "1":
            break
    # Return in root -> leaf order       
    return lineage[::-1] 


def format_lineage(lineage_info):
    """Format lineage
    - only the 7 main ranks
    - including empty ones
    - starting with 'tax='
    """
    main_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    rank_prefixes = {
        'kingdom': 'k',
        'phylum': 'p',
        'class': 'c',
        'order': 'o',
        'family': 'f',
        'genus': 'g',
        'species': 's'
    }
    rank_to_name = {entry['rank']: entry['name'] for entry in lineage_info}
    formatted_parts = []
    for rank in main_ranks:
        prefix = rank_prefixes[rank]
        name = rank_to_name.get(rank, '')  # Get '' if rank not in dict.keys()
        formatted_parts.append(f"{prefix}__{name}") 
    return 'tax=' + ';'.join(formatted_parts)


def process_input_file(input_file, nodes_info, scientific_names):
    """Process the input file
    Add formatted taxonomy information
    """
    results = []
    line_number = 0
    with open(input_file, 'r') as f:
        for line in f:
            line_number += 1
            line = line.strip()
            # Skip empty lines
            if not line:
                print(f"Warning: Empty line at line {line_number}.")
                continue

            parts = line.split('\t')
            if len(parts) != 2:
                print(f"Warning: Line {line_number} does not have exactly two columns: '{line}'")
                continue
            accession, taxon = parts
            # Check if taxon starts with "taxon:"
            if not taxon.startswith('taxon:'):
                print(f"Warning: Line {line_number} has invalid taxon format: '{taxon}'")
                continue
            tax_id = taxon.split(':')[1]
            # Verify tax_id exists in nodes_info
            if tax_id not in nodes_info:
                print(f"Warning: Tax ID {tax_id} not found in taxonomy database at line {line_number}")
                continue
            lineage_info = get_lineage_info(tax_id, nodes_info, scientific_names)
            formatted_lineage = format_lineage(lineage_info)
            results.append({
                'accession': accession,
                'lineage': formatted_lineage
            })
    return results


nodes_file = "nodes.dmp"
names_file = "names.dmp"
input_file = "ids_and_taxon.txt"

print("Parsing taxonomy files")
nodes_info, scientific_names = parse_taxonomy_files(nodes_file, names_file)

print("Processing input file")
results = process_input_file(input_file, nodes_info, scientific_names)

print("Writing file ncbi_taxa.tsv")
with open("ncbi_taxa.tsv", "w") as f:
    for entry in results:
        f.write(f"{entry['accession']}\t{entry['lineage']}\n")

print("DONE")
