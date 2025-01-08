def format_taxonomy(line):
    """Transformation of string structure
    
    One line at the time, add prefix where needed and leave empty
    space when information is not available.
    """
    parts = line.strip().split(',')
    taxid = parts[0]
    # parts[1]: Scientific name (ex. Malpighiales)
    kingdom, phylum, class_name, order, family, genus, species = parts[2:9]
    taxonomic_levels = [
        ('k', kingdom),
        ('p', phylum),
        ('c', class_name),
        ('o', order),
        ('f', family),
        ('g', genus),
        ('s', species)
    ]
    tax_parts = []
    for prefix, value in taxonomic_levels:
        if value == '-':
            tax_parts.append(f"{prefix}__")
        else:
            tax_parts.append(f"{prefix}__{value}")
    taxonomy = ';'.join(tax_parts)
    return f"{taxid}\ttax={taxonomy}"


with open('taxa_info.txt', 'r') as infile, open('formatted_taxa.tsv', 'w') as outfile:
    for line in infile:
        formatted_line = format_taxonomy(line)
        outfile.write(formatted_line + '\n')
