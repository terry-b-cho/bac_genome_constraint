import collections

def load_full_ontology(term_list=None):

    # Create empty GO term list if one is not passed
    if term_list is None:
        term_list = []

    # Dict to hold data for all GO terms
    all_terms_data = {}

    # Open GO ontology file
    with open("/n/scratch/users/a/aip485/bac_genome_constraint/data/go/go-basic.obo", "r") as rf:

        # Iterate over lines
        for line in rf:

            # Unpack line data
            row = line.strip().split(": ", 1)

            # If line contains data, record it
            if len(row) == 2:
                key = row[0]
                value = row[1]

                # If the data is the term ID, create a new container for the data in the following rows
                if key == "id":
                    term_data = collections.defaultdict(list)

                    # If the term is in the term list or no term list was passed, keep the data for this row
                    if len(term_list == 0) or (value in term_list):
                        all_terms_data[value] = term_data
                
                # Record the data for this row
                term_data[key].append(value)

    return all_terms_data


def load_ids_and_names(term_list=None):
    # Create empty GO term list if one is not passed
    if term_list is None:
        term_list = []

    # Dict to hold data for all GO terms
    all_terms_data = {}

    # Open GO ontology file
    with open("/n/scratch/users/a/aip485/bac_genome_constraint/data/go/go-basic.obo", "r") as rf:

        # Iterate over lines
        for line in rf:

            # Unpack line data
            row = line.strip().split(": ", 1)

            # If line contains data, record it
            if len(row) == 2:
                key = row[0]
                value = row[1]

                if key == "id":
                    id = value

                elif key == "name":
                    name = value
                    all_terms_data[id] = name

    return all_terms_data