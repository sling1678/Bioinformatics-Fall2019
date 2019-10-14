import csv

def ambiguity_codes(filename='ambiguity_codes.csv'):
    """
    Read in ambiguity_codes.csv
    """
    codes = {}
    with open(filename) as csvfile:
        csv_reader_handle = csv.reader(csvfile, delimiter='\t')
        header_row = True
        for row in csv_reader_handle:
            if header_row: # skip the header row in the csv file
                header_row = False
                continue
            row[0] = row[0].rstrip()
            row[1] = row[1].rstrip().replace(' or ', ',').split(',')
            row[1] = [item for item in row[1] if len(item)!= 0]
            key = row[0]
            value = row[1]
            codes[key] = value
    return codes

if __name__ == "__main__":
    codes = ambiguity_codes()
    print(codes)


