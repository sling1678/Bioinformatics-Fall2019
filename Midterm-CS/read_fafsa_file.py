#########################################################################
# Originally written by John Beck.
# FASTA files are text files in the following format
# 0. First line has header tag ">" as first character
# 1. A line starting with ">" has the name of the gene followed by "\n"
# 2. Next N lines of string+"\n" of one letter DNA, RNA, or polypeptide
# 3. Maybe one or more empty lines
# 4. Next DNA, RNA, or polypeptide sequence
# 5. Repeated till end of file
########################################################################

def process_header_line(line):
  """
  Parses header line for the gene and extracts name of the gene. The header line starts with a identifying character '>', followed by label and then other information, which we do not use but keep them around.
  Parameters:
    line : string - the header line identified by the first character
  Returns:
    name : string - the name of the gene
    header_line : string - the full content of the line after first character
  """
  label, header_line = "", ""
  if len(line) > 1:
    header_line = line[1:]
    label = header_line.split(' ')[0]
  return label, header_line

def read_file(filename):
  """
  Read in FASTA-format file containing one or more sequences

  Parameters:
    filename: path/to/file/containing/the/data
  Returns:
    results: a list of [label, full header_line, DNA or RNA or AA-seq string] list of every mRNA data
  """
  results = []
  with open(filename, "r") as infile:
    # find the first line that has label of the gene
    line = "*" # dummy char to get the while loop started
    while line[0] != ">": #skip till you find the line with label - should be first line
      line = infile.readline()
    line = line.rstrip()
    label, header_line = process_header_line(line)
    # Now we process lines that have the sequence
    # There may be blank lines till we reach the next gene
    # This means we need to check each line for ">" tag and continue if we find blank line
    sequence = "" # the sequence accumulator
    for line in infile:
      line = line.rstrip()
      if line == "":
        continue # go to next line of data
      if line[0] == ">":
        results.append([label, header_line, sequence])
        # start new gene
        label, header_line = process_header_line(line)
        sequence = ""
        continue # go to the next line which starts the sequence
      sequence += line.upper() # add to develong sequence
    # the last gene read must be recorded now
    results.append([label, header_line, sequence])
  return results

def contains_only_valid_chars(sequence, valid_chars):
  if not sequence:
    return False
  for letter in sequence:
    if letter not in valid_chars:
      return False
  return True

def get_valid_sequences(filename, valid_chars):
  """
  Reads fafsa file and tests the sequences for valid characters. 
  Parameters:
    filename: string - path-to-fafsa-file
    valid_chars: list of valid chars in the sequence
  Returns:
    valid_sequences: list of [name, header_line, sequence] of valid characters
    invalid_sequences: list of [name, header_line, sequence] of invalid characters
  """

  valid_sequences = []
  invalid_sequences = []

  results = read_file(filename)
  for item in results:
    if contains_only_valid_chars(item[2], valid_chars):
      valid_sequences.append(item)
    else:
      invalid_sequences.append(item)

  return valid_sequences, invalid_sequences
  

if __name__ == "__main__":
  # filename = "MidCS1.txt"
  filename = "test.txt"
  valid_chars = ['A', 'T', 'G', 'C']
  results = read_file(filename)
  print(results)
  # for items in results:
  #   valid = contains_only_valid_chars(items[2], valid_chars)
  #   print(valid)
  v, i = get_valid_sequences(filename, valid_chars)
  print(v)
  print(i)
