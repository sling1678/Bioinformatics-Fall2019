"""
  FASTA files are text files in the following format
  0. First line has header tag ">" as first character
  1. A line starting with ">" has the name of the gene followed by "\n"
  2. Next N lines of string+"\n" of one letter ribonucleotide (AUGC)
  3. Maybe one or more empty line
  4. Next gene
  5. Repeated till end of file

"""

def process_header_line(line):
  label, header_line = "", ""
  if len(line) > 1:
    header_line = line[1:]
    label = header_line.split(' ')[0]
  return label, header_line

def read_file(filename):
  """
  Read in FASTA-format file containing one or more mRNA.

  Parameters:
    filename: path/to/file/containing/the/data
  Returns:
    results: a list of [name, header_line, mRNA] list of every mRNA data
  """
  results = []
  with open(filename, "r") as infile:
    # First line is label
    line = infile.readline().rstrip()
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
      sequence += line # add to develong sequence
    # the last gene read must be recorded now
    results.append([label, header_line, sequence])
  return results

if __name__ == "__main__":
  filename = "data/Assignment1Sequences.txt"
  results = read_file(filename)
  print(results)
