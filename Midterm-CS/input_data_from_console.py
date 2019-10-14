def input_data_from_console( max_num=5):
  """
  Provides user to provide upto max_num of promoters.
  Parameters:
    max_num : number of different promoters to be provided
  Returns:
    promoters : list of list - each item [name, consensus, score] for each promoter

  """
  num_data_fields=3
  separation = "/"
  stop_char = '*'

  message_to_display = "Enter promoter_name{}consensus_sequence{}score separated by {}:\n".format(separation, separation, separation)
  error_message = "ERROR IN YOUR ENTRY: " + message_to_display + "\nTo stop type {} and enter".format(stop_char)
  starting_message_format = message_to_display + "(upto {} promoters; Type * as first character on a line and enter to end)"  


  print(starting_message_format.format(max_num))

  promoters = []
  i = 0
  while i < max_num:
    mydata = input()
    if len(mydata)==0:
      print(error_message)
      continue
    if mydata[0] == stop_char:
      break
    else:
      mydata = mydata.split(separation)

      if len(mydata)<num_data_fields:
          print(error_message)
          i -= 1
      else: # all good to record
        try:
          name = mydata[0].upper()
          consensus = mydata[1].upper()
          score = float(mydata[2])
          promoters.append( [name, consensus, score] )
        except ValueError:
          print("Score must be a numerical value that can be converted to a float! Try again")
          i -= 1          

    i += 1
  return promoters


if __name__ == "__main__":
  promoters = input_data_from_console()
  print(promoters)