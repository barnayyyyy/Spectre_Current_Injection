# By: John Barney, Lucas Nichols

#----------------------------------------------------------------------------------------------------------------------------------------------
# README
#----------------------------------------------------------------------------------------------------------------------------------------------
"""
This program simulates circuit behavior using Spectre, injecting a current source into specified nodes,
and compares the results to a baseline simulation. It is used to identify vulnerabilities to current 
injections in the circuit design.

The program works by:
1. Reading a netlist file and a settings file.
2. Running a baseline simulation.
3. Running simulations with current injections into specified nodes.
4. Comparing the results of the injected simulations to the baseline simulation.
5. Identifying any "upsets" where the circuit behavior deviated significantly from the baseline.

Prerequisites/Dependencies:
- Spectre circuit simulator
- Python libraries: numpy, itertools, os

Input: 
- Netlist file (circuit description in Spectre format)
- Settings file (simulation settings in Spectre format)
- Parameters for current injection and tolerances for comparison

Output:
- CSV files of the baseline and injected simulation results
- CSV files of any upsets identified
- Plots comparing the baseline and injected results
"""

#-----------------------------------------------------------------------
# Introduction
#-----------------------------------------------------------------------
# This script is for injecting a current source into a Cadence Spectre netlist
# and using that current source to find sensitive nodes, as well as the critical charge amount
# the script is called "run_spectre.py" on the lillywhite computer

# tested on:    Library name: scc9gena
#               Cell name: scc9gena_parts_dfr
# upset nodes:  ['I3.net34','I3.net42','I16.int1','I16.int2','I17.int1','M1','s0']

#-----------------------------------------------------------------------
# Generating Netlist & Obtaining _graphical_stimuli
#-----------------------------------------------------------------------
# 1: go to circuit schematic
# 2: go into the ADEL window
# 3: Setup >> Stimuli : set up simulation parameters
# 4: Analysis >> Choose : choose analysis type (tran is most likely best)
# 5: Simulation >> Netlist >> Create : generate the netlist and save it
# 6: Get the _graphical_stimuli.scs file from the following directory:
#       /home/YOURNAME/simulation/YOURPROJECT/spectre/schematic/netlist

#-----------------------------------------------------------------------                            
# INPUTS                         
#-----------------------------------------------------------------------                  
# MUST RUN IN TERMINAL:
#           "source /home/YOURNAME/SW_Run/bashrc.IC"

NETLIST                     = '/home/scratch/John/scc9gena_parts_dfr.net'   # STRING : setting template file
OUTPUT_NODE                 = 's0'                                          # STRING : node that is set as output
REFERENCE_NODE              = 'vgnd'                                        # STRING : node to take reference from for voltage
TRAN_TIME                   = '4n'                                          # STRING IN SPECTRE NOTATION : how long transient test is run

CURRENT_INIT                = 0                                             # %I1%   : Current starting amperage 
CURRENT_STOP_START          = 0.000325                                                                                                              
CURRENT_STOP_STOP           = 0.000425                                      # %I2%   : Current ending amperage                  Note: %I2% - %I1% = current peak
CURRENT_STOP_INCR           = 0.000000                                      # increment going from I1 to I2 (AMP)

DELAY_INIT                  = 0                                             # %DT0%  : Delay Time for start of current pulse to exist in netlist (default = 0)
DELAY_TIME_START            = 0.000000001                                   # %DT1%  : Delay Time for current pulse to start
DELAY_TIME_STOP             = 0.0000000011                                  # %DT2%  : Delay Time for current pulse to stop     Note: %DT2% - %DT1% = pulse width
DELAY_TIME_INCR             = 0.0000000000                                                                                         # increment going from D1 to D2 (SEC)

DAMPING_FACTOR_RISE         = 0.000000000001                                # %DF1%  : Damping Factor on rising edge
DAMPING_FACTOR_FALL_START   = 0.000000001500                                # %DF2%  : Damping Factor on falling edge
DAMPING_FACTOR_FALL_STOP    = 0.000000001500
DAMPING_FACTOR_FALL_INCR    = 0.000000000000                                                                            

THRESHOLD_ABOVE_Y           = 1                                             # FLOAT : difference between OUTPUT_NODE - REFERENCE_NODE to save
TOLERANCE_Y                 = THRESHOLD_ABOVE_Y                             # FLOAT : variance between OUTPUT_NODE - REFERENCE_NODE to save file
THRESHOLD_AFTER_TIME        = DELAY_TIME_START                              # FLOAT : wait how long to start comparing injected vs original data
TOLERANCE_TIME              = 0.000000005                                   # FLOAT : variance between injected time vs original time per element

TIME_INDEX                  = 1                                             # INTEGER : position number of time column in output CSV file
ONLY_SAVE_LOWEST_UPSET      = True                                          # BOOL : variable to only save lowest upset parameters or continue

#----------------------------------------------------------------------------------------------------------------------------------------------                            
# PROGRAM                         
#----------------------------------------------------------------------------------------------------------------------------------------------   
import os                                                                                                               # for operating system
import itertools
import numpy as np

#-----------------------------------------------------------------------
# Shorthands 
#-----------------------------------------------------------------------
# %NODE%            : variable for node current source injects to
# %REFERENCE_NODE%  : variable for node current source uses as reference

HOMEPATH = os.path.dirname(NETLIST)                                                                                     # finding working directory
TEMPLATE = HOMEPATH+'/TEMPLATE.tmpl' 

# spectre netlist settings
SETTINGS = 'Settings options rawfmt=nutascii\n\
include "'+HOMEPATH+'/_graphical_stimuli.scs"\n\
tran tran stop='+TRAN_TIME+' writefinal="spectre.fc" annotate=status maxiters=5\n'

# current source to be injected
ISOURCE = '// By: John Barney\n\
// Current Source templace to inject into Spectre Netlists\n\
// Library name: MAINLIB_TESTING\n\
// Cell name: CurrentSource_John\n\
// View name: schematic\n\
parameters DT0=%DT0% DT1=%DT1% DT2=%DT2% DF1=%DF1% DF2=%DF2% I1=%I1% I2=%I2%\n\
INJECTED_CURRENT (%REFERENCE_NODE% %NODE%) isource type=exp delay=DT0 val0=I1 val1=I2 td1=DT1 tau1=DF1 td2=DT2 tau2=DF2\n'

class OnlySaveLowestUpset(Exception):
    pass

#-----------------------------------------------------------------------
# Functions 
#-----------------------------------------------------------------------
# building block functions
def get_nodes_from_spectrefc(spectrefc_file=''):
    """
    This function takes a spectre.fc file as input and retrieves the nodes from it. If no file is provided, 
    the function prompts the user to input a file. The function parses the input file and collects node 
    names and final values. It then prepares a filtered list of node names for output.
    
    Args:
        spectrefc_file (str): String containing the path to the spectre.fc file to be processed. 
                              Defaults to an empty string.
    
    Returns:
        list: List of node names extracted from the spectre.fc file.
    """
    
    if spectrefc_file == '':
        print('Input a "spectre.fc" file')                                                                              # If no file is input, request for one.
    else:                                                               
        with open(spectrefc_file, 'r') as file:                                                                         # Open the specified file in read mode.
            lines = file.readlines()                                                                                    # Read all lines from the file into a list.

            data = {'Node Name': [], 'Final Value': []}                                                                 # Initialize a dictionary to store the node names and their final values.

            for line in lines:                                                              
                line = line.split('#')[0].strip()                                                                       # Remove comments from each line by splitting on '#' and taking the first part.

                if line:                                                                                                # If the line is not empty after removing comments.
                    parts = line.split('\t')                                                                            # Split the line into parts using tab as a delimiter.

                    if len(parts) == 2:                                                                                 # If there are exactly 2 parts (node name and final value).
                        data['Node Name'].append(parts[0].strip())                                                      # Strip whitespace and append node name to the data dictionary.
                        data['Final Value'].append(parts[1].strip())                                                    # Strip whitespace and append final value to the data dictionary.

        filtered_data = {'Node Name': [], 'Final Value': []}                                                            # Initialize another dictionary to store filtered node names and final values.

        for i in range(len(data['Node Name'])):                                                                         # For each node name in the data dictionary.
            node_name = data['Node Name'][i]                                                                            # Get the node name.
            final_value = data['Final Value'][i]                                                                        # Get the final value.
            
            filtered_data['Node Name'].append(node_name)                                                                # Append the node name to the filtered_data dictionary. Currently, no filtering is applied.
            filtered_data['Final Value'].append(final_value)                                                            # Append the final value to the filtered_data dictionary. Currently, no filtering is applied.

        nodes = filtered_data['Node Name']                                                                              # Get the list of node names from the filtered_data dictionary.
    
        return nodes                                                                                                    # Return the list of filtered node names.

def fixraw_trans(datfile):
    """
    This function converts the 'rawfmt=nutascii' formatted output file (from a SPICE-like simulator) to a more human-readable,
    space-delineated format. It reads the original file, processes each line, and writes the processed data to a new file.

    Args:
        datfile (str): The path of the input data file.

    Returns:
        str: The path of the output file containing the processed data.
    """

    num_var = 7                                                                                                         # Initialize num_var to 7 which represents the total number of variables plus 1.

    tmpstr = ''                                                                                                         # Initialize a temporary string to concatenate and store variable values.

    k = 0                                                                                                               # Initialize a counter variable.

    val = 0                                                                                                             # Initialize a flag to track the start of value lines.
    tran = 0                                                                                                            # Initialize a flag to track the start of Transient Analysis section.

    result = []                                                                                                         # Create an empty list to store the processed data.

    with open(datfile, 'r') as infile:                                                                                  # Open the input file in read mode.
        for line in infile:                                                                                             # Iterate over each line in the input file.
            if val and tran:                                                                                            # If we are in the section containing the transient analysis data.
                linlist = line.split()                                                                                  # Split the line into a list of words.
                for i in linlist:                                                                                       # For each word in the list.
                    tmpstr = tmpstr + i + ' '                                                                           # Append the word followed by a space to the temporary string.
                    if ((k + 1) % num_var) == 0 and not k == 0:                                                         # If the current word is the last variable of a group.
                        result.append(tmpstr)                                                                           # Append the temporary string to the result list.
                        tmpstr = ''                                                                                     # Reset the temporary string.
                    k += 1                                                                                              # Increment the counter.
            if line.find("Transient") >= 0:                                                                             # If the line marks the start of the Transient Analysis section.
                tran = 1                                                                                                # Set the tran flag to 1.
            if line.find("Value") >= 0 and tran == 1:                                                                   # If the line marks the start of the value lines.
                val = 1                                                                                                 # Set the val flag to 1.
            if line.find("No. Variables") >= 0 and tran == 1:                                                           # If the line specifies the number of variables.
                y = line.split()                                                                                        # Split the line into a list of words.
                num_var = int(y[2]) + 1                                                                                 # Set num_var to the number of variables plus 1.

    # Create a filename for the output file by appending '_processed.raw' to the name of the input file.
    output_file = os.path.splitext(datfile)[0] + '_processed.raw'

    with open(output_file, 'w') as outfile:                                                                             # Open the output file in write mode.
        outfile.write('\n'.join(result))                                                                                # Write the processed data to the output file.

    return output_file                                                                                                  # Return the path of the output file.

def ssv_to_csv(input_file, output_file, header_list=[]):
    """
    This function converts a space-separated value (SSV) file to a comma-separated value (CSV) file.
    If a header list is provided, it adds the header to the CSV file using the add_header_to_csv function.
    
    Args:
        input_file (str): Path to the input SSV file.
        output_file (str): Path to the output CSV file.
        header_list (list, optional): List of strings that forms the header of the CSV file. Defaults to an empty list.
        
    Returns:
        str: Path to the output CSV file.
    """
    
    with open(input_file, 'r') as file_in, open(output_file, 'w') as file_out:                                          # Open the input file in read mode and the output file in write mode.
        for line in file_in:                                                                                            # Iterate over each line in the input file.
            row = line.strip().split(' ')                                                                               # Remove leading/trailing whitespace from the line and split it into a list of words.
            csv_row = ','.join(row)                                                                                     # Combine the words into a single string with commas in between.
            file_out.write(csv_row + '\n')                                                                              # Write the comma-separated string to the output file followed by a newline character.

    if header_list != []:                                                                                               # If a header list is provided.
        output_file = add_header_to_csv(output_file, header_list)                                                       # Add the header to the CSV file using the add_header_to_csv function.
                                                        
    return output_file                                                                                                  # Return the path to the output CSV file.

def add_header_to_csv(input_file, header_list):
    """
    This function adds a header to an existing CSV file. It ensures that the header is properly aligned with the data rows
    by adding a necessary amount of leading commas to the header.
    
    Args:
        input_file (str): Path to the CSV file to be updated.
        header_list (list): List of strings that forms the header of the CSV file.
        
    Returns:
        str: Path to the updated CSV file.
    """
    
    with open(input_file, 'r') as file:                                                                                 # Open the input file in read mode.
        content = file.readlines()                                                                                      # Read all lines from the file into a list.

    # If the content list is empty, raise an exception because we can't add a header to an empty file.
    if not content:
        raise ValueError('Input file is empty')

    # Find the length of the longest row in the file.
    max_row_len = max(len(row.split(',')) for row in content)

    # Calculate the number of commas to add before the header list.
    num_commas = max_row_len - len(header_list)

    # Create a string of commas to add before the header list.
    comma_str = ',' * num_commas

    # Add the padded header row at the beginning of the file content.
    header_str = comma_str + ','.join(header_list) + '\n'
    new_content = header_str + ''.join(content)

    # Write the updated content back to the same file.
    with open(input_file, 'w') as file:
        file.write(new_content)
        
    return input_file                                                                                                   # Return the path to the updated CSV file.

def extract_column_from_csv(file_path, column_name=None, column_index=None):
    """
    This function extracts the elements of a specified column from a CSV file and returns them as a list. The column can 
    be specified either by its name or by its index.

    Args:
        file_path (str): The path to the CSV file.
        column_name (str, optional): The name of the column to be extracted. Defaults to None.
        column_index (int, optional): The index of the column to be extracted. Defaults to None.

    Raises:
        ValueError: If the column specified by column_name or column_index is not found in the CSV file.

    Returns:
        list: A list containing all the elements of the specified column.
    """

    with open(file_path, 'r') as file:                                                                                  # Open the CSV file in read mode.
        lines = file.readlines()                                                                                        # Read all lines from the file into a list.

        header = lines[0].strip().split(',')                                                                            # Get the header row from the file and split it into a list of column names.

        try:                                                                                                            # Try to find the specified column in the header row.
            if column_index is not None:                                                                                # If a column index is specified.
                column_name = header[column_index]                                                                      # Get the name of the column at the specified index.
            else:                                                                                                       # If a column name is specified.
                column_index = header.index(column_name)                                                                # Get the index of the specified column.
        except ValueError:                                                                                              # If the specified column is not found in the header row.
            raise ValueError("Column '{}' not found in the CSV file.".format(column_name))

        column_list = []                                                                                                # Initialize an empty list to store the elements of the specified column.

        for line in lines[1:]:                                                                                          # Iterate over each line in the file, skipping the header row.
            values = line.strip().split(',')                                                                            # Split the line into a list of values.
            column_element = values[column_index]                                                                       # Get the element at the specified column index.
            column_list.append(column_element)                                                                          # Add the element to the column list.

    column_list = [float(x) for x in column_list]                                                                       # Convert the elements of the column list to floats.

    return column_list                                                                                                  # Return the column list.

def find_closest_time(t, times_array):
    """
    This function finds the index of the time value in a given array that is closest to a specified time.
    
    Args:
        t (float): The specified time.
        times_array (numpy.ndarray): A numpy array of time values.
        
    Returns:
        int: The index of the time value in times_array that is closest to t.
    """
    
    # np.abs(times_array - t) creates a new array where each element is the absolute difference between t and the corresponding element in times_array.
    # np.argmin() then finds the index of the smallest value in this new array, which is the index of the time value in times_array that is closest to t.
    return np.argmin(np.abs(times_array - t))

# calls other functions
def raw_to_csv(datfile, csv_file, header_list=[]):
    """
    This function transforms a raw data file into a CSV file. It first converts the raw file into a more readable space-separated
    values (SSV) format, then converts the SSV file into a CSV file, and finally adds a header to the CSV file.

    Args:
        datfile (str): The path to the raw data file.
        csv_file (str): The path to the output CSV file.
        header_list (list, optional): A list of strings that forms the header of the CSV file. Defaults to an empty list.

    Returns:
        str: The path to the output CSV file.
    """
    
    # Step 1: Convert the raw file to a more readable format using the fixraw_trans function.
    # The fixraw_trans function returns the path to the fixed file, which is stored in the fixed_file variable.
    fixed_file = fixraw_trans(datfile)

    # Step 2: Convert the fixed space-separated values (SSV) file into a CSV file using the ssv_to_csv function.
    # The ssv_to_csv function takes the path to the SSV file and the path to the output CSV file as inputs, and also 
    # optionally a header list. The output CSV file will be written to the path specified by the csv_file variable.
    ssv_to_csv(fixed_file, csv_file, header_list)

    # The path to the output CSV file is returned by the function.
    return csv_file

def analyze_and_compare_data(original_file_path, injected_file_path, time_column, result_column, time_tolerance, result_tolerance, time_threshold=None, result_threshold=None):
    """
    This function compares data between two CSV files given a tolerance for time and result columns. It can filter data 
    based on time and result thresholds.

    Args:
        original_file_path (str): The path to the original CSV file.
        injected_file_path (str): The path to the CSV file with injected data.
        time_column (int): The index of the time column in the CSV files.
        result_column (str): The name of the result column in the CSV files.
        time_tolerance (float): The tolerance for comparing time values.
        result_tolerance (float): The tolerance for comparing result values.
        time_threshold (float, optional): The minimum time value to consider. Defaults to None.
        result_threshold (float, optional): The minimum result value to consider. Defaults to None.

    Returns:
        bool: True if all compared values are within the specified tolerances, False otherwise.
    """
    global times1, results1, times2, results2

    # Extract the time and result columns from the first file
    times1 = extract_column_from_csv(original_file_path, column_index=time_column)
    results1 = extract_column_from_csv(original_file_path, column_name=result_column)

    # Extract the time and result columns from the second file
    times2 = extract_column_from_csv(injected_file_path, column_index=time_column)
    results2 = extract_column_from_csv(injected_file_path, column_name=result_column)

    # Convert times1 to a numpy array for efficient computations
    times1_array = np.array(times1)

    # Filter the data from the second file based on the thresholds
    if time_threshold is not None and result_threshold is not None:
        times2, results2 = zip(*[(t, r) for t, r in zip(times2, results2) if t >= time_threshold and r >= result_threshold])

    # Compare the data
    for t2, r2 in zip(times2, results2):
        # Find the index of the closest time in times1 to t2
        idx = find_closest_time(t2, times1_array)

        # Get the corresponding time and result from the first file
        t1 = times1[idx]
        r1 = results1[idx]

        # Compare the values
        if abs(t1 - t2) > time_tolerance or abs(r1 - r2) > result_tolerance:
            return False                                                                                                # If either difference is greater than its respective tolerance, return False

    return True                                                                                                         # If all differences are within their respective tolerances, return True

# Main functions
def run_netlist_setup(netlist, settings, isource):
    """
    This function prepares and runs the Spectre simulation for a given netlist, settings and a current source to be injected. 

    Args:
        netlist (str): The path to the netlist file.
        settings (str): The Spectre simulation settings.
        isource (str): The injected current source for fault injection.

    Global Variables:
        NODES (list): The list of nodes obtained from the spectre.fc file.
        original_csv_path (str): The path to the original output file in CSV format.
        injected_template_data (str): Contents of netlist after settings and isource have been injected.
    """
    global NODES, original_csv_path, injected_template_data

    # ORIGINAL
    with open(netlist,'r') as f:                                                                                        # Open the netlist file for reading
        NETLIST_DATA = f.read()\
        .replace('"design_wrapper.lib"','"/home/PDKs/Skywater/c9pdk/MODELS/SPECTRE/c9fh-3r/Models/design_wrapper.lib"') # Replace the path to the library to make it runnable
    with open(TEMPLATE,'w') as f:                                                                                       # Open the template file for writing
        f.write(NETLIST_DATA+settings)                                                                                  # Write the combined netlist data and settings to the template
    os.system('spectre {} > {}/RESULTS.raw && mv TEMPLATE.raw {}/RESULTS.raw'.format(TEMPLATE, HOMEPATH, HOMEPATH))     # Run the spectre simulation and move the output to the results folder
    NODES = get_nodes_from_spectrefc(HOMEPATH+'/spectre.fc')                                                            # Extract the nodes from the spectre.fc file
    original_csv_path = raw_to_csv('RESULTS.raw',HOMEPATH+'/RESULTS_original.csv',NODES)                                # Convert the raw output to CSV format and save it

    # INJECTED CURRENT SOURCE
    command = 'if [ -d UPSETS ]; then rm -r UPSETS; fi && mkdir UPSETS'                                                 # Check if the UPSETS directory exists, if it does, remove it, then create it
    os.system(command)
    
    with open(TEMPLATE,'w') as f:
        f.write(NETLIST_DATA+settings+isource)                                                                          # Write the combined netlist data, settings and injected current source to the template
    with open(TEMPLATE,'r') as f:
        injected_template_data = f.read()                                                                               # Read the injected template

def run_injected_isource(tested_nodes):
    """
    This function simulates and analyzes the effect of injecting a current source at different nodes for different current amplitudes, delay times and damping factors.

    Args:
        tested_nodes (list): The list of nodes to test.

    Global Variables:
        Various parameters such as CURRENT_START, CURRENT_STOP, CURRENT_INCR, DELAY_TIME_START, DELAY_TIME_STOP, DELAY_TIME_INCR, DAMPING_FACTOR_FALL_START, DAMPING_FACTOR_FALL_STOP, DAMPING_FACTOR_FALL_INCR, REFERENCE_NODE, DELAY_INIT, DAMPING_FACTOR_RISE, TIME_INDEX, OUTPUT_NODE, TOLERANCE_TIME, TOLERANCE_Y, THRESHOLD_AFTER_TIME, THRESHOLD_ABOVE_Y are global constants defined elsewhere in the code.
    """
    
    # Create lists of current amplitudes, delay times, and damping factors to test
    if CURRENT_STOP_INCR == 0:
        currents = [CURRENT_STOP_STOP]
    else:
        currents = np.arange(CURRENT_STOP_START, CURRENT_STOP_STOP + CURRENT_STOP_INCR, CURRENT_STOP_INCR).tolist()

    if DELAY_TIME_INCR == 0:
        delay_times = [DELAY_TIME_STOP]
    else:
        delay_times = np.arange(DELAY_TIME_START, DELAY_TIME_STOP + DELAY_TIME_INCR, DELAY_TIME_INCR).tolist()

    if DAMPING_FACTOR_FALL_INCR == 0:
        damping_factors = [DAMPING_FACTOR_FALL_STOP]
    else:
        damping_factors = np.arange(DAMPING_FACTOR_FALL_START, DAMPING_FACTOR_FALL_STOP + DAMPING_FACTOR_FALL_INCR, DAMPING_FACTOR_FALL_INCR).tolist()

    # Converting elements in the lists to strings
    currents = [str(i) for i in currents]
    delay_times = [str(i) for i in delay_times]
    damping_factors = [str(i) for i in damping_factors]

    # Generate all combinations of nodes, current amplitudes, delay times and damping factors to test
    combinations = list(itertools.product(tested_nodes, currents, delay_times, damping_factors))
    total = len(combinations)

    print('\n# of total : node, current, delay_time, damping_factor')
    for i, (node, current, delay_time, damping_factor) in enumerate(combinations, start=1):
        print(i, "of", total, ":", node, current, delay_time, damping_factor)

        # Insert the current combination into the template
        injected_netlist_edit = injected_template_data\
                        .replace('%REFERENCE_NODE%', str(REFERENCE_NODE)).replace('%NODE%', str(node))\
                        .replace('%I1%', str(CURRENT_INIT)).replace('%I2%', str(current))\
                        .replace('%DT0%', str(DELAY_INIT        )).replace('%DT1%', str(DELAY_TIME_START)).replace('%DT2%', str(delay_time))\
                        .replace('%DF1%', str(DAMPING_FACTOR_RISE)).replace('%DF2%', str(damping_factor))

        # Write the edited template to a file
        with open(TEMPLATE, 'w') as f:
            f.write(injected_netlist_edit)

        try:
            # Run the spectre simulation with the edited template
            os.system('spectre {} > {}/RESULTS.raw && mv TEMPLATE.raw {}/RESULTS.raw'.format(TEMPLATE, HOMEPATH, HOMEPATH))
        except:
            continue

        # Convert the raw output to CSV format and save it
        injected_csv_path = raw_to_csv('RESULTS.raw',HOMEPATH+'/RESULTS_injected.csv',NODES)

        # Analyze and compare the injected data to the original data
        overlap = analyze_and_compare_data(original_csv_path, injected_csv_path, TIME_INDEX, OUTPUT_NODE, TOLERANCE_TIME, TOLERANCE_Y, THRESHOLD_AFTER_TIME, THRESHOLD_ABOVE_Y)
        
        # If an upset is detected, save the CSV file in the UPSETS folder
        if not overlap:
            # If the analysis of the injected data compared to the original data doesn't overlap within the specified tolerance, an upset is identified.
            print('UPSET')
            # Format the save name to include the node, current, delay time and damping factor values.
            savename = 'UPSETS/[node_{}][i2_{}][dt2_{}][df2_{}]'.format(node,current,delay_time,damping_factor)
            # Use a system call to copy the CSV file of the injected data to the UPSETS folder with the new save name.
            os.system('cp {} {}.csv'.format(injected_csv_path,savename))

            # plot()

            if ONLY_SAVE_LOWEST_UPSET == True:
                raise OnlySaveLowestUpset('ONLY_SAVE_LOWEST_UPSET == True')

#----------------------------------------------------------------------------------------------------------------------------------------------
# Run Procedure
#----------------------------------------------------------------------------------------------------------------------------------------------
# provide instructions for operation
print('CTRL + C : Stops current node iteration\n\
CTRL + Z : Stops entire program')

# 1: change to working directory
os.chdir(HOMEPATH)                          

# 2: run setup procedure
run_netlist_setup(NETLIST,SETTINGS,ISOURCE)

# OPTIONAL EDIT BLOCK
# ----------------------------------------------------------------------------------------------------------------------------------------------
# EXAMPLE OF NODE FILTERING
# NODES_filtered = [node for node in NODES if "_" not in node and ":" not in node]                                                # remove nodes that do not work with spectre (known errors)
# NODES_filtered.remove('D')
# NODES_filtered.remove('I16.int2')
NODES_filtered = NODES
#----------------------------------------------------------------------------------------------------------------------------------------------

# 3: run current injection test
for test_node in NODES_filtered:                                                                                                # suggested way to input nodes list
    try:
        run_injected_isource([test_node])
    except Exception as e:
        print(e)
        continue

# 4: cleanup procedure
commands = ['rm -r TEMPLATE.ahdlSimDB','rm -r TEMPLATE.log','rm -r RESULTS_processed.raw','rm -r TEMPLATE.raw.psf','rm -r RESULTS.raw','rm -r spectre.fc','rm -r RESULTS_injected.csv','rm -r TEMPLATE.tmpl', 'rm -r TEMPLATE.tmpl.tran.srf']
for command in commands:
    os.system(command)

print('End of program')