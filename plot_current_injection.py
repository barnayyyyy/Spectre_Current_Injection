# By: John Barney, Lucas Nichols

#-----------------------------------------------------------------------
# Introduction
#-----------------------------------------------------------------------
# This script is for plotting the results of run_current_injection.py
# The plot() function could be implemented into the run_current_injection.py if all modules were imported

import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------------------------------                            
# INPUTS                         
#-----------------------------------------------------------------------  
OUTPUT_NODE                 = 's0'
REFERENCE_NODE              = 'vgnd'
TIME_INDEX                  = 1

# EXAMPLE PATHS
original_csv_path = 'project\\John\\RESULTS_original.csv'
injected_csv_path = 'project\\John\\UPSETS\\[node_I16.int1][i2_0.000425][dt2_1.1e-09][df2_1.5e-09].csv'

#-----------------------------------------------------------------------                            
# INPUTS                         
#-----------------------------------------------------------------------  
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

def plot(times1, results1, times2, results2): 
    """
    Function to plot two data sets for comparison.

    Args:
        times1 (list): A list of time values for the first data set.
        results1 (list): A list of result values for the first data set.
        times2 (list): A list of time values for the second data set.
        results2 (list): A list of result values for the second data set.
    """
    # Create a new figure and axis using matplotlib's subplots function
    fig, ax = plt.subplots()

    # Plot the 'original' data (times1 vs results1) as a line plot
    ax.plot(times1, results1, '-',label='original')

    # Plot the 'injected' data (times2 vs results2) as a line plot with points
    ax.plot(times2, results2, '.-', label='injected')

    # Set the title and labels for the plot
    ax.set_title('Voltage vs Time')
    ax.set_xlabel('Time')
    ax.set_ylabel('Voltage')

    # Add a legend to the plot
    ax.legend()

    # Display the plot
    plt.show()

#-----------------------------------------------------------------------                            
# Run                         
#-----------------------------------------------------------------------  

# Extract the time and result columns from CSV files
times1 = extract_column_from_csv(original_csv_path, column_index=TIME_INDEX)
results1 = extract_column_from_csv(original_csv_path, column_name=OUTPUT_NODE)
times2 = extract_column_from_csv(injected_csv_path, column_index=TIME_INDEX)
results2 = extract_column_from_csv(injected_csv_path, column_name=OUTPUT_NODE)

# plot(x1,y1,x2,y2)
plot(times1,results1,times2,results2)