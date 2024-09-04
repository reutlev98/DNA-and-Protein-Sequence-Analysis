# DNA, Protein Sequence, and Player Data Analysis

## Introduction

This project includes Python scripts for analyzing DNA sequences, protein sequences, and player data. The main functionalities include finding Simple Sequence Repeats (SSRs) in DNA, identifying and visualizing PROSITE patterns in protein sequences, and performing various data operations on player datasets.

## Installation

### Prerequisites

- Python 3.x
- [Biopython](https://biopython.org/) (for handling biological sequences)
- [Matplotlib](https://matplotlib.org/) (for data visualization)
- [Pandas](https://pandas.pydata.org/) (for data manipulation and analysis)

### Installation Steps

1. Clone the repository or download the files.
2. Install the required Python packages:

   ```bash
   pip install biopython matplotlib pandas
   ```

## Usage

### `exercise1.py`

This script provides a function to find Simple Sequence Repeats (SSRs) in a given DNA sequence.

#### How to Use

1. Import the function `srr_find` from the script:

   ```python
   from exercise1 import srr_find
   ```

2. Pass a DNA sequence string to the `srr_find` function:

   ```python
   dna_seq = "AGATAGATAGATAGATAG"
   repeats = srr_find(dna_seq)
   print(repeats)
   ```

3. The function will return a list of SSRs and their counts.

### `exercise2.py`

This script defines a `Cell` class that models a biological cell with a genome. It includes methods for finding SSRs within the genome.

#### How to Use

1. Import the `Cell` class from the script:

   ```python
   from exercise2_207385741 import Cell
   ```

2. Create an instance of the `Cell` class with a name, genome, and reading frame:

   ```python
   my_cell = Cell("Cell1", "AGATAGATAGATAGATAG", 1)
   ```

3. Use the `srr_find` method to find SSRs within the genome:

   ```python
   repeats = my_cell.srr_find(0)
   print(repeats)
   ```

4. The method will return a list of SSRs found at the specified index.

### `exercise3.py`

This script defines a `myData` class that performs various data operations on player datasets, such as filtering by height, birth year, and specific columns.

#### How to Use

1. Import the `myData` class from the script:

   ```python
   from exercise3_207385741 import myData
   ```

2. Create an instance of the `myData` class by passing the path to a CSV file containing player data:

   ```python
   data = myData("path_to_your_data.csv")
   ```

### `exercise4.py`

This script is designed to analyze protein sequences and identify PROSITE patterns. It also generates visualizations of the occurrence of these patterns in protein sequences.

#### How to Use

1. Run the script to analyze protein sequences and generate a bar chart of pattern occurrences:

   ```python
   python exercise4_207385741.py
   ```

2. The script will generate a bar chart showing the number of occurrences of different PROSITE patterns in the analyzed protein sequences.


## Output

- **`exercise1.py`**: Returns a list of SSRs found in the DNA sequence.
- **`exercise2.py`**: Returns a list of SSRs found within the genome of the cell at the specified index.
- **`exercise3.py`**: Provides various data analysis results, including player counts, filtered DataFrames, sorted lists, and statistical summaries.
- **`exercise4.py`**: Generates a bar chart visualizing the occurrences of PROSITE patterns in protein sequences.

