# K-mer Parser and DNA Sequence Matcher

This program was written for DayZero Diagnostics Software Coding Assessment

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Overview](#overview)
3. [How to Run](#how-to-run)
4. [Testing and Efficiency](#testin-and-efficiency)
5. [Next Steps](#next-steps)

## System Requirements

 - Python3 

## Overview

This program does 2 things primarily:

1. Parse and process fastq data

    - This is done with a K-mer parser class which contains methods to clean an input fastq file, get all unique k-mers of a user specified length (the project description asked for all unique 21-mers in the sample 'SP1.fastq' file) along with their frequency count, and finally output said data to a database using sqlite3.

2. Provide a match function to to find highly matching k-mers in a DNA sequence

    - This is done using the k-mer parser class to get all unique k-mers that are the same length as the target sequence. Next, a simple matching cost algorithm is used to determine if they meet the proper matching criteria.

## How to Run

To use this program, download this repo and run kmers.py to generate the requested data.db and print a sample output to match

If one would like to use the functionality of the k-mer parser class, one can construct another parser object with different inputs, and or specify a different length for the get unique k-mers function. Additionally, one can output new data to another database using the toDatabase method.

Furthermore, the match function can be used with any DNA sequences (provided length of kseq <= seq) and can be modified to take in a fastq file for analysis (in constructor of the parser) which I did when testing it's functionality.

## Testing and Efficiency

To test the program, I manually validated that the fastqCleaner acquired all the DNA sequences and verified their lengths, the get unique kmers function outputted sequences of the correct length with no duplicates, and that the toDatabase method effectively created a database with the requested data. Furthermore, I tested the match function with the sample input provided and did some analysis at scale by comparing a target kseq with the SP1.fastq file. Overall, it worked as expected, however it could benefit from more through data validation and unit testing.

Efficiency wise, I believe the parser and match function could handle large inputs decently well as all of the algorithms used to parse and match run in polynomial time. 

## Next Steps

- Provide a UI

- Better database integration

- More testing

- Custom sequence matching threshold
