Program description
-------
dNWA.py is a command-line tool written in Python that accepts a specific DomainMapper output file with sequence information as input. It transforms the input format and aligns all given sequences pairwise using a modified Needleman-Wunsch algorithm. The tool then generates an output file with the specified filename and a timestamp, containing the alignments.  

How to use
-------
To use the program you have to provide an input file in the same directory or with the path specified. Example of use:

```sh
python dNWA.py inputfile.txt
```
Input format
-------
The input file must conform to a specific structure (for the script to parse DomainMapper output accordingly see the GitHub repository at https://github.com/bsarah/phylocog). The file format must be plain text. An example of a file with one cluster and one sequence is provided below:
```sh
>ClusID 0; #Accessions 1; Structural_annotation L=L
accessionname (0,0,START) (0,5,F-ID1) (10,15,F-ID2) (15,15,END)
```
The file contains one line for each cluster and one line for each accession. The line beginning with ">" contains information about the cluster, the number of accessions, and the structural annotation. The annotation for protein domains is separated by a "=" sign. All lines not beginning with a ">" contain protein sequence data along with the accession name. Reading data from a file depends on the structure of the input file. If the structure is not met, the function of the code cannot be ensured.

Output format
-------

### intermediate file
This document is used as input for the modified Needleman-Wunsch algorithm. It contains the mapped protein domain sequences from the provided file.

### output file 
The output file contains all alignments and the information about the sequences in the following representation: 
```sh
> \>ClusID of seq1,accession of seq1,>ClusID of seq2,accession of seq2,alignmentscore,length of alignment,normalized alignmentscore
> Alignment
```

The normalized alignment score is calculated by: alignment score / (length of alignment - number of zero matches).

Options
-------

| Option | Usage |
| ------ | ------ |
| filename | you need to provide an input file |
| -h, -\-help | displays help |
| -v, -\-verbose | writes a second output file that is more human-readable |
| -a, -\-noselfalignment| deactivates the computation of alignments like (a,a) |
| -l, -\-log | sets the log level to [DEBUG, INFO, WARNING, ERROR, CRITICAL] |
| -t, -\-temp | keeps the intermediate file |
| -g, -\-gapextension | different penalties for gap openings and gap extensions |
| -s, -\-score | pass 9 score values to overwrite default values |

Set the logging level option like -\-log=INFO.

To set the score provide 9 weights in the following order: 
1) "weight for matching f-group"
2) "weight for matching t-group"
3) "weight for matching h-group" 
4) "weight for matching x-group"
5) "weight for matching no f-id-label"
6) "weight for a mismatch of f-ids"
7) "weight for mismatch no-fid with fid-label"
8) "gap penalty for opening"
9) "gap penalty for extension"

Examples
-------
Some tests, along with their respective results, can be found in the repository's "example_data" directory.
 
Troubleshooting instructions
-------
For purposes of troubleshooting, the script should be executed with the "-l=DEBUG" option, which will result in the generation of a full log file. The log contains metrics that can assist in the troubleshooting process.

Credits and License
-------

Further information on the original Needleman-Wunsch algorithm can be found in the following reference: https://doi.org/10.1016/0022-2836(70)90057-4. The code in this repository is released under the MIT license. See [LICENSE](https://github.com/blaueste/bachelor_thesis/blob/main/LICENSE).
