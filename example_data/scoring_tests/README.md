# Example Data
## Python Implementation Tests
To test correctness, robustness, and handling of incorrect input, there are 21 test files in the python_tests folder.
* Test 01: Empty File - This test consists of a 0KB txt file used as input. The input file is not shown here for obvious reasons. The expected output is an error message.
* Test 02: No Information - This test contains a text file with no ClusID, no accession number, no structural annotation and no sequence ID.
* Test 03: Only 1 Sequence - This tests the behavior of the program when there is only one sequence, so that an alignment that is not a self-alignment cannot be calculated.
* Test 04: Not all information - This test is similar to test 2, but includes structural annotation.
* Test 05: Long sequences - This test contains 2 clusters. The first cluster has 2 accessions, the second cluster has 4 accessions. Unlike the other tests, this test contains real sequences that are particularly long compared to the other sequences provided.
* Test 06: Only 1 Line - This test contains normal cluster information but no sequence.
* Test 07: No meta information - This test contains no cluster or accession information, only a sequence.
* Test 08: 2 cluster - This test contains 2 artificial clusters with sequences of different lengths to ensure that the algorithm processes them correctly.
* Test 09: 2 sequences with one protein domain - This test contains 2 very small sequences with one protein domain.
* Test 10: Protein sequence range 0-0 - This test contains one sequence with one protein domain. As an edge case, the range of the protein domain is from 0 to 0.
* Test 11: Protein sequence range 0-2 - This test is similar to test 10, but the range of the protein domain is from 0 to 2.
* Test 12: Incorrect amount of accessions - In this test, the cluster does not have the correct number of accessions.
* Test 14: Sequence with only protein domains - This test is similar to test 13, but the protein domain range is given in a different way.
* Test 15: domain length 2 - This test contains one sequence with a domain of length 2. The test is useful to check if the program is mapping the sequence correctly and if 0 to 2 results in two or three domain labels.
* Test 16: 5 real sequences with one domain - This test contains a cluster with 5 sequences of similar but different length.
* Test 17: T-Match - Five sequences, two with an F-ID without an F-group.
* Test 18: H-Match - Similar to Test 17, but sequences two and four have no T-group.
* Test 19: X-Match - Similar to Test 18, but the 5th sequence consists of F-IDs containing only the X-group.
* Test 20: Empty Sequence - When specifying the last sequence, the F-ID is completely missing. This should not happen in the input file. Therefore, it should be intercepted with an error message.
* Test 21: Large Difference in Length - This test contains 2 sequences of different lengths and shows how the implementation handles sequences with very different lengths.

## Scoring Tests
12 test files to determine and calibrate the scoring weights for the algorithm can be found in the scoring_tests folder.
* Test 01: Domain Range - Test for an alignment of 2 short sequences with special artificial domain labels.The first sequence contains a protein domain that starts and ends at the same position as the sequence. The protein domain of the second sequence starts at position 1 and ends at position 20, the end of the sequence is at position 21.
* Test 02: Disjunct Domains - This test file contains two sequences with different domains. The domain is the same length as the sequence.
* Test 03: Different Number of Zeros - The expectation for the output is that the alignment score between $AA_01 and AA_02$ is just below 0 because there are gaps but no F-ID mismatches.
* Test 04: Many Sequences with One Domain - To test for different scores. The expectation is that no different labels are matched together. If one sequence has more F-ID labels than the other, no gaps are expected, but the F-ID is aligned to a zero.
* Test 05: Long Sequence with NC Domains - Test the alignment score of 4 long sequences with one NC domain.
* Test 06: Short Sequence with NC Domains - This test is similar to test 5, but with short sequences.
* Test 07: One Sequence Times 10 - Test if the score is the same in each alignment.
* Test 08: Only 0 - Test with two sequences that do not contain a protein domain.
* Test 09: Short and Long Sequences - Same protein domains but sequences 3 and 4 have length 110 instead of 10. The expectation is, that the algorithm inserts the gaps at positions where the other string contains 0.
* Test 10: Zero Score and Self-alignment Score - To test what the difference is between the self-alignment of a sequence without a protein domain and the self-alignment of a sequence with a protein domain.
* Test 11: Many Protein Domains - The test file contains four sequences, two sequences with two protein domains each and two sequences with three protein domains each.
* Test 12: Different Number of F-IDs - The expectation for the output is that the scores of the self-alignments are better than the scores between AA\_01 and AA\_02, because there are mismatches between 0 and the F-ID. The score of the alignment AA\_01 and AA\_02 should be just below 0.
