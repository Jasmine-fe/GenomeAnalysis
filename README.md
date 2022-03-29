# RepeatFi - "Repeat" identification based on "F"ragment "i"nteger Interval

### Introduction
A large proportion of the DNA sequences in many organisms are repeat sequences. The excessively repeat sequences in DNA are also closely related to diseases. However, the data generated by gene sequencing is very large. Finding repeated sequences in a large amount of sequence data will consume a lot of computer resources, including CPU power and memory. Therefore, we use CutMat method to cut the sequence. The concept of CutMat is to cut a sequence into fragments, then compare the length of the interval between fragments instead of the traditional word-by-word comparison method. This is an efficient method to find the repeated sequence in the genome.

---

### Repeat finding flowchart
<img width="527" alt="Screen Shot 2022-02-23 at 5 10 58 PM" src="https://user-images.githubusercontent.com/36653195/155289660-09ce9e12-220d-4b06-bbe7-ebf6c1fa559d.png">

**Choose cutter and Cutter a sequence to fragments by cutter**
  - Choose Cutter 1, 2, ..., n
  - Cut sequence into fragments each cutter, If there are 3 cutters, the tool will generate 3 different outcomes.
  - Group fragments by length for each cutter outcome. 
  
**Filter fragments**
  - Filtered out non-identical fragments and single fragments in the fragments’ length groups**
  - Generate repeat candidate table, including length, sequence and position information.

**Finding Repeat fragments**
  - For each cutter, annotate repeat fragments in sequence, mark 1 in repeat position, else mark 0, we called the marked result as SeqState.
  - If using multiple cutters, generate mergeState by sum up of all SeqStates which generated by each cutter.
  - Repeat identification
    single cutter: choose state == 1 in the SeqState 
    multiple cutter(N): choose state == N in the mergeState

**Generate output**
  - Generate Repeat Fragment Output.
