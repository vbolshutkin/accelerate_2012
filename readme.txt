Terms and conditions :
http://software.intel.com/fr-fr/articles/AYC-early2012_tcs/

The goal of the contest is to optimize and parallelize
a reference algorithm and code. The criteria : execution time
(various single memory machines/inputs/conditions)
Code, example input and expected output are provided.


Code is the reference, but to make it easier for you
to understand the code, here is an explanation :

Problem : Compare 1 reference sequence (in a file)
with a variable number of sequences (from 1 or more files).
Sequences are composed from 4 characters : A,T,C,G.

Comparaison : Find substrings of >= N characters
matching exactly.

Parameters :
  -  maximum number of worker threads to use.
  -  minimum matching size.
  -  reference sequence file (one sequence - one file)
  -  following arguments : sequence files to compare
     (one or more sequences per file)
Example : ./run 24 16 refseq.txt input.txt
See output.txt for the expected output.

Important : You can write your own code from scratch, but
result has to be exactly the same as this reference code
for any set of input-output-parameters.
More input-output files will be posted on the contest page
later to help you test your code.

Testing and benchmarking :
If you are registrered as a contestant, you can upload
a zip file containing your code and Makefile
(directly at the root of the zip file)
on our testing and benchmarking cluster.
http://intel-software-academic-program.com/ayc-upload/
A large input file is available from :
http://intel-software-academic-program.com/contests/ayc/early2012/test_input_1.tar.bz2

Advices : Think parallelism on a lot of cores,
think potential for execution time acceleration,
use efficient data structures,
not only algorithmic tricks.

Tips : For top performance, check vectorization
and STTNI instructions in SSE4.2
In order to process big files, you need to be
less greedy than our sample program.
In order to do that you can avoid storing non matching pairs,
by using for example hash tables instead of arrays.
