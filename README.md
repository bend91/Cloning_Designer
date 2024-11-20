## Cloning_Designer
Identifying and designing the optimal strategy for cloning DNA plasmids

Current strategies are optimised for NEBs HiFi cloning, but will work for any cloning that can be done with a 15bp overlap

Have implemented in C++ and python, the C++ has been compiled and the distributable is at CPP/bin/main.

Everything is command line at the moment but all works, just provide the sequences you want to clone together in fasta format, in the order you want them to be cloned, and the tool will provide optimal primers and a suggested annealing temperature.

Currently this only works for putting together the internal fragments, I haven't implemented a strategy for adding in a specific vector yet but that is likely the next thing to do.

# Example fasta file
>ScfV
NNNNNNNNNNNNNN
>Fc
NNNNNNNNNNNN....
>2A_RFP
NNNNNNN....


# Example output:
>ScFv
- Forward:
- Reverse: NNNNNNNNNNNNNNN

>Stalk
- Forward: NNNNNNNNNNNNNNN
- Reverse: NNNNNNNNNNNNNNN

>TM
- Forward: NNNNNNNNNNNNNNN
- Reverse: NNNNNNNNNNNNNNN

>Costim
- Forward: NNNNNNNNNNNNNNN
- Reverse: NNNNNNNNNNNNNNN

>ITAM
- Forward: NNNNNNNNNNNNNNN
- Reverse:

Optimal annealing temperature is: 68.8577oC


Obviously the Ns will be replaced by the actual primer sequences.


To run a test to see how it looks you can use the file in test_data, for the C++ just type test when it asks for the file name, for the python
