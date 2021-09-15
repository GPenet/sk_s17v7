# sk_s17v7
search 17 clues in a solution grid V7

This is an attempt to improve efficiency of the code.
The current version is V6 written for the search with the distribution 665 in bands and stacks.
V7 has the same limitation (this is the code needed to finish the scan in progress)
V6 performs significantly better than the previous version. 
The main changes in V7 currently in drafting mode have been partially implemented in the 18 search, we have

a main external loop band 1 band 2 changed to have in band 1 a re shaping of the still valid Uas
a switch to 128 uas in the main central loop 
the addition of stacks uas 
a siginficant change in the organization of the Gangster uas
an attempt to optimize the use of the gangster uas socket 2
a change in the guess sequence of the brute force (both bands 1+2 and final check for the whoe set of clus)
a change in the sequence of the next controls after the main loop testing the 128 first uas
  here the main goal is to reduce the calls to the brute force.
  the brute force for bands 1+2 is always done first
  the final check (if any) uses the fact that the bands 1+2 subset is valid to optimize the guess sequence
  
The guess sequence uses 2 facts to optimize the process
  the final true solution is known, the guess uses as much as possible the true value first
  the final check will nearly always fail and produce a us (usually no 17 clues)
  note : the guess "true" first increases the chance to get a small additionnal uas
  
And everywhere, the maximm is done to benefit of the fresh uas produced in he brute force check
