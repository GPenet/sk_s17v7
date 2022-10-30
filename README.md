# sk_s17v7
search 17/18 clues in a solution grid V7

This is an attempt to improve efficiency of the code, and an extension to the 18 search in the same model.
The version V6 is written for the search with the distribution 665 in bands and stacks.
V7 covers all the field for the 17 and the 18 clues search
pass1 with bands 3 having more than 6 clues
pass2 with not more than 6 clues in bands and stacks.
in the 17 clues search, pass2 is split in 2a,2b to avoid 5 clues in band 3


V7 performs significantly better than the previous version. 
The main changes in V7  have been partially implemented in the 18 search, we have

More UAs/GUAs produced at the start using a revised process to build them.
Potential valid bands 1+2 are searched trough a direct expansion of these UAs
Action on potential valid bands 1+2 has been improved compared to the V6

<uas first harvest>
known stack uas are used (internal table of band uas)
brute force variants have been used to speed up the process and extend the search ot more cases

<UA band 1+2 expansion> 
this is done in several steps with UAs clearing and restrucuting at the start of each step.
The last step procuces a table of potential valid bands 1+2 with common first 6clues (and more)

<gangster uas filter>
The use of gangster Uas has been improved 
The use of UAs 3 bands not gangster UAs has been widely improved, with a start benefiting from the known stack uas.

<use of fresh UAs>
As in V6, all fresh UAs are use as much as possible (in the limit of a reasonable size).
At the end, in the final brute force check, at least one false must be in band 3

<assigned cells in band 3>
The process in band 3 has been adjusted to benefit as much as possible from UAs reduced to one cell applying known constraints. 

detailed comments can be found at  https://gpenet.pagesperso-orange.fr/P17/p17_18.htm

