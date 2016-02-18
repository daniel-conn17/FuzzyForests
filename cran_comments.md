## Test environments
* local OS X install, R version 3.2.3 (2014-10-31)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 
There was one NOTE.
"* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Daniel Conn <djconn17@gmail.com>'
New submission"

Response to feedback by Professor Ripley:

"high-dimensional (with a hypthen)."
   A hyphen has been added.  
   
"Given a reference for this 'new' algorithm."
    A reference is given in inst/CITATION.
    References are also supplied in the documentation
    for appropriate functions.
  
"This runs its examples for far too long.  I stopped at

auk% cat *.timings
name    user    system  elapsed
WGCNA_control   0       0       0       
ff       31.657   0.591 150.787 
screen_control  0.001   0.000   0.000   
select_control  0.001   0.000   0.001
The policies require a few seconds per example, not a few minutes.  I also 
noticed multiple R processes running for minutes."

Some code from the examples is currently not checked because I 
believe that they take too long to run.
These examples should take between 5-10 seconds.
Now multiple R processes should not be running. 


