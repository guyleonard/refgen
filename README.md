REFGEN - REFormat GEne Names for Phylogenies

Leonard, G., et al. (2009). REFGEN and TREENAMER: Automated Sequence Data Handling for Phylogenetic Analysis in the Genomic Era. Evolutionary Bioinformatics Online, 5:1-4.

Install
-------

1) Copy refgen.html to your website directory e.g. /var/www/
2) Create a directory under /var/www called 'refgen' - make sure it is owned by the 'www-data' user and has 755 user permissions.
3) Place the 'refgen_usage.txt' and 'doe_jgi_list.csv' in the 'refgen' folder.
4) Put the 'refgen_1_beta.pl' file in to your /var/lib/cgi-bin directory


refgen_1_beta.pl - This is the current version running at http://gna-phylo.nhm.ac.uk/refgen.html - Much of the code has been completely rewritten from the publication version and new features added, unfortunately I did not track changes! From now on I will. Yes, this PERL script does have a bunch of HTML and JAVASCRIPT in it. Deal with it ;)

refgen_1_beta_5.pl - This should be moved to an experimental branch, read header comments for experiment.

rfs.cgi - Publication ready version of REFGEN - may not completely work with refgen.html in this repository.


REFGEN is released under a GPLv3 license.
