CLdb (CRISPR Loci Database) 
===========================

CLdb is toolset for organizing and analyzing large amounts of CRISPR data.

Existing webtools are easy but suffer from some limitatios:

* They can require A LOT of tedious clicking for large datasets.
* They are not very flexible.
 * Adding the tools to existing analysis workflows can be challenging.
 * For instance, it is usually hard to incorporate such tools into an IPython
notebook unless the webtool has a good API.
* They often limit tranparency and reproducibility of the analysis.
 * Lack of tranparency and reproducibility are a major issue in bioinformatics.


I'm developing this toolset in an attempt to remedy some of these issues.
This project is currently under heavy development and may change
dramatically in the future.


## Major Features

* Easy spacer Blasting
 * Filter out spacer blast hits to other CRISPR arrays
 * Get the protospacer of each blast hit
  * This includes the adjacent PAM region
 * Get the crRNA DNA template (crDNA) for each blast query
 * Make protospacer-crDNA alignments
 * Get summaries on protospacer-crDNA mismatches for the SEED sequence and entire protospacer
 * Get the PAM regions for each hit
* Make detailed comparative plots of CRISPR systems 
 * The plots can include information on:
  * CAS gene conservation among CRISPRs
  * Spacer conservation among CRISPRs
  * Location of the leader region
* Summarize your dataset quickly
 * Get the number of spacers shared among:
  * CRISPR loci
  * CRISPR subtypes
  * taxa
 * Make repeat consensus sequences
  * Use for making weblogos or trees



## Minor Features

* Organize and query subsets of your CRISPR dataset
 * Select by subtype, taxa, or individual CRISPR loci
* Make gff3 files of the CRISPR features



## INSTALLATION

	git clone https://github.com/nyoungb2/CLdb.git
	

You can either:

1. add the CLdb/bin/ to your $PATH
2. or install the perl/R scripts and packages

#### Adding CLdb/bin to your $PATH

See [this link](http://kb.iu.edu/data/acar.html) for more info.

#### Installing the scripts and packages

     cd CLdb/
     perl Build.PL
     ./Build
     ./Build test
     sudo ./Build install

You may need to use `Build install --install_base=/path/that/I/can/write/to/`
if you you lack root/superuser privilages.


## Documentation

See the [wiki](https://github.com/nyoungb2/CLdb/wiki).

## ISSUES/CONTACT

I'd love to here your feedback. 

Please provide it via [Issues](https://github.com/nyoungb2/CLdb/issues) on GitHub.


## LICENSE AND COPYRIGHT

Copyright (C) 2013 Nick Youngblut

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.