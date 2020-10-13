::

    .=============================================.
                                             _ 
    ._ _  ___ ._ _  ___  ___  ___ ._ _  ___ | |
    | ' |<_> || ' |/ . \| . \<_> || ' |/ ._>| |
    |_|_|<___||_|_|\___/|  _/<___||_|_|\___.|_| 2
                        |_|                    
    .=============================================.
    

Introduction
============

Nanopanel2 (np2) is a somatic variant caller for Nanopore panel sequencing data.
Np2 works directly on FAST5 files and outputs VCF v4.2 and TSV files containing 
variant calls. It also produces haplotype map PDFs that inform about the haplotypes
of called (PASS) variants.


Installation
============

The recommended way to use np2 is via the released `Singularity`_ v3.6.1 image which contains all required 3rd party tools in the supported software versions.

.. _Singularity: https://sylabs.io/docs/

Users that prefer to run the python code directly should make sure that the following
3rd party tools are available

3rd Party tools
---------------
The following 3rd party executables are called by the np2 python pipeline and are packaged in the singularity container:

* `samtools`_ v1.9
* `porechop`_ v0.2.4
* `minimap2`_ v2.17
* `ngmlr`_ v0.2.7
* `last`_ v1042
* bgzip v1.2.1++

Note that the call path of the respective tools can be configured in np2's JSON config file 
(section 'exe') which makes it possible to install these tools locally.

.. _samtools: https://github.com/samtools/samtools
.. _porechop: https://github.com/rrwick/Porechop
.. _minimap2: https://github.com/lh3/minimap2/
.. _ngmlr: https://github.com/philres/ngmlr
.. _last: http://last.cbrc.jp/


General usage
=============

::

   singularity run nanopanel2_XXX.sif call --conf config.json --out .

We recommend to run np2 with 128gb of RAM and 8 cores (configure via the 'threads' parameter in the JSON config file). 

Configuration file
------------------

Np2 is configured via a single JSON configuration file, `a commented example`_ can be found in the docs folder.
Please note that np2 uses the `_commentjson`_ package to parse input JSON files, so you can use Python/JavaScript style inline comments.
 
.. include:: docs/config.json.example

.. _`a commented example`: docs/config.json.example
.. _commentjson: https://github.com/vaidik/commentjson

