#include "Application.h"

const string Application::_help =
"NAME\n"
"     aodp - Automated Oligonucleotide Design Pipeline version \"2.5.0.1\"\n"
"\n"
"SYNOPSIS\n"
"    aodp [*options*] *output* *fasta-sequence-file...*\n"
"\n"
"DESCRIPTION\n"
"    \"aodp\" generates oligonucleotide signatures for sequences in FASTA\n"
"    format and for all groups in a phylogeny in the Newick tree format.\n"
"\n"
"INPUT\n"
"    Multiple input sequence files in the FASTA format can be specified\n"
"    anywhere on the command line. The usual file name wildcards can be used.\n"
"    Any command line argument without the prefix \"--\" will be treated as a\n"
"    file name.\n"
"\n"
"    Each FASTA file can contain multiple sequences. Sequence identifiers are\n"
"    read from the FASTA description lines (lines starting with \"\">).\n"
"    Everything on the description line following a space is ignored.\n"
"\n"
"    Multiple occurrences of sequences with the same identifier are\n"
"    considered as the same sequence, with several discountiguous sections.\n"
"\n"
"    If any of the specified *fasta-sequence-file*-s cannot be open or\n"
"    completely parsed, aodp will terminate with an error.\n"
"\n"
"OUTPUT\n"
"    At least one output option is required. For each output option, if\n"
"    file-name is not specified, the result is written to standard output.\n"
"\n"
"    --strings[=file-name]\n"
"        Writes all oligo strings grouped by sequence or group\n"
"\n"
"    --positions[=file-name]\n"
"        Output for Array Designer (tab-separated list of oligo sites)\n"
"\n"
"    --ranges[=file-name]\n"
"        All ranges of sites of oligos grouped by sequence or group\n"
"\n"
"    --fasta[=file-name]\n"
"        All oligos in FASTA format\n"
"\n"
"    --gff[=file-name]\n"
"        All oligos in GFF format\n"
"\n"
"    --newick[=file-name]\n"
"        Writes a phylogeny (Newick tree format), with exactly the same\n"
"        structure as the input tree, with generated labels for internal\n"
"        nodes.\n"
"\n"
"        Requires a --tree-file input\n"
"\n"
"    --node-list[=file-name]\n"
"        List of sequences for every internal node in the phylogeny\n"
"\n"
"        Requires a --tree-file input\n"
"\n"
"    --lineage[=file-name]\n"
"        Lineage of every sequence in the phylogeny\n"
"\n"
"        Requires a --tree-file input\n"
"\n"
"    --fold[=file-name]\n"
"        Predicted secondary structure and calculated two-state melting\n"
"        temperature for all oligonucleotide signature candidates discarded\n"
"        because of higher melting temperature than --max-melting\n"
"\n"
"        If the option --max-melting is not specified, will print the\n"
"        predicted secondary structure for all oligonucleotide candidates.\n"
"        Melting temperatures below 0C or over 100C will be indicated as \"*\".\n"
"\n"
"    --cladogram=(file-name)\n"
"        printout in the eps format of a cladogram associated with the\n"
"        annotated phylogeny tree. All nodes with identified oligos are\n"
"        marked with a \"*\" and these nodes and their descendents are coloured\n"
"        red.\n"
"\n"
"        Requires a --tree-file option and a file-name\n"
"\n"
"    --time[=file-name]\n"
"        Tab-separated \"user\", \"system\", \"elapsed\" time and maximum memory\n"
"        usage for various phases of processing\n"
"\n"
"    --basename=(name)\n"
"        All oligos in the following formats:\n"
"\n"
"        * \"name.oligo.strings\" strings grouped by sequence and/or group\n"
"\n"
"        * \"name.oligo.fasta\" FASTA\n"
"\n"
"        * \"name.oligo.gff\" GFF\n"
"\n"
"        * \"name.oligo.tab\" tab-separated\n"
"\n"
"        * \"name.oligo.positions\" positions\n"
"\n"
"        * \"name.oligo.ranges\" ranges\n"
"\n"
"        If the option --tree-file is specified (phylogeny tree), the\n"
"        following output files will also be generated:\n"
"\n"
"        * \"name.newick\" phylogeny tree with labeled internal nodes\n"
"\n"
"        * \"name.node-list\" list of sequences for every internal node\n"
"\n"
"        * \"name.lineage\" lineage of every sequence in the phylogeny\n"
"\n"
"        * \"name.cladogram.eps\" cladogram associated with the annotated\n"
"          phylogeny tree\n"
"\n"
"        The output options --strings, --fasta, --gff, --tab, --positions,\n"
"        --ranges --newick, --node-list, --lineage and --cladogram are\n"
"        incompatible with --basename.\n"
"\n"
"    --cluster-list[=file-name]\n"
"        Generates the list of clusters: groupings of sequences that can be\n"
"        identified by at least one oligonucleotide signature (cluster\n"
"        oligonucleotide signature).\n"
"\n"
"        The output contains the following tab-separated columns:\n"
"\n"
"        * Numeric identifier of the cluster. This is a generated value.\n"
"\n"
"        * Space-separated list of all identifiers of sequences contained in\n"
"          the cluster\n"
"\n"
"        * If the cluster matches exactly a target node (leaf node or\n"
"          internal phylogeny node), an additional column with the name of\n"
"          this node is included\n"
"\n"
"    --cluster-oligos[=file-name]\n"
"        Generates the list of all oligonucleotide signatures for all\n"
"        clusters identified. The output contains the following tab-separated\n"
"        columns:\n"
"\n"
"        * Numeric identifier of the cluster. This matches the value in\n"
"          --cluster-list.\n"
"\n"
"        * String representation of the cluster oligonucleotide signature\n"
"\n"
"    --clusters=(name)\n"
"        Generates a file of cluster nodes (\"name.cluster-list\") and a file\n"
"        of cluster oligonucleotide signatures (\"name.cluster-oligos\")\n"
"\n"
"        The output options --cluster-list and --cluster-oligos are\n"
"        incompatible with --clusters.\n"
"\n"
"    --match=(target-FASTA-file)\n"
"        Finds the closest matching source sequence for each sequence from\n"
"        the target-FASTA-file. Works as follows:\n"
"\n"
"        * Build the smallest set of sequences (\"minimum-set\") that can\n"
"          explain the largest portion of each target sequence based on\n"
"          matching cluster oligonucleotide signatures.\n"
"\n"
"        * Align (using a modified Needleman-Wunsch global alignment\n"
"          algorithm) each sequence in the \"minimum-set\" with the target\n"
"          sequence and calculate the percentage similarity. Multiple source\n"
"          sequences may have the same percentage similarity with the target\n"
"          sequence.\n"
"\n"
"        For each source sequence with maximum percentage similarity to the\n"
"        target sequence, prints to standard output or to the file specified\n"
"        by --match-output the following tab-delimited fields:\n"
"\n"
"        * Target sequence identifier\n"
"\n"
"        * Source sequence identifier or - if no matching source sequence was\n"
"          found\n"
"\n"
"        * Percentage similarity\n"
"\n"
"        * Length (bp) of the portion of the target sequence that matches\n"
"          perfectly the source sequence\n"
"\n"
"        * Length of the target sequence\n"
"\n"
"        * Size of the \"minimum-set\" (affects the running time)\n"
"\n"
"        * Size of the \"maximum-set\" of sequences contained in any clusters\n"
"          matched by the target sequence\n"
"\n"
"        Suggestion: --max-homolo=0 must be specified, otherwise matches\n"
"        containing homopolymers will be ignored and the match percentage\n"
"        will be lowered\n"
"\n"
"    --match-output=(output-file)\n"
"        Redirect the output from --match to the output file\n"
"\n"
"        Requires a --match input\n"
"\n"
"OPTIONS\n"
"    Other command line parameters are optional.\n"
"\n"
"    --help\n"
"        Display this command summary then exit.\n"
"\n"
"    --version\n"
"        Display version and copyright information\n"
"\n"
"    --oligo-size=(size[-size])\n"
"        Look for oligonucleotide signatures of sizes within the specified\n"
"        range\n"
"\n"
"    --tree-file=(relative-file-name)\n"
"        Use the phylogeny file in the Newick tree format that groups the\n"
"        sequences in the input; oligos will also be sought for all nodes in\n"
"        the phylogeny\n"
"\n"
"    --outgroup-file=(relative-file-name)\n"
"        The outgroup list is a case-sensitive list of species (sequence\n"
"        names) that are to be excluded from the final output. They will\n"
"        still be used in the generation of oligos, but oligos specific to\n"
"        them will be omitted.\n"
"\n"
"        Will terminate with an error if a sequence name in the\n"
"        --outgroup-file is not found in any *fasta-sequence-file*.\n"
"\n"
"    --isolation-file=(relative-file-name)\n"
"        A list of taxa or sequences to isolate (one item per line).\n"
"        Sequences whose name match (case-sensitive) as a complete substring\n"
"        any of the items in the --isolation-file will be the targets of the\n"
"        oligo search. Only oligos for sequences that match items on the list\n"
"        or nodes that have sequences that match items on the list will be\n"
"        sought. Individual entries in the --isolation-file may match more\n"
"        than one sequence. For example, \"carotovorum\" will match the\n"
"        following sequences:\n"
"\n"
"        * \"Pectobacterium_carotovorum_actinidae\"\n"
"\n"
"        * \"Pectobacterium_carotovorum_brasiliense\"\n"
"\n"
"        Will terminate with an error if an item in the --isolation-file does\n"
"        not match any sequence in any *fasta-sequence-file*.\n"
"\n"
"    --database=(file)\n"
"        Validate the resulting oligos against a reference database in the\n"
"        FASTA format. This option requires specifying a corresponding\n"
"        --taxonomy option.\n"
"\n"
"        Will terminate with an error if the database file contains sequence\n"
"        names that are not specified in the taxonomy file.\n"
"\n"
"        Requires a --taxonomy input.\n"
"\n"
"    --taxonomy=(file)\n"
"        Taxonomy information associated with the sequences in the reference\n"
"        database. Each line has tab-separated sequence name and lineage,\n"
"        ending in species name (\"s__Genus_species_...\").\n"
"\n"
"        The species name (Genus species) of all input sequences encoded as\n"
"        \"XX_9999_Genus_species_...\" will be extracted. Oligos for input\n"
"        sequences that don't match any sequence or match only sequences in\n"
"        the reference database for the same species will be kept\n"
"        (super-specific oligos); sequence oligos that match reference\n"
"        sequences with no species name or with another species name will be\n"
"        discarded.\n"
"\n"
"        Oligos for internal nodes (group oligos) that match reference\n"
"        sequences with a species name other than any of the species names in\n"
"        the group (from the input sequences) will also be discarded.\n"
"\n"
"        Requires a --database input.\n"
"\n"
"    --ignore-SNP\n"
"        Single polymorphism sites (SNP) will be ignored in the design of\n"
"        oligos. --ignore-SNP will generate less oligos in more time\n"
"\n"
"    --reverse-complement\n"
"        Will take into account the reverse complement of all sequences\n"
"        (default=no):\n"
"\n"
"        * SO will also be generated for the reverse complement of each\n"
"          sequence\n"
"\n"
"        * All generated SO will be validated against all direct and reverse\n"
"          complement of all sequences\n"
"\n"
"    --crowded=(yes/no)\n"
"        Indicates whether for the --ranges and --positions outputs an oligo\n"
"        range is populated with intermediary sites (default \"no\")\n"
"\n"
"              |<--range with signature oligos-->|\n"
"      \n"
"                        middle of range\n"
"         =====|================|================|=====  nucleotide sequence\n"
"\n"
"              *                *                *       --crowded=no\n"
"              * *    *    *    *    *    *    * *       --crowded=yes>\n"
"               |                      |\n"
"         --first-site-gap    --inter-site-gap\n"
"\n"
"    --first-site-gap=(gap-size)\n"
"        For a --crowded display, the size of the gap between the border of\n"
"        the range and the first interior site (default 5)\n"
"\n"
"    --inter-site-gap=(gap-size)\n"
"        For a --crowded display, indicates the size of the gap between sites\n"
"        inside an oligo range (default 5).\n"
"\n"
"        This parameter cannot be set to 0.\n"
"\n"
"    --ambiguous-sources=(yes/no)\n"
"        Whether sequences containing ambiguities are considered in the\n"
"        analysis. The names of these sequences will be written to a file\n"
"        called excluded.fasta in the current directory; default \"yes\"\n"
"\n"
"    --threads=(count)\n"
"        The maximum number of threads for multiprocessor systems. By\n"
"        default, aodp will detect the number of available processors/cores\n"
"        (\"n\") and will use \"n-1\" threads, or one thread on single processor\n"
"        systems.\n"
"\n"
"    --max-ambiguities=(count)\n"
"        Indicates the maximum number of ambiguous bases (default 5).\n"
"        Sequences with more than this number of ambiguous bases will not be\n"
"        included in further processing. Their names will be written to a\n"
"        file called excluded.fasta in the current directory. If this\n"
"        parameter is not specified, no sequences will be filtered out based\n"
"        on the number of their ambiguous bases.\n"
"\n"
"    --max-crowded-ambiguities=(count)\n"
"        Indicates the maximum number of ambiguous bases within an oligo\n"
"        size. Sequences with more than this number of ambiguous bases\n"
"        anywhere within a window will not be included in further pro\n"
"        cessing. Their names will be written to a file called excluded.fasta\n"
"        in the current directory. If this parameter is not specified, no\n"
"        sequences will be filtered out based on the number of their\n"
"        ambiguous bases.\n"
"\n"
"    --ambiguous-oligos=(yes/no)\n"
"        Whether oligos containing ambiguous bases will be sought in the\n"
"        analysis; default \"no\"\n"
"\n"
"    --max-homolo=(size)\n"
"        Maximum length of a homopolymer (e.g. only \"A\"s) in any oligo;\n"
"        default 4; 0 means no oligos with homopolymers will be removed\n"
"\n"
"    --max-melting=(temperature-C)\n"
"        Maximum melting temperature (Celsius) for any discovered oligo.\n"
"        Oligos with higher melting temperature will be removed from the\n"
"        result. If this option is not specified, all oligos are reported.\n"
"\n"
"        The two-state melting temperature is calculated using the NN model\n"
"        (SantaLucia and Hicks 2004), applied to the most stable\n"
"        single-strand self-folding configuration at temperature\n"
"        --max-melting.\n"
"\n"
"          ----------------------------------------------\n"
"         | Tm = DH x ( DS + R x ln ( CT / x )) + 273.15 |\n"
"          ----------------------------------------------\n"
"               (SantaLucia and Hicks 2004; eq. 3)\n"
"\n"
"        Use --salt to specify the salt (NaCl) concentration and --strand to\n"
"        specify the strand concentration (CT). In melting temperature\n"
"        calculations for oligonucleotide signatures, x is always 1.\n"
"\n"
"        Since a given strand may have a number of competing self-folding\n"
"        configurations, the actual melting temperature will be consistently\n"
"        lower for multi-state coupled equlibria than the calculated\n"
"        two-state melting temperature.\n"
"\n"
"        The options --max-melting and --ambiguous-oligos=yes are\n"
"        incompatible.\n"
"\n"
"    --salt=(Na+ concentration in M)\n"
"        Na+ concentration (default \"1M\"). Valid values are between \"0.05M\"\n"
"        and \"1.1M\"\n"
"\n"
"        The Na+ concentration is used in entropy calculations:\n"
"\n"
"          -----------------------------------------------\n"
"         | DS[Na+] = DS[1M NaCl] x 0.368 x N/2 x ln[Na+] |\n"
"          -----------------------------------------------\n"
"               (SantaLucia and Hicks 2004; eq. 5)\n"
"\n"
"        \"N\" is the total number of phosphates in the folded configuration.\n"
"        For self-folding configurations, \"N\" is the strand length in\n"
"        nucleotides minus 1.\n"
"\n"
"    --strand=(single strand concentration in \"mM\")\n"
"        strand concentration (default 0.1) in mM used in (SantaLucia and\n"
"        Hicks 2004; eq. 3). Valid values are between 0.01 and 100.\n"
"\n"
"        In melting temperature calculations for oligonucleotide signatures,\n"
"        x is ALWAYS 1.\n"
"\n"
"EXAMPLES\n"
"    By default, output from \"aodp\" is directed to the standard console:\n"
"\n"
"                 output type\n"
"                  =======\n"
"         $ aodp --strings Anisogramma.fasta\n"
"                          =================\n"
"                           input database\n"
"\n"
"    Wildcards are supported for the files in the input database:\n"
"\n"
"         $ aodp --strings fasta/*\n"
"\n"
"    By specifying a phylogeny tree in the Newick format, oligonucleotide\n"
"    signatures will also be sought for internal nodes of the phylogeny:\n"
"\n"
"         $ aodp --strings --tree-file=Anisogramma.tre Anisogramma.fasta\n"
"                                      ===============\n"
"                           phylogeny tree for input database\n"
"\n"
"    A number of output types (in separate files with the same prefix) will\n"
"    be automatically generated by using the option \"--basename\":\n"
"\n"
"                           output file prefix\n"
"                           ==================\n"
"         $ aodp --basename=Anisogramma.output Anisogramma.fasta\n"
"\n"
"    The target oligonucleotide signature length (in base pairs) can be\n"
"    specified as a number:\n"
"\n"
"         $ aodp --strings --oligo-size=32 Anisogramma.fasta\n"
"\n"
"    ...or as a range:\n"
"\n"
"         $ aodp --strings --oligo-size=24-32 Anisogramma.fasta\n"
"\n"
"    Oligonucleotide signatures containing ambiguous base pairs can be\n"
"    sought:\n"
"\n"
"         $ aodp --strings --ambiguous-oligos=yes Anisogramma.fasta\n"
"\n"
"    Sequences with ambiguities in the input database can be filtered out\n"
"    using the options --ambiguous-sources, --max-ambiguities and\n"
"    --max-crowded-ambiguities:\n"
"\n"
"         $ aodp --strings --max-crowded-ambiguities=8 Anisogramma.fasta\n"
"\n"
"    Specific sequences or groups of sequences can be excluded from the\n"
"    output using an *outgroup* file:\n"
"\n"
"         $ aodp --strings --outgroup-file=A.outgroup Anisogramma.fasta\n"
"\n"
"    The output can be restricted to only sequences specified within an\n"
"    *isolation* file:\n"
"\n"
"         $ aodp --strings --isolation-file=A.iso Anisogramma.fasta\n"
"\n"
"    The output can be filtered for maximum melting temperature (Celsius):\n"
"\n"
"         $ aodp --strings --oligo-size=48 --max-melting=45 Anisogramma.fasta\n"
"\n"
"    Additional parametes (strand concentration in mM and salt concentration\n"
"    in M) for the calculation of the melting temperature can be included:\n"
"\n"
"         $ aodp --strings --oligo-size=48 --max-melting=45 Anisogramma.fasta\n"
"                --salt=0.1 --strand=2\n"
"\n"
"    Actual melting temperature and folding configuration for the resulting\n"
"    oligonucleotide signatures can be displayed:\n"
"\n"
"         $ aodp --strings --oligo-size=48 --max-melting=45 Anisogramma.fasta\n"
"                --fold --salt=0.1 --strand=2\n"
"\n"
"    The output of \"aodp\" can be cross-referenced against a reference\n"
"    database, by using the --taxonomy and --database options:\n"
"\n"
"         $ aodp --tree-file=Anisogramma.tre Anisogramma.fasta\n"
"                --database=UNITE_public_mothur_full_10.09.2014\n"
"                --taxonomy=UNITE_public_mothur_full_10.09.2014_taxonomy.txt\n"
"\n"
"SEE ALSO\n"
"    clado\n"
"      generates an eps file from a phylogeny tree file in a Newick format.\n"
"      Nodes that have names ending in \"*\" are colored red.\n"
"\n"
"    clus\n"
"      interprets the results of matches of experimental samples against\n"
"      cluster oligo signatures generated by \"aodp\"\n"
"\n"
"    DNA thermodynamics\n"
"      SantaLucia J Jr, Hicks D. 2004. The Thermodynamics of DNA Structural\n"
"      Motifs. Annu. Rev. Biophys. 33:415-40\n"
"\n"
"COPYRIGHT\n"
"    This file is part of \"aodp\" (the Automated Oligonucleotide Design\n"
"    Pipeline)\n"
"\n"
"     (C) HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014-2018)\n"
"     (C) Manuel Zahariev mz@alumni.sfu.ca (2000-2008,2014-2018)\n"
"\n"
"    \"aodp\" is free software: you can redistribute it and/or modify it under\n"
"    the terms of version 3 of the GNU General Public License as published by\n"
"    the Free Software Foundation.\n"
"\n"
"    This program is distributed in the hope that it will be useful, but\n"
"    WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General\n"
"    Public License (version 3) for more details.\n"
"\n"
"    You should have received a copy of the GNU General Public License\n"
"    (version 3) along with this program. If not, see\n"
"    http://www.gnu.org/licenses/.\n"
"\n"
"AUTHORS\n"
"    Manuel Zahariev (mz@alumni.sfu.ca) Wen Chen (wen.chen@agr.gc.ca)\n"
"\n"
"\n";
