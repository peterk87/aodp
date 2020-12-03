#define __ApplicationTest_cpp__

#include "Application.h"

#include "Alignment.h"

/**
 * Test driver for features
 * 
 * call like so: assert( testDriver());
 * 
 * TODO: eventually, move to its own class
 * TODO: eventually, move components to their own methods
 */
bool Application::testDriver() {
// 	Reverse complement
	const vector<pair<string,string>> rc_data = {
		{ "ACGT", "ACGT" },
		{ "TCATCCTTTTCAGGTTGACCTC", "GAGGTCAACCTGAAAAGGATGA" },

// 	WARNING: generated values
		{"TVMARVCV", "BGBYTKBA"},
		{"VCVCGSGG", "CCSCGBGB"},
		{"GCSMCCRR", "YYGGKSGC"},
		{"SACMGRRG", "CYYCKGTS"},
		{"TRGRSSSR", "YSSSYCYA"},
		{"MRRTMMGS", "SCKKAYYK"},
		{"CCCASMSA", "TSKSTGGG"},
		{"VARMSRAC", "GTYSKYTB"},
		{"RASGVRVM", "KBYBCSTY"},
		{"MATVSTCV", "BGASBATK"},
		{"AGVVSVRG", "CYBSBBCT"},
		{"TCTCARSS", "SSYTGAGA"},
		{"RVATVVCA", "TGBBATBY"},
		{"CVGGVTMM", "KKABCCBG"},
		{"VRRAVAMA", "TKTBTYYB"},
		{"GAARVCGV", "BCGBYTTC"},
		{"TRCARSTMMMRMCAGR", "YCTGKYKKKASYTGYA"},
		{"SCRMTVVGMSGRMVCV", "BGBKYCSKCBBAKYGS"},
		{"CVARCAGGSMMMSGVV", "BBCSKKKSCCTGYTBG"},
		{"MVGAGCMVCAMSAVTA", "TABTSKTGBKGCTCBK"},
		{"GSGASCVRSTAMSTVC", "GBASKTASYBGSTCSC"},
		{"AVSGGCTGTVGVMSCV", "BGSKBCBACAGCCSBT"},
		{"GCATSVVTSASTTCGA", "TCGAASTSABBSATGC"},
		{"GTRVSTAAGCCSMSTR", "YASKSGGCTTASBYAC"},
		{"ASAMTRSRATVSASMV", "BKSTSBATYSYAKTST"},
		{"RTSTSVVVVVGTGATV", "BATCACBBBBBSASAY"},
		{"GVVGMTGSATTSCMVV", "BBKGSAATSCAKCBBC"},
		{"RGVTGRMCRVMCAAVG", "CBTTGKBYGKYCABCY"},
		{"VSARGAGMSGTRTVAG", "CTBAYACSKCTCYTSB"},
		{"GATVTMVTGVGTCMAT", "ATKGACBCABKABATC"},
		{"SVVSVGVGTCAMAMGC", "GCKTKTGACBCBSBBS"},
		{"AMRCVSCMSGAMTRVC", "GBYAKTCSKGSBGYKT"},
		{"MSARTMSRSVAMVMVCSVSGRAMVVASRGTCS", "SGACYSTBBKTYCSBSGBKBKTBSYSKAYTSK"},
		{"RTRCVGVTRVTTAGGSCGSMCATTARSRTGGS", "SCCAYSYTAATGKSCGSCCTAABYABCBGYAY"},
		{"SMMGTRVMTCTTATRVVGVSSSVGTTATTRTV", "BAYAATAACBSSSBCBBYATAAGAKBYACKKS"},
		{"MTTGSMASRVCSSVSSGMSAMCTMGMATSSMM", "KKSSATKCKAGKTSKCSSBSSGBYSTKSCAAK"},
		{"VGACCMVCVMMTTMTGCTATSTTRGRVTSRGM", "KCYSABYCYAASATAGCAKAAKKBGBKGGTCB"},
		{"RRCTCTAATSVCRSTATRTATMSMVSVSRAMR", "YKTYSBSBKSKATAYATASYGBSATTAGAGYY"},
		{"VSRTTVAGCASCMGSRGVCMMTAGRSTVCRCS", "SGYGBASYCTAKKGBCYSCKGSTGCTBAAYSB"},
		{"SRMAGSTCVGGGTGRVVSARAVCVVMACRVGS", "SCBYGTKBBGBTYTSBBYCACCCBGASCTKYS"},
		{"VMTVTMRTCMAVGGCTSTCCMRGMSATRVSRR", "YYSBYATSKCYKGGASAGCCBTKGAYKABAKB"},
		{"RCGSVCTASVGMTTCRMMGMASARTGRSATTM", "KAATSYCAYTSTKCKKYGAAKCBSTAGBSCGY"},
		{"RMCCVRGVASTGMMMVGAGTCSVACCVVTRRT", "AYYABBGGTBSGACTCBKKKCASTBCYBGGKY"},
		{"MCMVSGVSARTTGCTRRARTRVGCRCCSCCRC", "GYGGSGGYGCBYAYTYYAGCAAYTSBCSBKGK"},
		{"ATRGAMVMRVSRVRCCAMVTCMSVSARSAVAT", "ATBTSYTSBSKGABKTGGYBYSBYKBKTCYAT"},
		{"MSCRSMMGGCMVVGMCASTCMMSACMRGRSRM", "KYSYCYKGTSKKGASTGKCBBKGCCKKSYGSK"},
		{"GRCSMAAACMGGTRSAMGATTSVMVVVCCGST", "ASCGGBBBKBSAATCKTSYACCKGTTTKSGYC"},
		{"SMTMVGSAMMMTAMASMRVVMVTATVTAMATC", "GATKTABATABKBBYKSTKTAKKKTSCBKAKS"},
	};

	for( const auto& rc2: rc_data ) {
		assert( convertNu2Asc( convertNu2ReverseComplement( convertAsc2Nu( rc2.first ))) == rc2.second );
	}

// 	Cover flip
	assert(
		Cover<Position>( { 10000, 10000 }, {{ 10100, 100 }, { 10300, 300 }}).flip() ==
		Cover<Position>( { 10000, 10000 }, {{ 19400, 300 }, { 19800, 100 }}));

// 	Folding
	const vector<pair<string,array<double,5>>> fold_data = {
// // 		   sequence         DG    DH     DS    Tm    error
		{ "CGCAAAGCG",   { -0.9, -20.4, -62.9, 51.1, 0.05 }},  // calculated
// 		   (((...)))
// DINAMelt (2016-10-28)   -0.9  -20.4  -62.8  51.5

		{ "GCGAAACGC",   { -0.9, -20.4, -62.9, 51.1, 0.05 }},  // calculated
// DINAMelt (2016-10-28)   -0.9  -20.4  -62.8  51.5

		{ "GTAGAAACTAC", { 0.2,  -23.4, -76.0, 34.8, 0.05 }},  // calculated
// DINAMelt (2016-10-28)   0.2   -23.4  -76.1  34.4

		{ "ACGCAAAGCGT", { -2.3, -26.6, -78.4, 66.1, 0.05 }},  // calculated
// DINAMelt (2016-11-01)   -2.3, -26.6, -78.3, 66.4

		{ "ATGAAACAT",   { 1.2,  -13.5, -47.5, 11.1, 0.05 }},  // calculated
// 		                   1.2   -13.5  -47.5  11.3
// 
		{ "TAGAAACTA",   { 1.7,  -12.8, -46.7,  1.0, 0.05 }},  // calculated
// DINAMelt (2016-11-01)   1.1   -13.2  -46.0   13.5    TAGAAACTA
//                                                      .((...)).
// WARNING: discrepancy with DINAMelt http://unafold.rna.albany.edu/results2/twostate-fold/161030/135905/index.php

		{ "GCAAAATGC",   {-0.1,  -16.1, -51.5, 39.5, 0.05 }},  // calculated
// DINAMelt (2016-11-01)  -0.1   -16.1  -51.5  39.7

		{ "ACGCAAAGCGA", {-2.8,  -28.8, -84.0, 69.7, 0.05 }},  // calculated
// DINAMelt (2016-10-28)  -2.8   -28.8  -83.9  70.3

		{ "AACGCAAAGCGAA",{ -2.8, -28.8, -84.0, 69.7, 0.05 }}, // calculated
// DINAMelt (2016-10-28)    -2.8  -28.8  -83.9  70.3

		{ "CGGAATGCG",   {  -0.0, -19.6, -63.1, 37.2, 0.05 }}, // calculated
// DINAMelt (2016-11-01)     0.1  -15.9  -51.7  34.5

// WARNING: large discrepancies:
		{ "TCGGAATGCGT", { -0.9, -28.0, -87.2, 47.8, 0.05 }},  // calculated
// DINAMelt (2016-11-01)   -0.8  -24.3  -75.8  47.6
		
		{ "GCTCGGAATGCGTGTTCGGTACTG", { -0.9, -28.0, -87.2, 47.8, 0.05 }},  // calculated
// DINAMelt (2016-11-01)                -0.8  -24.3  -75.8  47.6

	};
	const Thermo th(
		Thermo::T37C,   // Tx = 37C
		0.01,           // cT = 0.01 mM
		1.0             // [Na+] = 1M
	);

	for( const auto& fo: fold_data ) {
		const string& s = fo.first;
		const auto& pa = fo.second;
		const double error = pa.at( 4 );

		const size_t l = s.size();
		deque<Length> d( l, l );
		const string nu = convertAsc2Nu( s );

		Fold f(
			nu,  // nucleotide string
			0,   // starting position in string
			l,   // length of string
			d,   // length deque
			l,   // min
			l,   // max (=min=length) for folding
			th,  // thermodynamic parameters
			onull // output
// 			cout
		);

		double dg, dh, ds, tm;
		
		f.fold();
		assert( f.foldParameters( 0, l-1, dg, dh, ds, tm ));

		assert( abs( dg - pa.at( 0 )) <= error );
		assert( abs( dh - pa.at( 1 )) <= error );
		assert( abs( ds - pa.at( 2 )) <= error );
		assert( abs( tm - pa.at( 3 )) <= error );
	};

// 	Align
	struct AlignData {
		string s1;
		string s2;
		Alignment::Element align;
	};
	const vector<struct AlignData> align_data = {
	{"AAAGTCC"
	,"GTCCAA"
	,{2,4,4}},


	{"GCATGCT"
	,"GATTACA"
	,{0,4,8}},


	{"CCAAAAGG"
	,"AATAA"
	,{3,4,5}},


	{"GAAAAAAT"
	,"GAAT",
	{2,3,4}},


	{"AGACTAGTTAC"
	,"CGAGACGT",
	{2,5,6}},

	{"CAAC"
	,"ATA"
	,{1,2,3}},


	{"CAGAAC"
	,"ATAA"
	,{2,3,4}},

	{"CCAAAAGG"
	,"AATAA"
	,{3,4,5}},

	{"TTTATTCCTGCGGGCGCATTTTCCCACACCATAAAAACTTTCCACGTGAACCGTTACAATTATGTTCTGTGCTCTCTCTCGGGAGGGCTGAACGAAGGTGGCGCTATATGTAAAAGTGTAGCGTTGCCGATGTACTTTTAAACCCATTACACTAATACTGAACTATACTCCGAGAACGAAAGTTTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTCGCACATCGATGAAGAACGCTGCGAACTGCGATACGTAATGCGAATTGCAGAATTCAGTGAGTCATCGAAATTTTGAACGCACATTGCACTTTCGGGATATTCCTGGAAGTATGCTTGTATCAGTGTCCGTACATCAAACTTGCCTTTCTTTTTTTGTGTAGTCAAGGAGAGAAATGGCAGAATGTGAGGTGTCTCGCTGACTCCCTCTTCGGAGGAGAAGACGCGAGTCCCTTTAAATGTACGTTCGCTCTTTCTTGTGTCTAAGATGAAGTGTGACTTTCGAACGCAGTGATCTGTTTAGATCGCTTTGCGCGAGTGGGCGACTTCGGTTAGAACATTAAAGGAAGCAACCTCTATTGGCGGTATGTTAGGCTTCGGCCCGACTTTGCAGCTGACAGTGTGTTGTTTTCTGTTCTTTCCTTGAGGTGTACCTGTCTTGTGTGAGGCAATGGTCTGGGCAAATGGTTATTGTGTAGTAGATTGTTGCTGCGCTTGGGCGCCCTACTTTATTGTGGGGTAAAGAAGGCAACACCAATTTGGGACTAGTCTGTGGGGATTTATTCCTGCGGGCGCATTTTC"
	,"CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCACTTAGTTGGGGGCTAGTCCCGGCGGCTGGCTGTCGATGTCAAAGTTGACGGCTGCTGCTGTGTGTCGGGCCCTATCATGGCGAGCGTTTGGGTCCCTCTCGGGGGAACTGAGCCAGTAGCCCTTTTCTTTTAAACCCATTCTTGAATACTGAATATACTGCGGGGACGAAAGTCTCTGCTTTTAACTAGATAG"
	,{54,148,242}},

	{"TTTATTCCTGCGGGCGCATTT"
	,"CCACACCTAAAAAACTTT"
	,{-4,8,20}},


	{"CCT"
	,"CCA"
	,{1,2,3}},

	{"CCACACCT"
	,"CCTGC"
	,{1,3,5}},

	{"CCACACCT"
	,"ACAC"
	,{4,4,4}},
	{"CATTACTGAGTTTATGCTC"
	 "TCACGAGCTAACCTCCCACCCGTGTCTATTACATCTTGTTGCTTCGGTGCGCAGCCCGCGGAGGTTTACCTCTAAAGGTCACGTGCCGAGGACGCCATTT"
	 "GAACTCTGTATTATATTGCAGTCTGAGAATATAACTTAATTAGTTAAAACTTTCAACAACGGATCTCTTGG"
	 "TTCCGGTATCGATGAAGAACGCAGCGAAATGCGATAAATAATGTGAATTGCAGAATTCAGTGAATCATCGAGTCTTTGAACGCACATTGCGCCCC"
	 "CTGGTATTCCGGGGGGCATGCCTGTCCGAGCGTCATTGCTGCCCTCAAGCCCGGCTTGTGTGTTGGGTCCTCGTCCCTCCGGGGACAGGCCCGAA"
	 "AGGCAATGGCAGTACCGCGTCCGGTCCTCGAGCGTATGGGGCTTTGTCACCCGCTCTGTAGGCCCGGCCGGCGCTCCGCCGACCAACCAAAAACT"
	 "ATTTTTCAGGTTGACCTCGGATCAGGTAGGG"
	,"TCACGAGCTAACCTCCCACCCGTGTCTATTACATCTTGTTGCTTCGGTGCGCAGCCCGCGGAGGTTTACCTCTAAAGGTCACGTGCCGAGGACGCCATTT"
	,{100,100,100}}  // from issue-102
	};

	Alignment a;
	for( const auto& a2: align_data ) {
		const string s1 = convertAsc2Nu( a2.s1 );
		const string s2 = convertAsc2Nu( a2.s2 );

		assert(( s1.size() < Alignment::max_length_smallest_sequence ) || ( s2.size() < Alignment::max_length_smallest_sequence ));
		auto al = a.align( s1, 0, s1.size(), s2, 0, s2.size());

		assert( al.S == a2.align.S );
		assert( al.M == a2.align.M );
		assert( al.L == a2.align.L );

		assert( al.L > 0 );
	}

	return true;
}

// This file is part of aodp (the Automated Oligonucleotide Design Pipeline)
// 
// (C)	HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014-2018)
// (C)	Manuel Zahariev mz@alumni.sfu.ca (2000-2008,2014-2018)
// 
// aodp is free software: you can redistribute it and/or
// modify it under the terms of version 3 of the GNU General Public
// License as published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License (version 3) for more details.
// 
// You should have received a copy of the GNU General Public License
// (version 3) along with this program. If not, see
// http://www.gnu.org/licenses/.
