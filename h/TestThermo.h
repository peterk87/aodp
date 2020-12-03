#ifndef __TestThermo_h__
#define __TestThermo_h__

#include "Thermo.h"

/**
 * Test driver for "nearest neighbor" model
 * 
 * Objective: test algorithm ~and NOT~ actual parameters or values in Thermo tables
 * 
 * Notes:
 *  - the Thermo Watson/Crick table (dimer matches) is taken from (Allawi and SantaLucia 1997); a newer version is provided in (SantaLucia and Hicks 2004)
 *  - Thermo parameters (initiation; A/T termination; symmetry) taken from (Allawi and SantaLucia 1997); a newer version is provided in (SantaLucia and Hicks 2004)
 * 
 * Verifies predicted values for:
 *  - DG37: perfect agreement
 *  - DH: almost perfect agreement
 *  - DS: perfect agreement for A/G and C/T mismatches; small disagreements for A/C and G/T mismatches
 *  - Tm: almost perfect agreement for A/G and C/T mismatches
 *        values ~look roughly ok~ (within +/- 1C for most entries) for A/C and G/T mismatches
 * 
 * Conclusions:
 * (1) High confidence in the algorithm itself; the newest tables and parameters should be used
 * (2) High confidence in the DG37 and DH tables and calculations
 * (3) Low confidence in the DS values for A/C and G/T mismatches: likely truncated/rounded values are reported; calculations with un-rounded values
 * (4) Low confidence in the Tm values for A/C and G/T mismatches: likely following DS errors and similar types of rounding errors
 * 
 * \see (Allawi and SantaLucia 1997), (Allawi and SantaLucia 1998a), (Allawi and SantaLucia 1998b) and (Allawi and SantaLucia 1998c)
 */
class TestThermo final : public Thermo {
public:
	TestThermo() : Thermo( Thermo::K + 37.0 ) {};
	void testDang() {
		map<string,array<double,3>> d;

		d.emplace( "AA/T" , array<double,3>{ DH.dangX( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DS.dangX( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DG37.dangX( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "AC/G" , array<double,3>{ DH.dangX( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DS.dangX( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DG37.dangX( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "AG/C" , array<double,3>{ DH.dangX( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DS.dangX( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DG37.dangX( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "AT/A" , array<double,3>{ DH.dangX( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DS.dangX( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DG37.dangX( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "CA/T" , array<double,3>{ DH.dangX( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DS.dangX( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DG37.dangX( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "CC/G" , array<double,3>{ DH.dangX( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DS.dangX( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DG37.dangX( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "CG/C" , array<double,3>{ DH.dangX( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DS.dangX( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DG37.dangX( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "CT/A" , array<double,3>{ DH.dangX( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DS.dangX( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DG37.dangX( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "GA/T" , array<double,3>{ DH.dangX( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DS.dangX( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DG37.dangX( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "GC/G" , array<double,3>{ DH.dangX( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DS.dangX( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DG37.dangX( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "GG/C" , array<double,3>{ DH.dangX( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DS.dangX( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DG37.dangX( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "GT/A" , array<double,3>{ DH.dangX( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DS.dangX( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DG37.dangX( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "TA/T" , array<double,3>{ DH.dangX( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DS.dangX( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DG37.dangX( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "TC/G" , array<double,3>{ DH.dangX( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DS.dangX( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DG37.dangX( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "TG/C" , array<double,3>{ DH.dangX( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DS.dangX( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DG37.dangX( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "TT/A" , array<double,3>{ DH.dangX( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DS.dangX( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DG37.dangX( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "A/AT" , array<double,3>{ DH.dangY( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DS.dangY( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 10.0, DG37.dangY( asc2nu.at( 'A' ), asc2nu.at( 'A' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "C/AG" , array<double,3>{ DH.dangY( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'G' )) / 10.0, DS.dangY( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'G' )) / 10.0, DG37.dangY( asc2nu.at( 'C' ), asc2nu.at( 'A' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "G/AC" , array<double,3>{ DH.dangY( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'C' )) / 10.0, DS.dangY( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'C' )) / 10.0, DG37.dangY( asc2nu.at( 'G' ), asc2nu.at( 'A' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "T/AA" , array<double,3>{ DH.dangY( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'A' )) / 10.0, DS.dangY( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'A' )) / 10.0, DG37.dangY( asc2nu.at( 'T' ), asc2nu.at( 'A' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "A/CT" , array<double,3>{ DH.dangY( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'T' )) / 10.0, DS.dangY( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'T' )) / 10.0, DG37.dangY( asc2nu.at( 'A' ), asc2nu.at( 'C' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "C/CG" , array<double,3>{ DH.dangY( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DS.dangY( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 10.0, DG37.dangY( asc2nu.at( 'C' ), asc2nu.at( 'C' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "G/CC" , array<double,3>{ DH.dangY( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'C' )) / 10.0, DS.dangY( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'C' )) / 10.0, DG37.dangY( asc2nu.at( 'G' ), asc2nu.at( 'C' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "T/CA" , array<double,3>{ DH.dangY( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'A' )) / 10.0, DS.dangY( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'A' )) / 10.0, DG37.dangY( asc2nu.at( 'T' ), asc2nu.at( 'C' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "A/GT" , array<double,3>{ DH.dangY( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'T' )) / 10.0, DS.dangY( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'T' )) / 10.0, DG37.dangY( asc2nu.at( 'A' ), asc2nu.at( 'G' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "C/GG" , array<double,3>{ DH.dangY( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'G' )) / 10.0, DS.dangY( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'G' )) / 10.0, DG37.dangY( asc2nu.at( 'C' ), asc2nu.at( 'G' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "G/GC" , array<double,3>{ DH.dangY( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DS.dangY( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 10.0, DG37.dangY( asc2nu.at( 'G' ), asc2nu.at( 'G' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "T/GA" , array<double,3>{ DH.dangY( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'A' )) / 10.0, DS.dangY( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'A' )) / 10.0, DG37.dangY( asc2nu.at( 'T' ), asc2nu.at( 'G' ), asc2nu.at( 'A' )) / 100.0 });

		d.emplace( "A/TT" , array<double,3>{ DH.dangY( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'T' )) / 10.0, DS.dangY( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'T' )) / 10.0, DG37.dangY( asc2nu.at( 'A' ), asc2nu.at( 'T' ), asc2nu.at( 'T' )) / 100.0 });
		d.emplace( "C/TG" , array<double,3>{ DH.dangY( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'G' )) / 10.0, DS.dangY( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'G' )) / 10.0, DG37.dangY( asc2nu.at( 'C' ), asc2nu.at( 'T' ), asc2nu.at( 'G' )) / 100.0 });
		d.emplace( "G/TC" , array<double,3>{ DH.dangY( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'C' )) / 10.0, DS.dangY( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'C' )) / 10.0, DG37.dangY( asc2nu.at( 'G' ), asc2nu.at( 'T' ), asc2nu.at( 'C' )) / 100.0 });
		d.emplace( "T/TA" , array<double,3>{ DH.dangY( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DS.dangY( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 10.0, DG37.dangY( asc2nu.at( 'T' ), asc2nu.at( 'T' ), asc2nu.at( 'A' )) / 100.0 });

		assert( d == nn_dang );
	};
	void testNN() {
		assert( -535 == DG37.calc(
			"CGTTGA",
// 			 |||*||
			"GCAACT"
		)); // (SantaLucia 1998; fig.1)
		assert( -832 == DG37.calc(
			"GGACTGACG",
// 			 ||||*||||
			"CCTGGCTGC"
		)); // (SantaLucia and Hicks 2004; eq.6)

// 	Adjust Thermo parameters to comply with (Allawi and SantaLucia 1997; table 1)
// NOTE: (Allawi and SantaLucia 1997; table 1) is used to calculate examples in:
//  - (Allawi and SantaLucia 1997)
//  - (Allawi and SantaLucia 1998a)
//  - (Allawi and SantaLucia 1998b)
//  - (Allawi and SantaLucia 1998c)
		DH.clear();
		DS.clear();
		DG37.clear();

		this->init( nn_parameters_correction, nn_wc_correction );

		for( const auto& e: test_seq ) {
			const string one = e.first.first;
			string two = e.second;

			const int dh_pred = int( round( - e.first.second.at( 3 ) * 10 ));
			const int ds_pred = int( round( - e.first.second.at( 5 ) * 10 ));
			const int dg_pred = int( round( - e.first.second.at( 1 ) * 100 ));

			assert( one.size() > 0 );
			assert( one.size() == two.size());

			int l = one.size();

			for( int i = 0 ; i < l ; ++i ) {
				if( two.at( i ) == ' ' ) { // replace spaces with complement from string one
					two.at( i ) = nu2asc.at( nu2compl.at( asc2nu.at( one.at( i ))));
					continue;
				}
			}

			const int dh_calc = DH.calc( one, two ) ;
			const int ds_calc = DS.calc( one, two ) ;
			const int dg_calc = DG37.calc( one, two ) ;

			assert( dh_calc == dh_pred );
			assert( ds_calc == ds_pred );
			assert( dg_calc == dg_pred );
		}
	};

	void test() {
		testNN();
		testDang();

		cout << "ok" << endl;
	};

private:
	const map<string,array<double,3>> nn_parameters_correction = { // ref: (Allawi and SantaLucia 1997; table 1)
// 	                                DH    DS   DG37
// 	                             kcal/mol  e.u. kcal/mol
		{ "init. w/term. G-C",   {  0.1,  -2.8,   0.98 }},
		{ "init. w/term. A-T",   {  2.3,   4.1,   1.03 }},

		{ "Initiation",          {  0.2,  -5.6,   1.96 }},
		{ "Terminal AT penalty", {  2.2,   6.9,   0.05 }}, // recalculated from GC, resp. AT initiation (in original)
		{ "Symmetry correction", {    0,  -1.4,   0.40  }},
	};

	const vector<pair<string,array<double,3>>> nn_wc_correction = { // ref: (Allawi and SantaLucia 1997; table 1)
// 	                   DH    DS   DG37
// 	                kcal/mol  e.u. kcal/mol
		{ "AA/TT", {  -7.9, -22.2, -1.00 }},
		{ "AT/TA", {  -7.2, -20.4, -0.88 }},
		{ "TA/AT", {  -7.2, -21.3, -0.58 }},
		{ "CA/GT", {  -8.5, -22.7, -1.45 }},
		{ "GT/CA", {  -8.4, -22.4, -1.44 }},
		{ "CT/GA", {  -7.8, -21.0, -1.28 }},
		{ "GA/CT", {  -8.2, -22.2, -1.30 }},
		{ "CG/GC", { -10.6, -27.2, -2.17 }},
		{ "GC/CG", {  -9.8, -24.4, -2.24 }},
		{ "GG/CC", {  -8.0, -19.9, -1.84 }},
	};

// 		{ "AA/TT", {  -7.6, -21.3, -1.00 }},
// 		{ "AT/TA", {  -7.2, -20.4, -0.88 }},
// 		{ "TA/AT", {  -7.2, -21.3, -0.58 }},
// 		{ "CA/GT", {  -8.5, -22.7, -1.45 }},
// 		{ "GT/CA", {  -8.4, -22.4, -1.44 }},
// 		{ "CT/GA", {  -7.8, -21.0, -1.28 }},
// 		{ "GA/CT", {  -8.2, -22.2, -1.30 }},
// 		{ "CG/GC", { -10.6, -27.2, -2.17 }},
// 		{ "GC/CG", {  -9.8, -24.4, -2.24 }},
// 		{ "GG/CC", {  -8.0, -19.9, -1.84 }}

	vector<pair<pair<string,array<double,8>>,string>> test_seq = {

// 		                               -DG37    |      -DH     |      -DS      |     Tm        | 
// 		                             kcal/mol   |    kcal/mol  |      e.u.     |     C         |
// 		      DNA duplex            expt   pred |  expt   pred |   expt   pred |  expt   pred  |  // counter
// 		   -------------------------------------+--------------+---------------+---------------+------------

// 		                            A/G mismatches (Allawi and SantaLucia 1998a; table 2)
// 		                               Molecules with Two-State Transitions
		                                                                                          {{
		     "CAAAAAAAG",        {  3.92,   3.89,  39.9,   45.3,  116.0,  133.0,  23.9,  26.3  }},		//	0
		     "    G    "                                                                        },{{
		     "CAAAGAAAG",        {  4.22,   4.33,  52.6,   47.7,  156.0,  139.7,  28.6,  28.7  }},
		     "    A    "                                                                        },{{
		     "GGACACTCG",        {  8.32,   7.86,  57.4,   51.2,  158.4,  139.6,  52.0,  51.1  }},
		     "    G    "                                                                        },{{
		     "GGACAGACG",        {  7.65,   7.39,  54.6,   56.3,  151.4,  157.4,  48.6,  47.3  }},
		     "    G    "                                                                        },{{
		     "GGACGCTCG",        {  7.49,   7.51,  54.7,   55.6,  152.2,  154.7,  47.7,  48.2  }},
		     "    A    "                                                                        },{{
		     "GGACGGACG",        {  7.38,   7.39,  52.7,   56.3,  146.2,  157.4,  47.4,  47.3  }},
		     "    A    "                                                                        },{{
		     "GGAGGCACG",        {  8.99,   8.29,  57.6,   51.4,  156.8,  138.8,  55.9,  54.0  }},
		     "    A    "                                                                        },{{
		     "CATGAAGCTAC",      {  8.46,   8.41,  71.1,   70.0,  201.9,  198.2,  49.6,  50.2  }},
		     "     G     "                                                                      },{{
		     "CATGAGGCTAC",      {  9.25,   8.61,  71.2,   66.8,  199.9,  187.3,  53.4,  51.8  }},
		     "     A     "                                                                      },{{
		     "CATGTAACTAC",      {  7.14,   6.84,  58.5,   60.5,  165.5,  172.8,  44.9,  43.4  }},
		     "     G     "                                                                      },{{
		     "GATCAATGTAC",      {  8.14,   7.72,  68.7,   65.7,  195.4,  186.5,  48.5,  47.6  }},
		     "     G     "                                                                      },{{
		     "GATCTATGTAC",      {  7.75,   7.27,  66.1,   63.7,  188.2,  181.8,  47.1,  45.2  }},
		     "     G     "                                                                      },{{
		     "GATCTGTGTAC",      {  6.85,   6.83,  61.8,   61.3,  177.1,  175.1,  43.0,  43.8  }},
		     "     A     "                                                                      },{{
		     "CCATCGCTACC",      {  9.49,   8.79,  74.5,   67.7,  209.4,  189.6,  53.7,  52.5  }},
		     "     A     "                                                                      },{{
		     "CCATTGCTACC",      {  8.18,   7.86,  67.1,   60.4,  189.8,  169.0,  49.0,  49.3  }},
		     "     A     "                                                                      },{{
		     "CCGACTCTAGCG",     { 10.15,  10.56,  64.7,   66.0,  176.0,  178.8,  60.0,  61.7  }},
		     "   G    G   "                                                                     },{{
		     "CGAGCATGATCG",     {  8.77,   8.82,  62.2,   64.2,  172.3,  178.2,  53.2,  53.6  }},
		     "   A    G   "                                                                     },{{
		     "CGCAAATTGGCG",     {  8.00,   7.80,  66.9,   59.0,  189.9,  164.8,  48.2,  49.1  }},
		     "   G    A   "                                                                     },{{
		     "CGCAAGAGACGG",     {  9.06,  10.32,  50.9,   60.4,  134.9,  161.2,  59.1,  63.3  }},
		     "   G    G   "                                                                     },{{
		     "CGTGGACCAACC",     {  7.46,   7.97,  52.0,   55.2,  143.5,  151.8,  48.0,  51.4  }},
		     "   A    G   "                                                                     },{{
		     "CTCACATGGGAG",     {  8.36,   7.56,  59.9,   56.4,  166.1,  157.4,  51.6,  47.8  }},
		     "   G    A   "                                                                     },{{
		     "CTCGACGTAGAG",     {  6.94,   6.79,  71.7,   65.8,  208.9,  190.4,  42.6,  42.1  }},
		     "   A    G   "                                                                     },{{
		     "GAGAACCTGCAG",     {  7.33,   6.93,  56.7,   51.5,  159.3,  143.0,  46.3,  46.1  }},
		     "   G    A   "                                                                     },{{
		     "GAGGACCTACAG",     {  7.83,   8.11,  48.5,   53.9,  131.0,  148.0,  51.4,  51.0  }},
		     "   A    G   "                                                                     },{{
		     "GCAACTCGGTAG",     {  8.78,   9.10,  60.7,   59.9,  167.5,  163.6,  53.7,  56.2  }},
		     "   G    A   "                                                                     },{{
		     "GCGATCTCAGCC",     {  9.37,  10.50,  57.9,   68.2,  156.5,  185.7,  58.1,  61.2  }},
		     "   G    G   "                                                                     },{{
		     "GGCAGAGAACGC",     { 10.61,  10.65,  73.5,   67.0,  202.6,  181.5,  59.3,  62.2  }},
		     "   G    G   "                                                                     },{{
// 		                               Molecules with Anomalous Two-State Transitions
		     "ATGAGCGCAT",       {  6.46,   4.72,  35.2,   45.8,   92.7,  132.2,  44.1,  31.2  }},
		     "   G  A   "                                                                       },{{
		     "ATGAGCGCAT",       {  6.69,   4.72,  46.7,   45.8,  129.0,  132.2,  43.9,  31.2  }},
		     "   G  A   "                                                                       },{{
// 		                               Molecules with Non-Two-State Transitions
		     "GGAGACACG",        {  9.30,   8.29,  61.4,   51.4,  167.8,  138.8,  56.5,  54.4  }},
// 			                                                              138.6
		     "    G    "                                                                        },{{
		     "ATGAGCTAAT",       {  4.23,   1.78,  29.3,   27.8,   80.7,   84.2,  22.4,  -1.9  }},
		     "  A    G  "                                                                       },{{
		     "ATGAGCTAAT",       {  3.43,   1.78,  36.0,   27.8,  105.0,   84.2,  18.8,  -1.9  }},
		     "  A    G  "                                                                       },{{
		     "CATGTGACTAC",      {  7.06,   6.84,  56.1,   60.5,  158.0,  172.8,  44.9,  43.4  }},
		     "     A     "                                                                      },{{
		     "CGTGTCGAAACG",     {  5.98,   7.87,  49.1,   60.0,  139.2,  167.6,  39.0,  49.6  }},
		     "   A    G   "                                                                     },{{


// 		                            C/T mismatches (Allawi and SantaLucia 1998b; table 2)
// 		                               Molecules with Two-State Transitions
		     "CAAACAAAG",        {  3.27,   3.38,  53.2,   46.0,  161.0,  137.2,  23.6,  22.7  }},		//	34
		     "    T    "                                                                        },{{
		     "CAAATAAAG",        {  3.17,   3.07,  50.0,   47.7,  151.0,  143.6,  22.2,  21.5  }},
		     "    C    "                                                                        },{{
		     "CGTCCGTCC",        {  6.70,   6.51,  58.3,   53.9,  166.4,  152.5,  42.6,  42.4  }},
		     "    T    "                                                                        },{{
		     "CGTGCCTCC",        {  6.75,   5.92,  58.1,   43.8,  165.4,  122.1,  42.9,  38.8  }},
		     "    T    "                                                                        },{{
		     "GGACCCTCG",        {  6.22,   5.77,  54.1,   46.6,  154.4,  131.5,  40.2,  37.9  }},
		     "    T    "                                                                        },{{
		     "GGACCGACG",        {  6.58,   6.51,  56.8,   53.9,  161.8,  152.5,  42.0,  42.4  }},
		     "    T    "                                                                        },{{
		     "GGAGCCACG",        {  6.58,   5.92,  58.0,   43.8,  165.7,  122.1,  41.9,  38.8  }},
		     "    T    "                                                                        },{{
		     "CACAGCAGGTC",      {  7.74,   8.15,  64.1,   62.1,  181.7,  173.8,  47.3,  50.1  }},
		     "     T     "                                                                      },{{
		     "CATGACGCTAC",      {  8.52,   7.62,  79.9,   66.2,  230.0,  188.6,  48.5,  46.8  }},
		     "     T     "                                                                      },{{
		     "CATGATGCTAC",      {  8.01,   7.31,  72.5,   67.4,  207.9,  193.4,  47.3,  45.2  }},
		     "     C     "                                                                      },{{
		     "CATGTCACTAC",      {  6.99,   6.28,  66.3,   62.0,  191.2,  179.5,  43.3,  40.3  }},
		     "     T     "                                                                      },{{
		     "CATGTTACTAC",      {  6.92,   6.28,  65.5,   62.0,  189.0,  179.5,  43.0,  40.3  }},
		     "     C     "                                                                      },{{
		     "GAACGCTGTCC",      {  8.37,   8.63,  63.7,   66.9,  178.4,  187.6,  50.7,  51.8  }},
		     "     T     "                                                                      },{{
		     "GACCTCCTGTG",      {  7.56,   7.57,  64.1,   59.0,  182.3,  165.7,  46.4,  47.5  }},
		     "     T     "                                                                      },{{
		     "GATCATTGTAC",      {  7.03,   6.51,  70.1,   64.9,  203.2,  187.9,  43.1,  41.6  }},
		     "     C     "                                                                      },{{
		     "GATCTCTGTAC",      {  6.59,   6.01,  67.6,   63.7,  196.8,  185.7,  41.3,  39.1  }},
		     "     T     "                                                                      },{{
		     "GCTAGCAATCC",      {  7.26,   7.07,  63.5,   60.4,  181.5,  171.9,  44.9,  44.4  }},
		     "     T     "                                                                      },{{
		     "CGCCAGAGCCGG",     {  6.71,   7.35,  46.8,   54.9,  129.1,  153.4,  44.0,  46.6  }},
		     "   T    T   "                                                                     },{{
		     "CGCTAGAGTCGG",     {  6.44,   7.35,  47.7,   55.4,  132.9,  155.0,  42.0,  46.5  }},
		     "   C    C   "                                                                     },{{
		     "GGCCGAGACCGC",     {  7.56,   7.77,  64.2,   58.6,  182.6,  163.8,  46.4,  48.6  }},
		     "   T    T   "                                                                     },{{
		     "GGCTGAGATCGC",     {  7.34,   8.04,  58.6,   63.4,  165.3,  178.3,  46.1,  49.3  }},
		     "   C    C   "                                                                     },{{
		     "CGACCATATGTTCG",   {  6.33,   6.58,  52.1,   64.2,  147.5,  185.9,  40.9,  41.2  }},
		     "   T      C   "                                                                    },{{
		     "CGTCTCATGATACG",   {  7.25,   7.84,  74.3,   78.4,  216.0,  227.4,  43.7,  45.9  }},
		     "   T      C   "                                                                    },{{
		     "CTCCACATGTTGAG",   {  6.78,   6.72,  67.0,   72.4,  194.2,  211.6,  42.2,  41.8  }},
		     "   T      C   "                                                                    },{{
		     "CTCTCATATGCGAG",   {  6.51,   6.00,  65.7,   68.8,  190.9,  202.3,  41.0,  38.7  }},
		     "   C      T   "                                                                    },{{
		     "CAACTTGATATTAATA", {  9.70,  10.13,  98.4,   99.4,  286.0,  287.4,  50.2,  52.0  }},
		     "        C       "                                                                 },{{
// 		                                Molecules with Non-Two-State Transitions
		     "CGAGCGTCC",        {  6.29,   6.35,  58.4,   50.2,  168.0,  141.2,  40.3,  41.6  }},
		     "    T    "                                                                        },{{
		     "GAACGCAGTCC",      {  7.12,   8.44,  74.2,   64.0,  216.2,   179.0, 43.2,  40.6  }},
// 		incorrect:                          6.32           62.0            179.3      data entry error in (Allawi and SantaLucia 1998b): switched lines
		     "     T     "                                                                      },{{
		     "GATCTTTGTAC",      {  6.56,   6.32,  32.3,   62.0,   83.1,   179.3, 45.7,  51.2  }}, 
// 		incorrect:                          8.44           64.0            179.0
		     "     C     "                                                                      },{{
		     "CTCTATGGTACTGC",   {  7.48,   7.76,  78.6,   74.2,  229.3,  214.0,  44.3,  46.3  }},
		     "   C      T   "                                                                   },{{
		     "GCATCTGCGGCTAG",   { 10.03,   9.87,  42.9,   75.6,  105.9,  211.8,  72.1,  55.4  }},
		     "   C      T   "                                                                   },{{

// 		                   A/C mismatches (Allawi and SantaLucia 1998c; table 2)
// NOTE: entries with pH 5.0 are ignored
// WARNING: there are small, systematic errors in predicted DS
//          These are assumed to be calculation errors in the reference.
//          the calculated value is entered in the test array; the original value (reference) is included in comments below each line
		     "CAAAAAAAG",        {  2.92,   2.56,   35.8,  37.8,  106.0,  113.3,  14.9,  14.1  }},		//	65
		     "    C    "                                                                        },{{
		     "CAAACAAAG",        {  3.08,   3.08,   40.3,  39.0,  120.0,  115.5,  18.2,  17.5  }},
// 			                                                              115.9
		     "    A    "                                                                        },{{
		     "CGAGCGTCC",        {  6.98,   6.15,   54.1,  49.8,  152.1,  140.6,  44.6,  40.5  }},
// 			                                                              140.5
		     "    A    "                                                                        },{{
		     "CGTCCGTCC",        {  6.58,   5.99,   57.8,  49.1,  165.2,  138.8,  41.9,  39.2  }},
// 			                                                              138.9
		     "    A    "                                                                        },{{
		     "CGTGCCTCC",        {  7.07,   6.24,   53.5,  46.8,  149.7,  130.6,  45.2,  41.4  }},
// 			                                                              130.5
		     "    A    "                                                                        },{{
		     "GGACCCTCG",        {  6.06,   5.77,   49.1,  45.2,  138.8,  126.9,  39.4,  37.9  }},
// 			                                                              127.0
		     "    A    "                                                                        },{{
		     "GGACCGACG",        {  6.08,   5.99,   50.3,  49.1,  142.6,  138.8,  39.5,  39.2  }},
// 			                                                              138.9
		     "    A    "                                                                        },{{
		     "GGAGCCACG",        {  6.96,   6.24,   56.0,  46.8,  158.2,  130.6,  44.3,  41.4  }},
// 			                                                              130.5
		     "    A    "                                                                        },{{
		     "CATGAAGCTAC",      {  7.24,   6.99,   67.8,  62.5,  195.4,  178.7,  44.3,  43.9  }},
// 			                                                              178.8
		     "     C     "                                                                      },{{
		     "CCACACCAGAG",      {  8.49,   7.95,   70.5,  55.0,  199.9,  151.4,  49.9,  50.6  }},
// 			                                                              151.6
		     "     A     "                                                                      },{{
		     "CCACACGAGAG",      {  7.89,   7.86,   65.8,  58.0,  186.7,  161.4,  47.8,  49.2  }},
// 			                                                              161.6
		     "     A     "                                                                      },{{
		     "GAGAGCACACC",      {  8.14,   7.98,   66.0,  62.4,  186.4,  175.2,  49.0,  49.2  }},
// 			                                                              175.3
		     "     A     "                                                                      },{{
		     "GATCAATGTAC",      {  6.56,   6.23,   66.3,  56.8,  192.7,  162.7,  41.2,  40.3  }},
// 			                                                              162.9
		     "     C     "                                                                      },{{
		     "GATCTATGTAC",      {  6.11,   6.02,   61.6,  55.0,  178.8,  157.6,  39.2,  38.8  }},
// 			                                                              158.0
		     "     C     "                                                                      },{{
		     "GTACAAAGATC",      {  5.85,   5.50,   59.5,  53.8,  173.1,  155.4,  37.9,  36.6  }},
		     "     C     "                                                                      },{{
		     "GTAGCATCATG",      {  7.77,   7.14,   68.1,  58.2,  194.5,  164.4,  46.8,  45.1  }},
// 			                                                              164.6
		     "     C     "                                                                      },{{
		     "GTAGTAACATG",      {  6.33,   5.75,   60.5,  53.2,  174.7,  152.7,  40.4,  37.6  }},
// 			                                                              152.9
		     "     C     "                                                                      },{{
		     "GTAGTCACATG",      {  6.21,   5.75,   55.5,  53.2,  159.0,  152.7,  40.0,  37.6  }},
// 			                                                              152.9
		     "     A     "                                                                      },{{
		     "CCGACTCTAGCG",     {  6.75,   7.33,   52.1,  54.1,  146.2,  150.7,  43.5,  46.6  }},
// 			                                                              150.9
		     "   C    C   "                                                                     },{{
		     "CGCAAATTCGCG",     {  5.37,   5.18,   48.6,  44.6,  139.4,  127.2,  35.0,  33.4  }},
		     "   C    A   "                                                                     },{{
		     "CGCAAGAGACGG",     {  6.30,   6.96,   50.4,  48.6,  142.1,  134.2,  40.9,  45.7  }},
// 			                                                              134.1
		     "   C    C   "                                                                     },{{
		     "GCGCTCTCCGCC",     {  6.98,   7.52,   50.6,  58.1,  140.5,  163.0,  45.2,  47.3  }},
		     "   A    A   "                                                                     },{{
		     "GGCCGAGACCGC",     {  7.16,   7.29,   51.8,  49.2,  143.9,  135.0,  46.2,  47.2  }},
// 			                                                              135.3
		     "   A    A   "                                                                     },{{
		     "CGACCATATGATCG",   {  6.28,   6.66,   39.1,  55.0,  105.8,  155.7,  41.8,  42.2  }},
// 			                                                              156.1
		     "   A      C   "                                                                   },{{
		     "CGTCGAGGACAACC",   {  7.20,   7.93,   42.6,  56.8,  114.3,  157.4,  48.5,  50.1  }},
		     "   A      C   "                                                                   },{{
		     "CGTCTCATGAAACG",   {  5.92,   6.82,   61.3,  58.6,  178.6,  166.8,  38.2,  43.4  }},
		     "   A      C   "                                                                   },{{
		     "CTCACATATGCGAG",  {  7.62,    5.60,   77.1,  68.0,  224.0,  201.1,  45.0,  37.1  }},
// 			                                                              200.9
		     "   C      A   "                                                                   },{{
		     "CTCCACATGTAGAG",  {  5.97,    6.04,   59.1,  64.8,  171.2,  189.2,  38.5,  38.2  }},
// 			                                                              189.8
		     "   A      C   "                                                                   },{{
		     "GAGAAGCGGTCCAG",  {  8.94,    8.04,   59.1,  51.1,  161.6,  138.6,  55.1,  52.5  }},
		     "   C      A   "                                                                   },{{
		     "GCAACTCCGGCTAG",  {  9.41,    9.32,   66.9,  72.5,  185.3,  203.6,  55.4,  53.9  }},
// 			                                                              203.4
		     "   C      A   "                                                                   },{{
		     "CAACTTGATATTAATA",{  8.71,    9.53,   90.9,  87.9,  265.0,  252.2,  47.7,  51.6  }},
// 			                                                              252.4
		     "       C        "                                                                 },{{

// 		                            G/T mismatches (Allawi and SantaLucia 1997; table 3)
// WARNING: there are small, systematic errors in predicted DS and, in some cases, DH
//          These are assumed to be calculation errors in the reference.
//          the calculated value is entered in the test array; the original value (reference) is included in comments below each line

// 		                               Molecules with Two-State Transitions
		     "CAAAGAAAG",       {  4.21,    3.72,   46.7,  48.0,  137.0,  142.5,  27.6,  24.5  }},		//	96
// 			                                               46.1           136.6
		     "GTTTTTTTC"                                                                       },{{
		     "CAAATAAAG",       {  4.43,    4.27,   55.6,  50.3,  165.0,  148.1,  30.2,  28.5  }},
// 			                                               48.4           142.2
		     "GTTTGTTTC"                                                                       },{{
		     "CGTGTCTCC",       {  7.77,    8.03,   52.0,  52.4,  142.7,  142.9,  50.0,  52.1  }},
// 			                                                              142.8
		     "GCACGGAGG"                                                                       },{{
		     "CGAGTGTCC",       {  8.47,    8.43,   62.5,  59.5,  174.1,  164.5,  51.6,  51.8  }},
// 			                                                              164.8
		     "GCTCGCAGG"                                                                       },{{
		     "GGACTCTCG",       {  7.47,    7.61,   57.9,  50.5,  162.7,  138.1,  46.9,  49.1  }},
// 			                                                              138.4
		     "CCTGGGAGC"                                                                       },{{
		     "GGACTGACG",       {  8.37,    8.32,   62.6,  58.5,  174.9,  161.6,  51.0,  50.9  }},
// 			                                                              162.2
		     "CCTGGCTGC"                                                                       },{{
		     "GGAGTCACG",       {  8.82,    8.03,   65.4,  52.4,  182.4,  142.9,  52.6,  52.1  }},
// 			                                                              142.8
		     "CCTCGGTGC"                                                                       },{{
		     "CATGAGGCTAC",     {  8.61,    8.27,   69.9,  67.2,  197.6,  189.8,  50.6,  50.1  }},
// 			                                                              189.6
		     "GTACTTCGATG"                                                                     },{{
		     "CATGTGACTAC",     {  7.72,    7.23,   64.2,  65.6,  182.1,  187.9,  47.2,  44.8  }},
// 			                                                              188.0
		     "GTACATTGATG"                                                                     },{{
		     "CCATCGCTACC",     {  10.42,   9.71,   79.8,  71.6,  223.8,  199.4,  56.6,  55.7  }},
		     "GGTAGTGATGG"                                                                     },{{
		     "CCATTGCTACC",     {  9.06,    8.51,   75.7,  67.3,  214.8,  189.4,  51.5,  51.0  }},
// 			                                                              189.3
		     "GGTAATGATGG"                                                                     },{{
		     "GATCATTGTAC",     {  7.81,    7.10,   69.3,  65.9,  198.1,  189.3,  46.9,  43.8  }},
// 			                                                              189.6
		     "CTAGTGACATG"                                                                     },{{
		     "GATCTTTGTAC",     {  7.31,    6.66,   67.6,  64.0,  194.5,  184.6,  44.6,  41.8  }},
// 			                                                              184.9
		     "CTAGAGACATG"                                                                     },{{
		     "GTAGCGTCATG",     {  9.80,    9.06,   76.6,  72.0,  215.5,  202.7,  54.7,  52.8  }},
// 			                                                              202.6
		     "CATCGTAGTAC"                                                                     },{{
		     "GTAGTGACATG",     {  7.88,    7.23,   68.3,  65.6,  194.7,  187.9,  47.3,  44.8  }},
// 			                                                              188.0
		     "CATCATTGTAC"                                                                     },{{
		     "CCATGCGTAACG",    {  8.98,    9.31,   71.2,  70.0,  200.7,  195.6,  52.1,  54.1  }},
// 			                                               72.0
		     "GGTATGCGTTGC"                                                                     },{{
		     "CGAGACGTTTCG",    {  6.98,    7.53,   61.0,  65.4,  174.1,  186.6,  43.7,  45.7  }},
// 			                                                              186.8
		     "   T    G   "                                                                     },{{
		     "CGAGCATGTTCG",    {  7.23,    8.12,   59.8,  68.4,  169.6,  194.4,  45.3,  48.7  }},
// 			                                                              194.2
		     "   T    G   "                                                                     },{{
		     "CGCGAATTTGCG",    {  9.78,    9.60,   79.3,  74.4,  224.1,  209.0,  53.9,  54.5  }},
// 			                                                              208.8
		     "   T    G   "                                                                     },{{
		     "CGTGACGTTACG",    {  8.13,    8.37,   73.3,  68.0,  210.0,  192.2,  47.7,  49.3  }},
// 			                                                              192.6
		     "   T    G   "                                                                     },{{
		     "CGTGTCGATACG",    {  8.39,    8.63,   73.8,  70.0,  210.9,  197.8,  48.8,  49.9  }},
// 			                                                              198.4
		     "   T    G   "                                                                     },{{
		     "CGTTACGTGACG",    {  7.87,    8.37,   64.9,  68.0,  184.0,  192.2,  47.8,  49.3  }},
// 			                                                              192.6
		     "   G    T   "                                                                     },{{
		     "CTCGGATCTGAG",    {  8.43,    7.86,   75.0,  69.2,  214.5,  197.6,  48.8,  46.2  }},
// 			                                                              198.4
		     "   T    G   "                                                                     },{{
		     "CTCTCATGGGAG",    {  6.47,    7.06,   50.4,  55.0,  141.5,  154.4,  41.9,  45.3  }},
		     "   G    T   "                                                                     },{{
		     "CTCTGATCGGAG",    {  7.52,    7.86,   60.3,  69.2,  170.3,  197.6,  46.7,  46.2  }},
// 			                                                              198.4
		     "   G    T   "                                                                     },{{
		     "CTGTCATGGCAG",    {  8.30,    7.90,   59.2,  58.8,  164.2,  164.0,  51.4,  50.8  }},
// 			                                                              163.2
		     "   G    T   "                                                                     },{{
		     "CTGTGATCGCAG",    {  8.65,    8.70,   67.0,  73.0,  188.2,  207.2,  51.4,  50.6  }},
		     "   G    T   "                                                                     },{{
		     "CTTGGATCTAAG",    {  5.91,    5.46,   64.0,  60.6,  187.3,  177.6,  38.1,  35.2  }},
// 			                                                              178.2
		     "   T    G   "                                                                     },{{
		     "CAACTTGATATTAATA",{  9.40,    9.77,   91.3,  98.5,  264.0,  285.7,  47.1,  50.1  }},
// 			                                                              286.4
		     "GTTGAATTATAATTAT"                                                                 },{{
		     "CAACTTGATATTAATA",{ 10.10,   10.97,   92.6, 102.8,  266.0,  295.7,  49.4,  53.4  }},
// 			                                                              296.5
		     "GTTGAGCTATAATTAT"                                                                 },{{
		     "CAACTTGATATTAATA",{ 10.50,   10.94,   95.5, 100.1,  274.0,  287.1,  50.5,  54.0  }},
// 			                                                              287.7
		     "GTTGAACTATAGTTAT"                                                                 },{{

// 		                             Molecules with Marginally Non-Two-State NNParametersdynamics
		     "CGTCTGTCC",       {  8.04,    8.32,   56.5,  58.5,  156.5,  161.6,  50.1,  50.9  }},
// 			                                                              162.2
		     "GCAGGCAGG"                                                                        },{{
		     "GATCTGTGTAC",     {  7.81,    7.21,   70.6,  66.3,  202.5,  190.2,  46.7,  44.1  }},
// 			                                                              190.7
		     "CTAGATACATG"                                                                      },{{
		     "CGAGTCGATTCG",    {  7.65,    7.79,   69.6,  67.4,  199.7,  192.2,  46.1,  46.4  }},
// 			                                                              192.6
		     "   T    G   "                                                                     },{{
		     "CTTGCATGTAAG",    {  6.14,    6.30,   56.7,  64.4,  162.9,  187.2,  39.5,  40.5  }},
// 			                                                              187.0
		     "   T    G   "                                                                     },{{
		     "CGTGTCTAGATACG",  {  9.22,    9.60,   78.3,  82.2,  222.6,  233.9,  51.7,  52.1  }},
// 			                                                              234.4
		     "   T      G   "                                                                   },{{

// 		                                Molecules with Non-Two-State Transitions
		     "CGTTGCGTAACG",    {  7.89,   9.35,   58.7,   73.2,  163.1, 205.8,  49.2,  57.5  }},
// 		 incorrect:                       10.33            73.3          203.4                likely incorrect row in original
		     "    T  G    "                                                                     },{{
		     "CTCGCATGTGAG",    {  8.19,    8.70,   64.6,  73.0,  181.9,  207.2,  49.7,  50.6  }},
		     "   T    G   "                                                                     },{{
		     "CTGGCATGTCAG",    {  8.54,    7.90,   67.8,  58.8,  191.2,  164.0,  50.7,  50.8  }},
// 			                                                              163.2
		     "   T    G   "                                                                     },{{
		     "CTGGGATCTCAG",    {  7.66,    7.06,   73.8,  55.0,  213.4,  154.4,  45.6,  45.3  }},
		     "   T    G   "                                                                     },{{
		     "GCGTACGCATGCG",   { 12.65,   15.16,   80.9,  97.3,  220.0,  264.7,  66.3,  71.0  }},		//	136
// 			                                                              264.4
		     "CGCATGTGTACGC"                                                                     },

			 
	};
};

#endif

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
