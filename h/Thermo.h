#ifndef __Thermo_h__
#define __Thermo_h__

#include <iostream>

#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>

#include <cmath>
#include <cassert>

#include "Types.h"
#include "array2.h"
#include "util.h"
#include "NNParameters.h"

using namespace std;

class Thermo {
public:
	Thermo( const double _Tx, const double _cT, const double _Na ) :
		Tx ( _Tx ),
		cT( _cT ),
		rho( Thermo::R * log( cT / Thermo::x / 1000 )), // 1000 converts from mM to M
		Na( _Na ),
		lambda( 0.184 * log( Na )),
		_lambda_rho( int( round( lambda / rho * 1000 )))
	{
		this->init( nn_parameters, nn_wc );
	};

	const static double R; // gas constant (cal/K/mol)
	const static double K; // 0C (K)
	const static double T37C; // 37C (K)
	const static double x; // strand concentration divisor (SantaLucia and Hicks 2004; eq.3); always 1 for oligonucleotides

	const double Tx;     // experimental temperature (K)
	const double cT;     // strand concentration (mM)
	const double rho;    // pseudo-melting temperature divisor
	const double Na;     // [Na+] concentration (M)
	const double lambda; // salt adjusment factor 0.184 * ln( [Na+] )

	const int _lambda_rho;   // salt adjusment factor (1000x)

	NNParameters DH;   // enthalpy
	NNParameters DS;   // entropy
	NNParameters DG37; // free energy at 37C
	NNParameters DG;   // free energy at experimental temperature
	/**
	 * Pseudo-melting temperature in K; 1M NaCl; 1000x multiplier (Tp)
	 * 
	 * Maximizing this value is equivalent to minimizing the free energy at experimental temperature
	 * 
	 * This value is ON THE SAME SIDE of the experimental temperature (Tx) as the melting temperature (Tm)
	 * 
	 *  +---------------------------------------+
	 *  | Tx - Tp = ( 1 + alpha ) * ( Tx - Tm ) | (1)
	 *  +---------------------------------------+
	 * 
	 *  +----------------------+
	 *  | Tx > Tp <=> Tx > Tm  | (2)
	 *  +----------------------+
	 * 
	 * where:
	 *  - Tx = experimental temperature; design target: this needs to be higher than the melting temperature
	 *  - Tp = pseudo-melting temperature (this value)
	 *  - Tm = actual melting temperature
	 *  - alpha = DS / r
	 *  - rho   = R ln ( cT / x )
	 *  - DS = entropy (from tables)
	 *  - R  = gas constant
	 *  - cT = strand concentration (mM)
	 *  - x  = 4 for non-self complementary duplexes / 1 for self-complementary duplexes (always 1 for oligonucleotides)
	 * 
	 *  +---------------------------------------+
	 *  |  Tp = DG0 / r = ( DH - Tx * DS ) / r  | (3)
	 *  +---------------------------------------+
	 * 
	 * where:
	 *  - DG0 = free energy at experimental temperature (Tp)
	 *  - DH  = entalpy
	 * 
	 * NOTE:
	 *  +-------------------+
	 *  | ( 1 + alpha ) > 0 | (4)
	 *  +-------------------+
	 * 
	 *  - Tp is additive in NN model (similarly with DG0 or DG37) (5)
	 */
	NNParameters Tp;

	/**
	 * Calculates free energy at a given temperature (in Celsius) from enthalpy and entropy
	 * 
	 * \param dh enthalpy ( kcal / mol )
	 * \param ds entropy ( cal / K / mol )
	 * \param t  temperature ( K )
	 */
	inline double dg( const double dh, const double ds, const double t ) {
		return dh - t * ds / 1000.0; // convert ds to kcal / mol
	}
	
	/**
	 * Pseudo-melting temperature (K) @ 1M NaCl; 1000x multiplier
	 * 
	 * Use this to fill in Tp with calculations using values from tables
	 * 
	 * \param dh enthalpy ( kcal / mol )
	 * \param ds entropy ( cal / K / mol)
	 * \param Tx experimental temperature (K)
	 * \param cT strand concentration (mM)
	 */
	inline int pseudoT( const double dh, const double ds ) const {
// 	NOTE: convert all values to mol
		const double dg0 = ( dh * 1000 ) - Tx * ds; // 1000 converts DH from kcal/mol to cal/mol
		const double _Tp = dg0 / rho;

		const int Tp = int( round( _Tp * 1000 )); // 1000x multiplier

		return Tp;
	};
	/**
	 * Salt-dependent pseudo-melting temperature (K) 1000x
	 * 
	 * \param Tp pseudo-melting temperature (K) 1000x
	 * \param n  oligonucleotide length
	 * 
	 * Use this to compare with experimental temperature
	 */
	inline int pseudoTSalt( const int Tp, const Length n ) const {
		const int _Tp = Tp - _lambda_rho * ( n-1 );

		return _Tp;
	};
	/**
	 * Actual melting temperature, from pseudo-melting temperature and entropy
	 * dynamic programming terms
	 * 
	 * \param Tp pseudo-melting temperature (K) 1000x
	 * \param ds entropy ( cal / K / mol ) 10x
	 */
	inline double tm( const int Tp, const int ds, const Length n ) const {
		const double alpha_n = ds / 10.0 / rho + _lambda_rho / 1000.0 * ( n-1 );
		const double Tp_n = Tp / 1000.0 - _lambda_rho / 1000.0 * ( n-1 );

		assert( 1 + alpha_n > 0 );

		const double _tm = Tx - ( Tx - Tp_n ) / ( 1 + alpha_n );

		return _tm;
	};

	/**
	 * Salt-corrected free energy
	 */
	inline double dgSalt( const int dg, const Length n ) const {
		return dg / 100.0 - lambda * n * Tx ;
	};
	/**
	 * Melting temperature (Kelvin) from free energy, enthalpy, by applying salt correction
	 */
	inline double dg2tm( const int dg, const int ds, const Length n ) const {
		const double dss = ds / 10.0 + lambda * n;
		const double dgs = dgSalt( dg, n ) * 1000.0;

		return dgs / dss + Tx;
	};
	inline double dgds2dh( const int dg37, const int ds ) const {
		return dg37 / 100.0 + Tx * ds / 10.0 / 1000.0;
	};
protected:
	/**
	 * Initializes nearest-neighbor parameters and dimer/dimer values from published tables
	 * 
	 * NOTE: Not all dimer/dimer values are published; they need to be inferred from their reciprocal; example:
	 *       TG/GA is inferred from AG/GT 
	 */
	virtual void init(
		const map<string,array<double,3>>& with_nn_parameters,
		const vector<pair<string,array<double,3>>>& with_nn_wc
	) {
		const double dhi   = with_nn_parameters.at( "Initiation" ).at( 0 );
		const double dsi   = with_nn_parameters.at( "Initiation" ).at( 1 );
		const double dgi37 = with_nn_parameters.at( "Initiation" ).at( 2 );

		const double dgi0  = dg( dhi, dsi, Tx );
		assert( abs( dgi37 - dg( dhi, dsi, T37C )) < 0.05 );

		DH._initiation =          int( round( dhi   * 10 ));
		DS._initiation =          int( round( dsi   * 10 ));
		DG37._initiation =        int( round( dgi37 * 100 ));
		DG._initiation =          int( round( dgi0   * 100 ));

		const double dht   = with_nn_parameters.at( "Terminal AT penalty" ).at( 0 );
		const double dst   = with_nn_parameters.at( "Terminal AT penalty" ).at( 1 );
		const double dgt37 = with_nn_parameters.at( "Terminal AT penalty" ).at( 2 );

		const double dgt0  = dg( dht, dst, Tx );
		assert( abs( dgt37 - dg( dht, dst, T37C )) < 0.01 );

		DH._terminal_at_penalty   = int( round( dht   * 10 ));
		DS._terminal_at_penalty   = int( round( dst   * 10 ));
		DG37._terminal_at_penalty = int( round( dgt37 * 100 ));
		DG._terminal_at_penalty   = int( round( dgt0  * 100 ));

		const double dhs   = with_nn_parameters.at( "Symmetry correction" ).at( 0 );
		const double dss   = with_nn_parameters.at( "Symmetry correction" ).at( 1 );
		const double dgs37 = with_nn_parameters.at( "Symmetry correction" ).at( 2 );

		const double dgs0  = dg( dhs, dss, Tx );
		assert( abs( dgs37 - dg( dhs, dss, T37C )) < 0.05 );

		DH._symmetry_correction   = int( round( dhs   * 10 ));
		DS._symmetry_correction   = int( round( dss   * 10 ));
		DG37._symmetry_correction = int( round( dgs37 * 100 ));
		DG._symmetry_correction   = int( round( dgs0  * 100 ));

// 	TODO: make initialization prettier
		Tp._initiation          = Thermo::pseudoT( dhi, dsi );
		Tp._terminal_at_penalty = Thermo::pseudoT( dht, dst );
		Tp._symmetry_correction = Thermo::pseudoT( dhs, dss );

// 	Watson/Crick pairs (perfect duplex matches)
		for( const auto& e: with_nn_wc ) { // generate reciprocal
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.02 );

			DH._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dh * 10 ));
			DS._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( ds * 10 ));
			DG37._nn.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg37 * 100 ));
			DG._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg0  * 100 ));

			Tp._nn.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = Thermo::pseudoT( dh, ds );
		}

		for( const auto& e: with_nn_wc ) { // overwrite reciprocal, where applicable
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.02 );

			DH._nn.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dh * 10 ));
			DS._nn.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( ds * 10 ));
			DG37._nn.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dg37 * 100 ));
			DG._nn.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dg0  * 100 ));

			Tp._nn.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = Thermo::pseudoT( dh, ds );
		}

// 	mismatches
		for( const auto& e: nn_mismatch ) { // first, generate the reverse
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.1 );

			DH._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dh * 10 ));
			DS._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( ds * 10 ));
			DG37._nn.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg37 * 100 ));
			DG._nn.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg0  * 100 ));

			Tp._nn.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = Thermo::pseudoT( dh, ds );
		}
		for( const auto& e: nn_mismatch ) { // overwrite reciprocal, where applicable
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.1 );

			DH._nn.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dh * 10 ));
			DS._nn.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( ds * 10 ));
			DG37._nn.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dg37 * 100 ));
			DG._nn.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 )   = int( round( dg0  * 100 ));

			Tp._nn.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = Thermo::pseudoT( dh, ds );
		}

// 	terminal mismatches
		for( const auto& e: terminal_mismatch ) { // first, generate the reverse
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.1 );

			DH._terminal_mismatch.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dh * 10 ));
			DS._terminal_mismatch.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( ds * 10 ));
			DG37._terminal_mismatch.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg37 * 100 ));
			DG._terminal_mismatch.at(   4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = int( round( dg0  * 100 ));

			Tp._nn.at( 4 * b2 + 1 * b1 + 64 * a2 + 16 * a1 ) = Thermo::pseudoT( dh, ds );
		}
		for( const auto& e: terminal_mismatch ) { // write direct
			const Symbol a1 = nu2pre.at( asc2nu.at( e.first.at( 0 )));
			const Symbol a2 = nu2pre.at( asc2nu.at( e.first.at( 1 )));

			const Symbol b1 = nu2pre.at( asc2nu.at( e.first.at( 3 )));
			const Symbol b2 = nu2pre.at( asc2nu.at( e.first.at( 4 )));

			const double dh = e.second.at( 0 );
			const double ds = e.second.at( 1 );
			const double dg37 = e.second.at( 2 );

			const double dg0 = dg( dh, ds, Tx );

			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.1 );

			DH._terminal_mismatch.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dh * 10 ));
			DS._terminal_mismatch.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( ds * 10 ));
			DG37._terminal_mismatch.at( 4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dg37 * 100 ));
			DG._terminal_mismatch.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = int( round( dg0  * 100 ));

			Tp._terminal_mismatch.at(   4 * a1 + 1 * a2 + 64 * b1 + 16 * b2 ) = Thermo::pseudoT( dh, ds );
		}

// 	dangling ends
		for( const auto& e: nn_dang ) {
			if( e.first.at( 2 ) == '/' ) { // dangling ab/c
				const Symbol a = nu2pre.at( asc2nu.at( e.first.at( 0 )));
				const Symbol b = nu2pre.at( asc2nu.at( e.first.at( 1 )));

				const Symbol c = nu2pre.at( asc2nu.at( e.first.at( 3 )));

				const double dh = e.second.at( 0 );
				const double ds = e.second.at( 1 );
				const double dg37 = e.second.at( 2 );

				const double dg0 = dg( dh, ds, Tx );
				assert( abs( dg37 - dg( dh, ds, T37C )) < 0.06 );

				DH._dangx.at(   a * 4 + b + c * 16 ) = int( round( dh * 10 ));
				DS._dangx.at(   a * 4 + b + c * 16 ) = int( round( ds * 10 ));
				DG37._dangx.at( a * 4 + b + c * 16 ) = int( round( dg37 * 100 ));
				DG._dangx.at(   a * 4 + b + c * 16 ) = int( round( dg0  * 100 ));

				Tp._dangx.at( a * 4 + b + c * 16 ) = Thermo::pseudoT( dh, ds );

				continue;
			}

			if( e.first.at( 1 ) == '/' ) { // dangling a/bc
				const Symbol a = nu2pre.at( asc2nu.at( e.first.at( 0 )));

				const Symbol b = nu2pre.at( asc2nu.at( e.first.at( 2 )));
				const Symbol c = nu2pre.at( asc2nu.at( e.first.at( 3 )));

				const double dh = e.second.at( 0 );
				const double ds = e.second.at( 1 );
				const double dg37 = e.second.at( 2 );

				const double dg0 = dg( dh, ds, Tx );
				assert( abs( dg37 - dg( dh, ds, T37C )) < 0.5 );

				DH._dangy.at(   a + b * 16 + c * 4 ) = int( round( dh * 10 ));
				DS._dangy.at(   a + b * 16 + c * 4 ) = int( round( ds * 10 ));
				DG37._dangy.at( a + b * 16 + c * 4 ) = int( round( dg37 * 100 ));
				DG._dangy.at(   a + b * 16 + c * 4 ) = int( round( dg0  * 100 ));

				Tp._dangy.at( a + b * 16 + c * 4 ) = Thermo::pseudoT( dh, ds );

				continue;
			}

			assert( false ); // make sure all entries are either "AB/C" or "A/BC"
		}

// 	loops
		for( int i = 1 ; i < 512 ; i++ ) { // initialize all loop elements
			if( nn_loop.count( i )) {
				const array<double,3>& lo = nn_loop.at( i );

				const double gl = lo.at( 0 );
// 	(SantaLucia and Hicks 2004; Table4; a)
// 	WARNING: reference table has incorrect formula (MISSING minus!)
// 	DS = - DG37 × 1000 / 310.15 = - DG37 × 1000 / ( 273.15 + 37 )
// 	DH = 0
				const double sl = - gl * 1000 /( T37C );
				const double hl = 0;

				const double gl0 = dg( hl, sl, Tx );
				assert( abs( gl - dg( hl, sl, T37C )) < 0.00001 );

				DG37._loop.at( i ) = int( round( gl * 100 ));
				DG._loop.at( i )   = int( round( gl0 * 100 ));
				DH._loop.at( i )   = int( round( hl * 10 ));
				DS._loop.at( i )   = int( round( sl * 10 ));

				Tp._loop.at( i ) = Thermo::pseudoT( hl, sl );

				if( i >= 256 ) continue; // no bulges or hairpins of larger sizes

				const double gb = lo.at( 1 );
				const double sb = - gb * 1000 /( T37C );
				const double hb = 0;

				const double gb0 = dg( hb, sb, Tx );
				assert( abs( gb - dg( hb, sb, T37C )) < 0.00001 );

				DG37._bulge.at( i ) = int( round( gb * 100 ));
				DG._bulge.at( i )   = int( round( gb0 * 100 ));
				DH._bulge.at( i )   = int( round( hb * 10  ));
				DS._bulge.at( i )   = int( round( sb * 10  ));

				Tp._bulge.at( i ) = Thermo::pseudoT( hb, sb );

				const double gh = lo.at( 2 );
				const double sh = - gh * 1000 /( T37C );
				const double hh = 0;

				const double gh0 = dg( hh, sh, Tx );
				assert( abs( gh - dg( hh, sh, T37C )) < 0.00001 );

				DG37._hairpin.at( i ) = int( round( gh * 100 ));
				DG._hairpin.at( i )   = int( round( gh0 * 100 ));
				DH._hairpin.at( i )   = int( round( hh * 10  ));
				DS._hairpin.at( i )   = int( round( sh * 10  ));

				Tp._hairpin.at( i ) = Thermo::pseudoT( hh, sh );

				continue;
			}

// Here, there is no entry for loop of length i
// Use Jacobson-Stockmayer extrapolation (SantaLucia and Hicks 2004; eq. 7)
			const auto& phi = *( nn_loop.crbegin());
			const array<double,3>& hi = phi.second; // DG37 for the last element
			const unsigned int x = phi.first;       // length of last element

			const double gl = hi.at( 0 ) + 2.44 * R / 1000 * ( T37C ) * log ( double( i ) / x );
			const double sl = - hi.at( 0 ) * 1000 / ( T37C ) - 2.44 * R * log ( double( i ) / x );
			const double hl = 0;

			const double gl0 = dg( hl, sl, Tx );
			assert( abs( gl - dg( hl, sl, T37C )) < 0.00001 );

// DG37 (loop-n) = DG37 (loop-x) + 2.44 × R × 310.15 × ln(n/x) (Jacobson-Stockmayer)
			DG37._loop.at( i ) = int( round( gl * 100 ));
			DG._loop.at( i )   = int( round( gl0 * 100 ));
			DH._loop.at( i )   = int( round( hl * 10  ));
			DS._loop.at( i )   = int( round( sl * 10  ));

			Tp._loop.at( i ) = Thermo::pseudoT( hl, sl );

			if( i >= 256 ) continue; // no bulges or hairpins of larger sizes

			const double gb = hi.at( 1 ) + 2.44 * R / 1000 * ( T37C ) * log ( double( i ) / x );
			const double sb = - hi.at( 1 ) * 1000 / ( T37C ) - 2.44 * R * log ( double( i ) / x );
			const double hb = 0;

			const double gb0 = dg( hb, sb, Tx );
			assert( abs( gb - dg( hb, sb, T37C )) < 0.00001 );

			DG37._bulge.at( i ) = int( round( gb * 100 ));
			DG._bulge.at( i )   = int( round( gb0 * 100 ));
			DH._bulge.at( i )   = int( round( hb * 10  ));
			DS._bulge.at( i )   = int( round( sb * 10  ));

			Tp._bulge.at( i ) = Thermo::pseudoT( hb, sb );

			const double gh = hi.at( 2 ) + 2.44 * R / 1000 * ( T37C ) * log ( double( i ) / x );
			const double sh = - hi.at( 2 ) * 1000 / ( T37C ) - 2.44 * R * log ( double( i ) / x );
			const double hh = 0;

			const double gh0 = dg( hh, sh, Tx );
			assert( abs( gh - dg( hh, sh, T37C )) < 0.00001 );

			DG37._hairpin.at( i ) = int( round( gh * 100 ));
			DG._hairpin.at( i )   = int( round( gh0 * 100 ));
			DH._hairpin.at( i )   = int( round( hh * 10  ));
			DS._hairpin.at( i )   = int( round( sh * 10  ));

			Tp._hairpin.at( i ) = Thermo::pseudoT( hh, sh );
		}

		for( const auto& h: hairpin_increments ) {
			const string s = convertAsc2Nu( h.first );

			const double dg37 = h.second.at( 0 );
			const double dh   = h.second.at( 1 );
			const double ds   = - ( dg37 - dh ) * 1000 / ( T37C );

			const double dg0 = dg( dh, ds, Tx );
			assert( abs( dg37 - dg( dh, ds, T37C )) < 0.00001 );

			DG37._hairpin_increments.emplace( s, int( round( dg37 * 100 )));
			DG._hairpin_increments.emplace(   s, int( round( dg0  * 100 )));
			DH._hairpin_increments.emplace(   s, int( round( dh   * 10  )));
			DS._hairpin_increments.emplace(   s, int( round( ds   * 10  )));

			Tp._hairpin_increments.emplace(   s, Thermo::pseudoT( dh, ds ));
		}

		for( const auto& ml: nn_multiloop ) {
			for( unsigned int i = 0 ; i < 256 ; i++ ) {
				if( i < ml.second.size()) { // actual values in initializer

// 	The value +2.6 kcal mol-1 is added to the values in the table to derive the total multiloop penalty.
// 	(SantaLucia and Hicks 2004; table S4; a)
					const double dg37 = ml.second.at( i ) + 2.6;
					const double dh   = 0;
					const double ds   = - dg37 * 1000 / ( T37C );

					const double dg0 = dg( dh, ds, Tx );
					assert( abs( dg37 - dg( dh, ds, T37C )) < 0.00001 );

// 	(SantaLucia and Hicks 2004; Table4; a)
// 	DS = DG37 × 1000 / 310.15 = DG37 × 1000 / ( 273.15 + 37 )
					DG37._multiloop.at( i, ml.first -3 ) = int( round( dg37 * 100 ));
					DG._multiloop.at(   i, ml.first -3 ) = int( round( dg0  * 100 ));
					DH._multiloop.at(   i, ml.first -3 ) = int( round( dh   * 10  ));
					DS._multiloop.at(   i, ml.first -3 ) = int( round( ds   * 10  ));

					Tp._multiloop.at(   i, ml.first -3 ) = Thermo::pseudoT( dh, ds );

					continue;
				}

				const int x     = ml.second.size() - 1;
				const double hi = ml.second.at( x );

// DG37 (loop-n) = DG37 (loop-x) + 2.44 × R × 310.15 × ln(n/x) (Jacobson-Stockmayer)
				const double dg37 = hi + 2.44 * R / 1000 * ( T37C ) * log ( double( i ) / x );
				const double dh   = 0;
				const double ds   = - dg37 * 1000 / ( T37C );

				const double dg0 = dg( dh, ds, Tx );
				assert( abs( dg37 - dg( dh, ds, T37C )) < 0.00001 );

				DG37._multiloop.at( i, ml.first -3 ) = int( round( dg37 * 100 ));
				DG._multiloop.at(   i, ml.first -3 ) = int( round( dg0  * 100 ));
				DH._multiloop.at(   i, ml.first -3 ) = int( round( dh   * 10  ));
				DS._multiloop.at(   i, ml.first -3 ) = int( round( ds   * 10  ));

				Tp._multiloop.at(   i, ml.first -3 ) = Thermo::pseudoT( dh, ds );
			}
		}
	};

	friend ostream& operator<< ( ostream& o, const Thermo& t ) { // TEST
		o << "DG37" << endl << t.DG37 << endl;
		o << "DG"   << endl << t.DG   << endl;
		o << "DH"   << endl << t.DH   << endl;
		o << "DS"   << endl << t.DS   << endl;

		for( const auto& h: t.hairpin_increments ) {
			o
				<< h.first
				<< '\t' << t.DG37._hairpin_increments.at( convertAsc2Nu( h.first ))
				<< '\t' << t.DG._hairpin_increments.at( convertAsc2Nu( h.first )) << '\t' << t.hairpin_increments.at( h.first ).at( 0 )
				<< '\t' << t.DH._hairpin_increments.at( convertAsc2Nu( h.first )) << '\t' << t.hairpin_increments.at( h.first ).at( 1 )
				<< '\t' << t.DS._hairpin_increments.at( convertAsc2Nu( h.first ))
				<< endl;
		}

		return o;
	};

	/**
	 * Table of published thermodynamic parameters for sequence matching
	 * 
	 * \see (SantaLucia and Hicks 2004; table 1)
	 */
	const map<string,array<double,3>> nn_parameters = {
// 	                                DH     DS   DG37
// 	                             kcal/mol  e.u. kcal/mol
		{ "Initiation",          {  +0.2,  -5.7, +1.96 }},
		{ "Terminal AT penalty", {  +2.2,  +6.9, +0.05 }},
		{ "Symmetry correction", {   0.0,  -1.4, +0.43 }},
	};

	/**
	 * Table of published Thermo perfect WC matches
	 * 
	 * \see (SantaLucia and Hicks 2004; table 1)
	 */
	const vector<pair<string,array<double,3>>> nn_wc = {
// 		               DH     DS   DG37
// 		           kcal/mol  e.u. kcal/mol
		{ "AA/TT", {  -7.6, -21.3, -1.00 }},
		{ "AT/TA", {  -7.2, -20.4, -0.88 }},
		{ "TA/AT", {  -7.2, -21.3, -0.58 }},
		{ "CA/GT", {  -8.5, -22.7, -1.45 }},
		{ "GT/CA", {  -8.4, -22.4, -1.44 }},
		{ "CT/GA", {  -7.8, -21.0, -1.28 }},
		{ "GA/CT", {  -8.2, -22.2, -1.30 }},
		{ "CG/GC", { -10.6, -27.2, -2.17 }},
		{ "GC/CG", {  -9.8, -24.4, -2.24 }},
		{ "GG/CC", {  -8.0, -19.9, -1.84 }}
	};

/**
 * Table of published single Thermo mismatches
 * 
 * \see (Allawi and SantaLucia 1997), (Allawi and SantaLucia 1998a), (Allawi and SantaLucia 1998b) and (Allawi and SantaLucia 1998c)
 * 
 * The values for DG37 correspond to \see (SantaLucia and Hicks 2004; table 2)
 */
	const vector<pair<string,array<double,3>>> nn_mismatch = {
// 		              DH     DS   DG37
// 		            kcal/mol  e.u. kcal/mol
		{ "AA/TG", {  -0.6,  -2.3,  0.14 }},  // G/A mismatches (Allawi and SantaLucia 1998a; table 4)
		{ "AG/TA", {  -0.7,  -2.3,  0.02 }},
		{ "CA/GG", {  -0.7,  -2.3,  0.03 }},
		{ "CG/GA", {  -4.0, -13.2,  0.11 }},
		{ "GA/CG", {  -0.6,  -1.0, -0.25 }},
		{ "GG/CA", {   0.5,   3.2, -0.52 }},
		{ "TA/AG", {   0.7,   0.7,  0.42 }},
		{ "TG/AA", {   3.0,   7.4,  0.74 }},

		{ "AG/TT", {   1.0,   0.9,  0.71 }},  // G/T mismatches (Allawi and SantaLucia 1997; table 5)
		{ "AT/TG", {  -2.5,  -8.3,  0.07 }},
		{ "CG/GT", {  -4.1, -11.7, -0.47 }},
		{ "CT/GG", {  -2.8,  -8.0, -0.32 }},
		{ "GG/CT", {   3.3,  10.4,  0.08 }},
// 		{ "GG/TT", {   5.8,  16.3,  0.74 }},  // ignore terms that are not specifed in (SantaLucia and Hicks 2004; table 2) -- double mismatches
		{ "GT/CG", {  -4.4, -12.3, -0.59 }},
// 		{ "GT/TG", {   4.1,   9.5,  1.15 }},
		{ "TG/AT", {  -0.1,  -1.7,  0.43 }},
// 		{ "TG/GT", {  -1.4,  -6.2,  0.52 }},
		{ "TT/AG", {  -1.3,  -5.3,  0.34 }},

		{ "AC/TT", {   0.7,   0.2,  0.64 }},  // C/T mismatches (Allawi and SantaLucia 1998b; table 4)
		{ "AT/TC", {  -1.2,  -6.2,  0.73 }},
		{ "CC/GT", {  -0.8,  -4.5,  0.62 }},
		{ "CT/GC", {  -1.5,  -6.1,  0.40 }},
		{ "GC/CT", {   2.3,   5.4,  0.62 }},
		{ "GT/CC", {   5.2,  13.5,  0.98 }},
		{ "TC/AT", {   1.2,   0.7,  0.97 }},
		{ "TT/AC", {   1.0,   0.7,  0.75 }},

		{ "AA/TC", {   2.3,   4.6,  0.88 }},  // A/C mismatches (Allawi and SantaLucia 1998c; table 4)
		{ "AC/TA", {   5.3,  14.6,  0.77 }},
		{ "CA/GC", {   1.9,   3.7,  0.75 }},
		{ "CC/GA", {   0.6,  -0.6,  0.79 }},
		{ "GA/CC", {   5.2,  14.2,  0.81 }},
		{ "GC/CA", {  -0.7,  -3.8,  0.47 }},
		{ "TA/AC", {   3.4,   8.0,  0.92 }},
		{ "TC/AA", {   7.6,  20.2,  1.33 }},

// 		{ "AA/TC", {  -0.8,  -3.8,  0.39 }},  // ignore terms for pH 5.0
// 		{ "AC/TA", {  -6.3, -20.2, -0.02 }},
// 		{ "CA/GC", {  -4.2, -13.6,  0.02 }},
// 		{ "CC/GA", {  -1.3,  -4.9,  0.23 }},
// 		{ "GA/CC", {  -3.3, -10.3, -0.10 }},
// 		{ "GC/CA", {  -4.9, -14.7, -0.33 }},
// 		{ "TA/AC", {  -2.1,  -7.6,  0.26 }},
// 		{ "TC/AA", {   2.2,   4.8,  0.70 }},

		{ "AA/TA", {   1.2,   1.7,  0.61 }},  // AA CC GG TT mismatches (Peyret etal 1999)
		{ "CA/GA", {  -0.9,  -4.2,  0.43 }},
		{ "GA/CA", {  -2.9,  -9.8,  0.17 }},
		{ "TA/AA", {   4.7,  12.9,  0.69 }},
		{ "AC/TC", {   0.0,  -4.4,  1.33 }},
		{ "CC/GC", {  -1.5,  -7.2,  0.70 }},
		{ "GC/CC", {   3.6,   8.9,  0.79 }},
		{ "TC/AC", {   6.1,  16.4,  1.05 }},
		{ "AG/TG", {  -3.1,  -9.5, -0.13 }},
		{ "CG/GG", {  -4.9, -15.3, -0.11 }},
		{ "GG/CG", {  -6.0, -15.8, -1.11 }},
		{ "TG/AG", {   1.6,   3.6,  0.44 }},
		{ "AT/TT", {  -2.7, -10.8,  0.69 }},
		{ "CT/GT", {  -5.0, -15.8, -0.12 }},
		{ "GT/CT", {  -2.2,  -8.4,  0.45 }},
		{ "TT/AT", {   0.2,  -1.5,  0.68 }},
	};

	/**
	 * Table of terminal mismatches
	 * 
	 * NOTE: referenced as "(Varma and SantaLucia) manuscript in preparation", but not available in literature
	 * NOTE: values here are obtained from DINAMelt and DINAFold
	 */
	const vector<pair<string,array<double,3>>> terminal_mismatch = {
// 		              DH       DS   DG37
// 		            kcal/mol  e.u. kcal/mol
		{ "AC/AG",  {  -8.4, -21.1, -1.9 }},
		{ "AC/CG",  {  -6.5, -17.2, -1.2 }},
		{ "AC/GG",  { -10.2, -28.4, -1.4 }},
		{ "CC/AG",  {  -6.5, -16.4, -1.5 }},
		{ "CC/CG",  {  -4.6, -12.5, -0.8 }},
		{ "CC/TG",  {  -8.8, -25.6, -0.9 }},
		{ "GC/AG",  {  -7.2, -18.0, -1.6 }},
		{ "GC/GG",  {  -9.0, -25.3, -1.2 }},
		{ "GC/TG",  {  -9.5, -27.2, -1.1 }},
		{ "TC/CG",  {  -4.2, -11.0, -0.8 }},
		{ "TC/GG",  {  -7.9, -22.2, -1.0 }},
		{ "TC/TG",  {  -8.4, -24.1, -0.9 }},

		{ "AG/AC",  {  -9.6, -26.5, -1.4 }},
		{ "AG/CC",  {  -6.3, -17.5, -0.9 }},
		{ "AG/GC",  {  -6.9, -20.4, -0.6 }},
		{ "CG/AC",  {  -9.9, -28.2, -1.2 }},
		{ "CG/CC",  {  -6.6, -19.2, -0.7 }},
		{ "CG/TC",  {  -9.2, -26.9, -0.9 }},
		{ "GG/AC",  {  -9.8, -27.2, -1.4 }},
		{ "GG/GC",  {  -7.1, -21.1, -0.6 }},
		{ "GG/TC",  {  -9.1, -25.9, -1.1 }},
		{ "TG/CC",  {  -7.5, -21.3, -0.9 }},
		{ "TG/GC",  {  -8.1, -24.2, -0.6 }},
		{ "TG/TC",  { -10.1, -29.0, -1.1 }},

		{ "AA/AT", {   -0.5,   1.6, -1.0 }},
		{ "AA/CT", {    4.6,  17.1, -0.7 }},
		{ "AA/GT", {   -1.4,  -1.2, -1.0 }},
		{ "CA/AT", {   -0.1,   2.6, -0.9 }},
		{ "CA/CT", {    5.0,  18.1, -0.6 }},
		{ "CA/TT", {    3.5,  13.6, -0.7 }},
		{ "GA/AT", {   -1.8,  -2.2, -1.1 }},
		{ "GA/GT", {   -2.7,  -5.1, -1.1 }},
		{ "GA/TT", {    1.8,   8.8, -0.9 }},
		{ "TA/CT", {   -2.5,  -5.1, -0.9 }},
		{ "TA/GT", {   -8.5, -23.5, -1.2 }},
		{ "TA/TT", {   -4.0,  -9.6, -1.0 }},

		{ "TA/AA",  {   1.1,   5.5, -0.6 }},
		{ "TA/AC",  {  -0.1,   0.0, -0.1 }},
		{ "TA/AG",  {   4.0,  13.2, -0.1 }},
		{ "TC/AA",  {   1.6,   6.7, -0.5 }},
		{ "TC/AC",  {   0.0,   0.0,  0.0 }}, // does not fold in DINAMelt
		{ "TC/AT",  {   4.3,  14.1, -0.1 }},
		{ "TG/AA",  {  -2.5,  -6.5, -0.5 }},
		{ "TG/AG",  {   0.0,   0.0,  0.0 }}, // does not fold in DINAMelt
		{ "TG/AT",  {   0.2,   0.9, -0.1 }},
		{ "TT/AC",  {   0.0,   0.0,  0.0 }}, // does not fold in DINAMelt
		{ "TT/AG",  {   0.0,   0.0,  0.0 }}, // does not fold in DINAMelt
		{ "TT/AT",  {   4.3,  14.1, -0.1 }},
	};
	/**
	 * Dangling end contributions
	 * 
	 * \see (Bommarito et al 2000)
	 */
	const map<string,array<double,3>> nn_dang = {
// 		              DH       DS   DG37
// 		            kcal/mol  e.u. kcal/mol
		{ "AA/T", {    0.2,    2.3, -0.51 }},
		{ "AC/G", {   -6.3,  -17.1, -0.96 }},
		{ "AG/C", {   -3.7,  -10.0, -0.58 }},
		{ "AT/A", {   -2.9,   -7.6, -0.50 }},

		{ "CA/T", {    0.6,    3.3, -0.42 }},
		{ "CC/G", {   -4.4,  -12.6, -0.52 }},
		{ "CG/C", {   -4.0,  -11.9, -0.34 }},
		{ "CT/A", {   -4.1,  -13.0, -0.02 }},

		{ "GA/T", {   -1.1,   -1.6, -0.62 }},
		{ "GC/G", {   -5.1,  -14.0, -0.72 }},
		{ "GG/C", {   -3.9,  -10.9, -0.56 }},
		{ "GT/A", {   -4.2,  -15.0,  0.48 }},

		{ "TA/T", {   -6.9,  -20.0, -0.71 }},
		{ "TC/G", {   -4.0,  -10.9, -0.58 }},
		{ "TG/C", {   -4.9,  -13.8, -0.61 }},
		{ "TT/A", {   -0.2,   -0.5, -0.10 }},

		{ "A/AT", {   -0.7,   -0.8, -0.48 }},
		{ "C/AG", {   -2.1,   -3.9, -0.92 }},
		{ "G/AC", {   -5.9,  -16.5, -0.82 }},
		{ "T/AA", {   -0.5,   -1.1, -0.12 }},

		{ "A/CT", {    4.4,   14.9, -0.19 }},
		{ "C/CG", {   -0.2,   -0.1, -0.23 }},
		{ "G/CC", {   -2.6,   -7.4, -0.31 }},
		{ "T/CA", {    4.7,   14.2,  0.28 }}, 

		{ "A/GT", {   -1.6,   -3.6, -0.50 }},
		{ "C/GG", {   -3.9,  -11.2, -0.44 }},
		{ "G/GC", {   -3.2,  -10.4, -0.01 }},
		{ "T/GA", {   -4.1,  -13.1, -0.01 }},

		{ "A/TT", {    2.9,   10.4, -0.29 }},
		{ "C/TG", {   -4.4,  -13.1, -0.35 }},
		{ "G/TC", {   -5.2,  -15.0, -0.52 }},
		{ "T/TA", {   -3.8,  -12.6,  0.13 }},
	};
	const map<unsigned int, array<double,3>> nn_loop = {
// 	 	Loop  Internal  Bulge   Hairpin
// 	 	size    loop    loop     loop
		{ 1,  {  0.0,    4.0,    0.0 }},
		{ 2,  {  0.0,    2.9,    0.0 }},
		{ 3,  {  3.2,    3.1,    3.5 }},
		{ 4,  {  3.6,    3.2,    3.5 }},
		{ 5,  {  4.0,    3.3,    3.3 }},
		{ 6,  {  4.4,    3.5,    4.0 }},
		{ 7,  {  4.6,    3.7,    4.2 }},
		{ 8,  {  4.8,    3.9,    4.3 }},
		{ 9,  {  4.9,    4.1,    4.5 }},
		{ 10, {  4.9,    4.3,    4.6 }},
		{ 12, {  5.2,    4.5,    5.0 }},
		{ 14, {  5.4,    4.8,    5.1 }},
		{ 16, {  5.6,    5.0,    5.3 }},
		{ 18, {  5.8,    5.2,    5.5 }},
		{ 20, {  5.9,    5.3,    5.7 }},
		{ 25, {  6.3,    5.6,    6.1 }},
		{ 30, {  6.6,    5.9,    6.3 }},
	};

	const map<string,array<double,2>> hairpin_increments = {
		{ "AGAAT", { -1.5, -1.5 }},
		{ "AGCAT", { -1.5, -1.5 }},
		{ "AGGAT", { -1.5, -1.5 }},
		{ "AGTAT", { -1.5, -1.5 }},
		{ "CGAAG", { -2.0, -2.0 }},
		{ "CGCAG", { -2.0, -2.0 }},
		{ "CGGAG", { -2.0, -2.0 }},
		{ "CGTAG", { -2.0, -2.0 }},
		{ "GGAAC", { -2.0, -2.0 }},
		{ "GGCAC", { -2.0, -2.0 }},
		{ "GGGAC", { -2.0, -2.0 }},
		{ "GGTAC", { -2.0, -2.0 }},
		{ "TGAAA", { -1.5, -1.5 }},
		{ "TGCAA", { -1.5, -1.5 }},
		{ "TGGAA", { -1.5, -1.5 }},
		{ "TGTAA", { -1.5, -1.5 }},

		{ "AAAAAT", {  0.7,  0.5 }},
		{ "AAAACT", {  0.2,  0.7 }},
		{ "AAACAT", {  0.5,  1   }},
		{ "ACTTGT", { -1.3,  0   }},
		{ "AGAAAT", { -1.6, -1.1 }},
		{ "AGAGAT", { -1.6, -1.1 }},
		{ "AGATAT", { -2,   -1.5 }},
		{ "AGCAAT", { -2.1, -1.6 }},
		{ "AGCGAT", { -1.6, -1.1 }},
		{ "AGCTTT", { -0.3,  0.2 }},
		{ "AGGAAT", { -1.6, -1.1 }},
		{ "AGGGAT", { -1.6, -1.1 }},
		{ "AGGGGT", {  0.3,  0.5 }},
		{ "AGTAAT", { -2.1, -1.6 }},
		{ "AGTGAT", { -1.6, -1.1 }},
		{ "AGTTCT", {  0.3,  0.8 }},
		{ "ATTCGT", { -0.7, -0.2 }},
		{ "ATTTGT", { -0.5,  0 }},
		{ "ATTTTT", { -1,   -0.5 }},
		{ "CAAAAG", {  0.9,  0.5 }},
		{ "CAAACG", {  0.7,  0.7 }},
		{ "CAACAG", {  1,    1   }},
		{ "CAACCG", {  0,    0   }},
		{ "CCTTGG", { -0.8,  0   }},
		{ "CGAAAG", { -1.1, -1.1 }},
		{ "CGAGAG", { -1.1, -1.1 }},
		{ "CGATAG", { -1.5, -1.5 }},
		{ "CGCAAG", { -1.6, -1.6 }},
		{ "CGCGAG", { -1.1, -1.1 }},
		{ "CGCTTG", {  0.2,  0.2 }},
		{ "CGGAAG", { -1.1, -1.1 }},
		{ "CGGGAG", { -1,   -1   }},
		{ "CGGGGG", {  0.8,  0.5 }},
		{ "CGTAAG", { -1.6, -1.6 }},
		{ "CGTGAG", { -1.1, -1.1 }},
		{ "CGTTCG", {  0.8,  0.8 }},
		{ "CTTCGG", { -0.2, -0.2 }},
		{ "CTTTGG", {  0,    0   }},
		{ "CTTTTG", { -0.5, -0.5 }},
		{ "GAAAAC", {  1.5,  0.5 }},
		{ "GAAACC", {  0.7,  0.7 }},
		{ "GAACAC", {  1,    1   }},
		{ "GCTTGC", { -0.8,  0   }},
		{ "GGAAAC", { -1.1, -1.1 }},
		{ "GGAGAC", { -1.1, -1.1 }},
		{ "GGATAC", { -1.6, -1.6 }},
		{ "GGCAAC", { -1.6, -1.6 }},
		{ "GGCGAC", { -1.1, -1.1 }},
		{ "GGCTTC", {  0.2,  0.2 }},
		{ "GGGAAC", { -1.1, -1.1 }},
		{ "GGGGAC", { -1.1, -1.1 }},
		{ "GGGGGC", {  0.8,  0.5 }},
		{ "GGTAAC", { -1.6, -1.6 }},
		{ "GGTGAC", { -1.1, -1.1 }},
		{ "GGTTCC", {  0.8,  0.8 }},
		{ "GTTCGC", { -0.2, -0.2 }},
		{ "GTTTGC", {  0,    0   }},
		{ "GTTTTC", { -0.5, -0.5 }},
		{ "GAAAAT", {  1.5,  0.5 }},
		{ "GAAACT", {  1,    1   }},
		{ "GAACAT", {  1,    1   }},
		{ "GCTTGT", { -0.5,  0   }},
		{ "GGAAAT", { -1.1, -1.1 }},
		{ "GGAGAT", { -1.1, -1.1 }},
		{ "GGATAT", { -1.6, -1.6 }},
		{ "GGCAAT", { -1.6, -1.6 }},
		{ "GGCGAT", { -1.1, -1.1 }},
		{ "GGCTTT", { -0.1, -0.1 }},
		{ "GGGAAT", { -1.1, -1.1 }},
		{ "GGGGAT", { -1.1, -1.1 }},
		{ "GGGGGT", {  0.8,  0.5 }},
		{ "GGTAAT", { -1.6, -1.6 }},
		{ "GGTGAT", { -1.1, -1.1 }},
		{ "GTATAT", { -0.5, -0.5 }},
		{ "GTTCGT", { -0.4, -0.4 }},
		{ "GTTTGT", { -0.4, -0.4 }},
		{ "GTTTTT", { -0.5, -0.5 }},
		{ "TAAAAA", {  0.4,  0.5 }},
		{ "TAAACA", {  0.2,  0.7 }},
		{ "TAACAA", {  0.5,  1   }},
		{ "TCTTGA", { -1.3,  0   }},
		{ "TGAAAA", { -1.6, -1.1 }},
		{ "TGAGAA", { -1.6, -1.1 }},
		{ "TGATAA", { -2.1, -1.6 }},
		{ "TGCAAA", { -2.1, -1.6 }},
		{ "TGCGAA", { -1.6, -1.1 }},
		{ "TGCTTA", { -0.3,  0.2 }},
		{ "TGGAAA", { -1.6, -1.1 }},
		{ "TGGGAA", { -1.6, -1.1 }},
		{ "TGGGGA", {  0.3,  0.5 }},
		{ "TGTAAA", { -2.1, -1.6 }},
		{ "TGTGAA", { -1.6, -1.1 }},
		{ "TGTTCA", {  0.3,  0.8 }},
		{ "TTTCGA", { -0.7, -0.2 }},
		{ "TTTTGA", { -0.5,  0   }},
		{ "TTTTTA", { -1,   -0.5 }},
		{ "TAAAAG", {  1,    0.5 }},
		{ "TAAACG", {  0.5,  1   }},
		{ "TAACAG", {  0.5,  1   }},
		{ "TCTTGG", { -1,    0   }},
		{ "TGAAAG", { -1.5, -1   }},
		{ "TGAGAG", { -1.5, -1   }},
		{ "TGATAG", { -2,   -1.5 }},
		{ "TGCAAG", { -2,   -1.5 }},
		{ "TGCGAG", { -1.5, -1   }},
		{ "TGCTTG", { -0.6, -0.1 }},
		{ "TGGAAG", { -1.5, -1   }},
		{ "TGGGAG", { -1.5, -1   }},
		{ "TGGGGG", {  0.3,  0.5 }},
		{ "TGTAAG", { -2,   -1.5 }},
		{ "TGTGAG", { -1.5, -1   }},
		{ "TTTCGG", { -0.9, -0.4 }},
		{ "TTTTAG", { -1.5, -1   }},
		{ "TTTTGG", { -0.9, -0.4 }},
		{ "TTTTTG", { -1,   -0.5 }},
	};

	/**
	 * The value +2.6 kcal mol -1 is added to the values in the table to derive the total
	 * multiloop penalty. For lengths greater than those shown in the table, a Jacobson-
	 * Stockmayer extrapolation is used.
	 */

	const map<int,array<double,7>> nn_multiloop = {
// 	                                # of SS residues
// 		                 0     1     2     3     4     5     6
// 		# of helices
		{ 3,         {  2.0,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6 }}, 
		{ 4,         { -1.0, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6 }},
		{ 5,         {  2.0,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0 }},
		{ 6,         {  2.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2 }},
		{ 7,         {  2.0,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4 }},
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
