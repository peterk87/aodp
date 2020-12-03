#ifndef __array2_h__
#define __array2_h__

#include <iostream>
#include <iomanip>

#include <array>
#include <limits>

using namespace std;

/**
 * Two-dimensional array; horizontal dimension wraps around (WARNING)
 * 
 * Use as a dynamic programming buffer, e.g. for thermodynamic calculations
 */
template <typename T, size_t h, size_t v, const T zero = 0>
class array2: public array<array<T,v>,h> {
public:
	array2( const T& f = zero ) {
		for( size_t li = 0 ; li < h ; li++ ) {
			this->array<array<T,v>,h>::at( li ).fill( f );
		}
	};
	inline T& at( size_t i, size_t j ) {
		return this->array<array<T,v>,h>::at( i % h ).at( j );
	};
	inline const T& at( size_t i, size_t j ) const {
		return this->array<array<T,v>,h>::at( i % h ).at( j );
	};

	inline bool zeroAt( const size_t i, const size_t j ) {
		return ( this->at( i, j ) == zero );
	};

	/**
	 * Write a matrix with single-characters column headings
	 * 
	 * TEST
	 */
	void put( ostream& o, const int n, const int i, const string& s, const int w = 3 ) {
		o << string( w+1, ' ' );
		for( int ii=0 ; ii < n ; ii++ ) {
			o << setw( w ) << ( i+ii ) % h << ' ';
		}
		o << endl;
		o << string( w+1, ' ' );
		for( int ii=0 ; ii < n ; ii++ ) {
			o << setw( w ) << s.at( i+ii ) << ' ';
		}
		o << endl;

		for( int j = 0 ; j < n ; j++ ) {
			o << setw( w ) << j << ' ';

			for( int ii=0 ; ii < ( n-j ) ; ii++ ) {
				o << setw( w );

				if( this->zeroAt(( i+ii ) % h, j )) o << '.';
				else o << this->at(( i+ii ) % h, j );

				o << ' ';
			}
			o << endl;
		}
	};

	inline void set( const size_t i, const size_t j, const T& y ) {
		T& x = this->at( i, j );

		if( y == zero ) return;
		x = y;
	};
	template<typename... Args>
	void set( const size_t i, const size_t j, const T& y1, const T& y2, Args... args ) {
		if( y1 == zero ) {
			set( i, j, y2, args... );
			return;
		}

		if( y2 == zero ) {
			set( i, j, y1, args... );
			return;
		}

		set( i, j, y1+y2, args... );
	};

	/**
	* Replace element x at position (\param i, \param j) with \param y if y > x
	* 
	* \returns true if the copy happened
	*/
	inline bool copyMore( const size_t i, const size_t j, const T& y ) {
		T& x = this->at( i, j );

		if( y == zero ) return false; // WARNING: if y==zero, the replacement will not happen

		if( x == zero ) {
			x = y;
			return true;
		}

		if( x >= y ) return false;

		x = y;
		return true;
	};
	/**
	 * Replace element x at position (\param i, \param j) with the sum of the remaining arguments,
	 * if sum > x
	 * 
	 * NOTE: recursive definition
	 */
	template<typename... Args>
	inline bool copyMore( const size_t i, const size_t j, const T& y1, const T& y2, Args... args ) {
		if( y1 == zero ) return copyMore( i, j, y2, args... );
		if( y2 == zero ) return copyMore( i, j, y1, args... );

		return copyMore( i, j, y1+y2, args... );
	};

	/**
	* Replace element x at position (\param i, \param j) with \param y if y < x
	* 
	* \returns true if the copy happened
	*/
	inline bool copyLess( const size_t i, const size_t j, const T& y ) {
		T& x = this->at( i, j );

		if( y == zero ) return false; // WARNING: if y==zero, the replacement will not happen

		if( x == zero ) {
			x = y;
			return true;
		}

		if( x <= y ) return false;

		x = y;
		return true;
	};
	/**
	 * Replace element x at position (\param i, \param j) with the sum of the remaining arguments,
	 * if sum < x
	 * 
	 * NOTE: recursive definition
	 */
	template<typename... Args>
	inline bool copyLess( const size_t i, const size_t j, const T& y1, const T& y2, Args... args ) {
		if( y1 == zero ) return copyLess( i, j, y2, args... );
		if( y2 == zero ) return copyLess( i, j, y1, args... );

		return copyLess( i, j, y1+y2, args... );
	};
};


/**
 * "Square" array
 */
template <typename T, size_t l, const T zero = 0>
using sq = array2<T, l, l, zero>;

/**
 * "Square" array with different zero
 */
template <typename T, size_t l>
using sqm = array2<T, l, l, numeric_limits<T>::max()>;

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
