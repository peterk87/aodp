#ifndef __Range_h__
#define __Range_h__

#include <iostream>
#include <set>
#include <vector>
#include <stack>
#include <functional>
#include <limits>

using namespace std;

template<typename T> class Cover;

/**
 * Range of values of an integer and unsigned type
 * 
 * Supports the following operations:
 * !         : the range is empty
 * & T       : range contains element (set containment)
 * & Range   : range overlaps another range
 * <= Range  : range is included in another range (set inclusion)
 * == Range  : two ranges are equal
 * 
 * - Range   : complement of the range as a Cover<T>
 * += T      : amplification by a number (increase the size left and right)
 * += Range  : smallest container range
 * &= Range  : set intersection
 * >> T      : right shift by a number
 * << T      : left shift by a number
 * 
 * Range < Range : compare operator used in set<Range>
 * 
 * NOTE: prevent initialization Range< type > unless type is integer and unsigned
 * through hack in constructor
 * 
 */
template<typename T>
class Range
{
	friend class Cover<T>;

	static_assert(
		numeric_limits<T>::is_integer && !numeric_limits<T>::is_signed,
		"For Range<T>: T must be an unsigned integer type"
	);

private:
	T l;
	T h;

public:
	static const Range<T> zero;
	static const Range<T> universe;

	/**
	 * Constructor, also overloads default
	 * 
	 * \param lo low bound of range (type T)
	 * \param le length of range (type T)
	 * 
	 * NOTE: since le is unsigned, an invalid range (l > h) cannot be created
	 */
	Range(
		const T& lo,
		const T& le
	
	): l( lo ), h( lo + le ){};

// 	OPERATORS
	/**
	 * Test whether a Range is empty
	 */
	inline bool operator! () const {
		return ( l==h );
	};

	/**
	 * Test whether a Range is contains element \param e
	 * Set containment
	 */
	inline bool operator& ( T& e ) const {
		return ( l <= e ) && ( e < h );
	};

	/**
	 * Test whether two Ranges are equal
	 */
	inline bool operator== ( const Range<T>& r ) const {
		return ( l == r.l ) && ( h == r.h );
	};

	/**
	 * Set inclusion
	 */
	inline bool operator<= ( const Range<T>& r ) const {
		return ( r.l <= l ) && ( h < r.h );
	};

	/**
	 * Complement for ranges
	 */
	inline Cover<T> operator- () const {
		Cover<T> result;

		if( l > 0 )
			result += Range<T>{ 0, l };

		if( h < numeric_limits<T>::max()) {
			result += Range<T>{ h, T( numeric_limits<T>::max()-h )};
		}

		return result;
	};

	/**
	* Smallest range containing both parameter ranges
	* 
	*              (1)        |    (2)     |   (3)     |   (4)
	* ------------------------+------------+-----------+--------
	*                         |            |           |
	* r1      [-----)         | [-----)    | [------)  |  [---)
	* r2                [===) |    [=====) |  [===)    | [=====)
	* r1+=r2  [-------------) | [--------) | [------)  | [-----)
	* 
	* Of course, extrapolate for edge cases
	* 
	*/
	inline Range<T>& operator+=( const Range<T>& right ){
		this->l = min( this->l, right.l );
		this->h = max( this->h, right.h );
		return *this;
	};
	/**
	* Set intersection for Ranges
	*/
	inline Range<T>& operator&=( const Range<T>& right ){
		this->l = max( this->l, right.l );
		this->h = min( this->h, right.h );

		if( this->h < this->l ){
// 	WARNING (not thread-safe): before the next line, the range is invalid
			this->h = this->l;
		}

		return *this;
	};
	/**
	* Shift a range to the right
	*/
	inline Range<T>& operator>>( const T& right ){
		if( l > T( numeric_limits<T>::max() - right )) {
			*this = zero;
			return *this;
		}

		this->l += right;

		if( h > T( numeric_limits<T>::max() - right ))
			this->h = numeric_limits<T>::max();
		else
			this->h += right;

		return *this;
	};
	/**
	* Shift a range to the left
	*/
	inline Range<T>& operator<<( const T& right ){
		if( h < right ) {
			*this = zero;
			return *this;
		}

		this->h -= right;

		if( l < right )
			this->l = 0;
		else
			this->l -= right;

		return *this;
	};
	/**
	 * Amplify a range. Will increase the size of the range by \param a
	 * Will obey numeric_limits.
	 * 
	 * range     [----)    += 4
	 *       [3210----0123)
	 */
	inline Range<T>& operator+=( const T& a ){
		this->l = ( this->l < a ) ? 0 : ( this->l - a );
		this->h = ( this->h < ( numeric_limits<T>::max() - a )) ? ( this->h + a ) : numeric_limits<T>::max();

		return *this;
	};

	/**
	 * Whether this Range overlaps another Range
	 */
	inline bool operator& ( Range<T> r ) const {
		return
			(( this->l <= r.l ) && ( this->h > r.l ))
			|| (( r.l <= this->l ) && ( r.h > this->l ))
			;
	};

// 	METHODS
	inline T size() const {
		return h-l;
	}

	/**
	 * Low bound of the range
	 */
	inline const T& lo() const { return l; }
	/**
	 * High bound of the range
	 */
	inline const T& hi() const { return h; }

	/**
	 * Fill a range with intermediary "sites", like so:
	 * 
	 * lo               mid                  hi
	 * [     *     * ... * ...   *     *     ]
	 *  < f >|< i >| ... |       |< i >|< f >
	 *
	 * where:
	 * - lo: lower bound
	 * - hi: higher bound
	 * - mid: midpoint of range
	 * - f: \param first_site_gap
	 * - i: \param inter_site_gap
	 *
	 * \return vector of zero-based sites
	 * 
	 * NOTE:
	 *  - the low and high range will not be included in the result
	 *  - adjacent sites that are too close will be ommitted
	 *  - the return result should not be copied on return (return value optimization)
	 *    reference: http://en.wikipedia.org/wiki/Return_value_optimization
	 */
	inline vector< T > fill( const T& first_site_gap, const T& inter_site_gap ) const {
		vector< T > result;
		stack< T > s;

		T mid = ( T )(( l+h )/2 );
		T half = ( T )(( h-l )/2 );

// 	build first half of result
		for( T i=first_site_gap ; i < half ; i+= inter_site_gap ){
			result.push_back( l + i );
			s.push( i );
		}

// 	add mid point
		result.push_back( mid );

// 	build second half of result
		while( !s.empty()){
			result.push_back( h - s.top());
			s.pop();
		}

// 	WARNING: result object will live until released by caller if compiler has 
// 	used return value optimization
		return result;
	};

	/**
	* Call the function f of object o for every sub-range of size between m and M
	* 
	* NOTE: protect against size of subsequence Length being larger than size of Position
	*/
	template<typename L, class C, typename... Args>
	inline void cover( C& o, void ( C::*f )( T, L, Args... ), L m, L M, Args... args ) const {

		static_assert( 
			numeric_limits<L>::is_integer
			&& !numeric_limits<L>::is_signed
			&& ( numeric_limits<L>::digits <= numeric_limits<T>::digits ),
			"For Range<T>::cover<L...>: L must be an unsigned integer type of at most the bit size of type T"
		);

		T p = l;

		if( p > ( hi() - m ))
// 	The minimum size will not fit in the range
			return;


		while( p <= ( hi() -M )){
// 	Range size is greater than maximum size
			((o).*(f))( p, M, args... );

			p++;
		}

		while( p <= ( hi() - m )){
			((o).*(f))( p, hi()-p, args... );
			p++;
		}
	};

	/**
	 * Return the length of the largest Range:
	 *  - of length between m and M inclusive,
	 *  - starting at position \param p
	 *  - that is fully contained in this Range
	 * 
	 * If there is no Range that satisfies that, return 0.
	 * 
	 * Rationale: used in covering an area of a sequence with oligo
	 * candidates of length between m and M
	 */
	template<typename L>
	inline L cover( T p, L m, L M ) const {
		static_assert(
			numeric_limits<L>::is_integer
			&& !numeric_limits<L>::is_signed
			&& ( numeric_limits<L>::digits <= numeric_limits<T>::digits ),
			"For Range<T>::cover<L...>: L must be an unsigned integer type of at most the bit size of type T"
		);

		if( p > h ) {
			return 0;
		}

		if(( h-p ) < m ) {
			return 0;
		}

		if(( h-p ) > M ) {
			return M;
		}

		return ( h-p );
	};
	/**
	* Compare operator used in set< Range< T >>
	* 
	* Lexicographic order of ranges (lo/hi)
	*/
	inline friend bool operator<( const Range<T>& left, const Range<T>& right ){
		if( left.l < right.l ) {
			return true;
		}
		if( left.l > right.l ) {
			return false;
		}

		assert( left.l == right.l );

		return ( left.h < right.h );
	};
	/**
	* Print a range (TEST)
	*/
	inline friend ostream& operator<<( ostream& left, const Range<T>& right ){
		left << "[" << right.lo() << ":" << right.hi() << ")";
		return left;
	}
};

template<typename T> const Range<T> Range<T>::zero( 0, 0 );
template<typename T> const Range<T> Range<T>::universe( 0, numeric_limits<T>::max());


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
