#ifndef __Cover_h__
#define __Cover_h__

#include <cassert>

#include "Range.h"

using namespace std;

template<typename T> class Cover;

/**
 * Set of Range<T>
 * 
 * Supports the following operations:
 * 
 * & T      : whether a Range in the set contains the element
 * 
 * !        : Cover is empty (contains no ranges)
 * += Range : add a range to Cover
 * -        : set complement
 * &= Range : change the universe of the Cover
 * 
 * >> T     : shift all ranges in the Cover to the right
 * += T     : amplify all ranges in the Cover
 * 
 * clear
 * flip
 */
template<typename T> class Cover: public set<Range<T>> {
public:
	Cover<T>(): universe( Range<T>::universe ) {
	};

	Cover<T>( const Cover<T>& c ): set<Range<T>>( c ), universe( c.universe ) {
	};
	/**
	 * Constructor that uses parameters of the "universe" range to be created
	 */
	Cover<T>( T lo, T le ) : universe( Range<T>{ lo, le }) {};

	/**
	 * Constructor from Range components
	 * used in TEST drivers
	 */
	Cover<T>( const Range<T>& ra, const set<Range<T>>& se ) : set<Range<T>>( se ), universe( ra ) {};

	/**
	 * Return the sum of the lengths of all Ranges in the Cover
	 * 
	 * WARNING: linear complexity
	 */
	T length() const {
		T l = 0;
		for( const Range<T>& r: *this ) {
			l += r.size();
		}
		return l;
	}

	/**
	 * Return the maximum size of all ranges contained in a window of size w
	 */
	template<typename L>
	inline L window( L w ) const {
		L re = 0;

		if( !*this ) {
			return re;
		}
		for( auto i = this->cbegin() ; i != this->cend() ; ++i ) {
			if( i->size() >= w ) {
				return w;
			}

			L _re = i->size();

			auto j = i;
			for(  ++j ; j != this->end() ; ++j ) { // measure all portions of ranges that fit in a window of size w after the right edge of *i
				Range<T> h{ i->hi(), w };

				h &= *j;
				if( h.size()) {

					_re += h.size();
					if( _re >= w ) {
						return w;
					}

					continue;
				}

				break;
			}

			if( _re > re ) {
				re = _re;
			}
		}

		return re;
	};

	/**
	 * Change the universe of the Cover
	 * 
	 * WARNING: some of the Ranges in the Cover may end up outside the universe
	 */
	inline Cover<T>& operator&= ( const Range<T>& u ) {
		universe = u;
		return *this;
	};

	
	inline void clear() {
		universe = Range<T>::universe;
		set<Range<T>>::clear();
	};

	/**
	 * Whether a Range in the set contains the element
	 * 
	 * WARNING: linear complexity
	 */
	inline bool operator& ( T e ) const {
		for( const Range<T>& r : *this ) {
			if( r & e )
				return true;
		}

		return false;
	};

	/**
	 * Whether the cover is empty
	 */
	inline bool operator! () const {
		return this->empty();
	};

	/**
	* Shift all ranges in a cover and the universe of the cover to the right
	*/
	inline Cover<T>& operator>> ( const T& right ){
		const Cover<T> c = *this;
		this->clear();

		this->universe &= c.universe;
		this->universe >> right;

		for( Range<T> r: c ) {
			r >> right;
			r &= universe;
			*this += r;
		}

		return *this;
	}
	/**
	* Amplify all ranges in a cover by \param a
	* 
	* Support for covering original Cover by ranges of length a+1 (maximum):
	* - chop off at left edge of universe
	* - do not chop off at right edge of universe (support for lengths down to minimum)
	* - combine amplified ranges with each other (if the un-amplified range is overlapping an amplified range)
	* 
	* WARNING: some of the Ranges in the Cover may end up outside the left edge of the universe
	* WARNING: some of the Ranges may overlap
	* WARNING: linear complexity
	* 
	* universe    --[                                 )----
	* *this       ------------[     )----------------------
	* *this += 4  --------[3210     0123]------------------
	* 
	* universe    --[                                 )----
	* *this       -----[ )---------------------------------
	* *this += 4  --[210 0123)-----------------------------    chop left
	* 
	* universe    --[                                 )----
	* *this       --------------------------------[ )------
	* *this += 4  ----------------------------[3210 0123)--    do not chop right
	* 
	*                             r1     r2
	* *this       ---------------[   )-[    )--------------
	* *this += 4  -----------[3210   0123   0123)----------    combine
	*                                  ==                        r1+=4 overlaps with r2
	* 
	*                             r1        r2
	* *this       ---------------[   )------[ )-------------
	* *this += 4  -----------[3210   0123)-----------------
	*             ----------------------[3210 0123)---------   do not combine
	*                                ====                        match r1+=4
	*                                 ~~~~                       skip (between r1+=4 and r2+=4)
	*                                  ~~~~                      skip (between r1+=4 and r2+=4)
	*                                   ====                     match r2+=4
	* 
	* universe    --[                                    )-
	* *this       ------------------------------------[ )--
	* *this += 4  --------------------------------[3210 01)    obey numeric limits!
	* 
	* \see cover
	*/
	inline Cover<T>& operator+= ( const T& a ){
		if( !*this ) {
			return *this;
		}

		Cover<T> c = *this;
// 	remove the ranges from the cover; keep the universe
		this->set<Range<T>>::clear();

		typename Cover<T>::const_iterator i = c.begin();
		while( true ) {
			Range<T> r = *i;

			r += a;
			if( r.l < universe.lo()) {
				r.l = universe.lo();
			}

			while( true ) {
				if( ++i == c.end()) {
					*this += r;
					return *this;
				}

				if( r & *i ) { // combine
					Range<T> q = *i;
					q += a;
					if( q.l < universe.lo()) {
						q.l = universe.lo();
					}

					r += q;
					continue;
				}

				*this += r;
				break;
			}
		}

		return *this;
	}

	/**
	* Add a Range to a Cover
	*/
	inline Cover<T>& operator+=( const Range<T>& right ){
		this->emplace( right );
		return *this;
	};
	/**
	* Complement of a cover
	* 
	* Restrictions: non-overlapping ranges; all ranges within the universe
	* WARNING: linear complexity
	*/
	inline Cover<T> operator- () const {
		Cover<T> result;
		result &= this->universe;

		T l = this->range().lo();

		for( const Range<T>& r: *this ) { // WARNING: iterates in the order of ranges
			assert( l <= r.lo());

			Range<T> _r( l, r.lo()-l );

			if( !!_r ) {
				result += _r;
			}

			l = r.hi();
		}

		assert( l <= universe.hi());

		Range<T> _r( l, this->universe.hi()-l );
		if( !!_r ) {
			result += _r;
		}

		return result;
	};

	/**
	* Add a Range to a set of non-overlapping Ranges
	* 
	* If the collection contains a Range overlapping the right-side parameter, combine it
	* repeatedly until the result in non-overlapping with any of the Ranges in the set
	* 
	* TODO: Cover<T>::combineWithRange
	*/
	inline Cover<T>& combineWithRange( const Range<T>& right ){
		Range<T> e = right;
		e &= universe;

		if( !e ) {
			return *this;
		}

		Cover<T> c = *this;
		this->set<Range<T>>::clear(); // do not reset the universe

// 	repeatedly combine ranges until they are non-overlapping
		for( Range<T> r: c ) {
			if( r & e ) {
				e += r;
			}
		}

		for( Range<T> r: c ) {
			if(!( r & e )) {
				*this += r;
			}
		}
		*this += e;

		return *this;
	};

	/**
	 * Reverse ("flip") the Cover so that it matches the reverse string
	 * Used for maintaining Covers for the reverse complement
	 */
	Cover<T>& flip() {
		const Cover<T> c = *this;
		this->clear();
		this->universe = c.universe;

		for( const Range<T>& r: c ) {
			this->emplace( c.universe.hi() - r.hi() + c.universe.lo(), r.size());
		}
		return *this;
	};

	/**
	* Call the function f of object o for every sub-range of size between m and M of c
	*/
	template<typename L, class C, typename... Args >
	void cover( C& o, void ( C::*f )( T lo, L le, Args... args ), L m, L M, Args... args ) const {
		for( const Range<T>& r : *this ) {
			r.cover( o, f, m, M, args... );
		}
	}

	/**
	 * Return the maximum length of a subsequence fully contained in a Range in the Cover
	 * 
	 * \param p : the position where to search for the subsequence
	 * \param m : minimum length of the subsequence being sought
	 * \param M : maximum length of the subsequence being sought
	 * \param r : iterator pointing to the Range in the Cover
	 * 
	 * \return The length of a subsequence starting at position \param p fully contained
	 *         in the range pointed to by r and the universe:
	 * 
	 * WARNING: If r points to a Range above the location of p, the return will be 0
	 *          and r will remain unchanged
	 * WARNING: not thread safe!
	 */
	template<typename L>
	inline L cover( T p, L m, L M, typename Cover<T>::iterator& r ) const {
		assert( r != this->end());

		T rh = r->hi();
		T rl = r->lo();
		T ul = universe.lo();
		T uh = universe.hi();

		if(( p < ul ) || ( p < rl )) {
			return 0;
		}

		if( p <= ( rh-M )) {
			if( p <= ( uh-M )) {
				return M;
			}

			if( p <= ( uh-m )) {
				return ( uh-p );
			}
		}

		return 0;
	}

	/**
	* Print a set of ranges (TEST)
	* 
	* WARNING: this is not an operator of the class, but an overloaded std::operator<<
	*/
	friend inline ostream& operator<<( ostream& left, const Cover<T>& right ) {
		left << right.universe << "{";
		for( Range<T> r: right ) {
			left << r;
		}
		left << "}";

		return left;
	};

	inline const Range<T>& range() const { return universe; };
private:
	/**
	 * Universe of values that applies to the Cover
	 */
	Range<T> universe;
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
