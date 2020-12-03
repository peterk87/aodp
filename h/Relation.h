#ifndef __Relation_h__
#define __Relation_h__

#include <map>
#include <unordered_set>

#include <cmath>

// TEST
#include <bitset>

using namespace std;

/**
 * Binary relation between two sets of types T1 and T2
 * 
 * M1 and M2 must be mapping types, such as map, multimap, unsorted_map, unsorted_multimap
 * In its most general form, this is a many-to-many relation.
 * 
 * \see Many2Many, One2One and One2Many for derived relations
 */
template <typename T1, typename T2, class M1 = multimap<T1,T2>, class M2 = multimap<T2,T1>> class _Relation {
public:
	M1 from;
	M2 to;

	typedef M1 type_from;
	typedef M2 type_to;

	/**
	 * \return  0 : ok
	 *          1 : fail because of k1 (implemented in derived classes)
	 *          2 : fail because of k2 (implemented in derived classes)
	 */
	virtual int emplace( const T1& k1, const T2& k2 ) {
		this->from.emplace( k1, k2 );
		this->to.emplace( k2, k1 );

		return 0;
	};

// 	TEST
	friend ostream& operator<< ( ostream& o, const _Relation& r ) {
		for( const auto& e: r.from ) {
			cout << e.first << "->" << e.second << " " << endl;
		}
		o << "---------------------" << endl;
		for( const auto& e: r.to ) {
			cout << e.second << "<-" << e.first << " " << endl;
		}
		o << "=====================";
		return o;
	};
};

/**
 * Many-to-many relation
 */
template <typename T1, typename T2, typename M1, typename M2> class Relation : public _Relation<T1, T2, M1, M2> {
public:
	inline bool has( const T1& k1 ) const {
		return ( this->from.find( k1 ) != this->from.end());
	};

	inline bool has( const T2& k2 ) const {
		return ( this->to.find( k2 ) != this->to.end());
	};

	inline pair<typename M1::const_iterator, const typename M1::const_iterator> equal_range( const T1& k1 ) const {
		return this->from.equal_range( k1 );
	};

	inline pair<typename M2::const_iterator, const typename M2::const_iterator> equal_range( const T2& k2 ) const {
		return this->to.equal_range( k2 );
	};
};

/**
 * Skip implementation of "has" and "equal_range" for relations where the "from" and "to" are the same type
 * 
 * Reasoning: has and equal_range will end up defining exactly the same function twice (for the same type)
 */
template <typename T, typename M1, typename M2> class Relation<T,T,M1,M2> : public _Relation<T, T, M1, M2> {
};

/**
 * One-to-many relation: restriction of a many-to-many relation
 */
template <typename T1, typename T2, class M1 = multimap<T1,T2>> class One2Many: public Relation<T1, T2, M1, map<T2,T1>> {
public:
	int emplace( const T1& k1, const T2& k2 ) {
// 	fails if key k2 is already contained by the "to" side of the relation
		if( this->to.find( k2 ) != this->to.end()) {
			return 2;
		}

		return Relation<T1,T2,M1,map<T2,T1>>::emplace( k1, k2 );
	};

	inline const T1& at( const T2& k2 ) const {
		return this->to.at( k2 );
	};

	/**
	 * Erase all relations containing the key on the from side
	 */
	inline void erase( const T1& k1 ) {
		for( auto i2from = this->from.equal_range( k1 ) ; i2from.first != i2from.second ; ++i2from.first ) {
			this->to.erase( i2from.first->second );
		}
		this->from.erase( k1 );
	};
};

/**
 * One-to-one relation: restriction on a one-to-many relation
 */
template <typename T1, typename T2, class M1 = map<T1,T2>> class One2One: public One2Many<T1, T2, M1> {
public:
	virtual int emplace( const T1& k1, const T2& k2 ) {
		if( this->from.find( k1 ) != this->from.end()) {
// 	fails if k1 is already contained by the "from" side of the relation
			return 1;
		}

// 	defer to the One2Many emplace; may still fail is k2 is already in the to side of the relation
		return One2Many<T1,T2,M1>::emplace( k1, k2 );
	};

	using One2Many<T1, T2, M1>::at;
	inline const T2& at( const T1& k1 ) const {
		return this->from.at( k1 );
	};

	/**
	 * Erase the relation containing the key on the right side
	 */
	inline void erase( const T2& k2 ) {
		if( this->to.find( k2 ) == this->to.end()) {
			return;
		}
		this->from.erase( this->to.at( k2 ));
		this->to.erase( k2 );
	};
};

/**
 * Many-to-many relation for unsigned numeric types
 */
template <typename T1, typename T2, typename T = unsigned long long > class Many2Many : public unordered_set<T> {
	static_assert( numeric_limits<T1>::is_integer && !numeric_limits<T1>::is_signed,
				   "For Many2Many<T1,T2,T>: T1 must be an unsigned integer type" );
	static_assert( numeric_limits<T2>::is_integer && !numeric_limits<T1>::is_signed,
				   "For Many2Many<T1,T2,T>: T2 must be an unsigned integer type" );
	static_assert( numeric_limits<T>::is_integer &&  !numeric_limits<T>::is_signed,
				   "For Many2Many<T1,T2,T>: T  must be an unsigned integer type" );

	static_assert( numeric_limits<T1>::digits + numeric_limits<T2>::digits <= numeric_limits<T>::digits,
				   "For Many2Many<T1,T2,T>: size of type T must be larger than the combined sizes of T1 and T2" );

protected:
	static const size_t s1 = sizeof( T1 ) * 8;
	static const T f1 = pow( 2, s1 ) - 1;

	static T combine( const T1& k1, const T2& k2 ) {
		return ( T{ k1 } << s1 ) | k2;
	};
public:
	static pair<T1,T2> split( const T& k ) {
		return make_pair( T1( k >> s1 ), T2( k & f1 ));
	};
	inline bool emplace( const T1& k1, const T2& k2 ) {
		if( has( k1, k2 )) {
			return false;
		}
		this->unordered_set<T>::emplace( combine( k1, k2 ));
		return true;
	};

	inline bool has( const T1& k1, const T2& k2 ) const {
		return this->count( combine( k1, k2 ));
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
