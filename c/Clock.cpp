#define __Clock_cpp__

#include "Clock.h"

const long Clock::ticks = sysconf(_SC_CLK_TCK);

/**
 * Print "maximum resident set size"
 * 
 * Reference: man getrusage
 * WARNING: may not display correct information on all architectures
 */
inline ostream& operator<< ( ostream& o, struct rusage& r ) {
	getrusage( RUSAGE_SELF, &r );

	if( r.ru_maxrss > 1073741824L ){
		o << r.ru_maxrss / 1073741824L << " TiB";
		return o;
	}
	if( r.ru_maxrss > 1048576L ){
		o << r.ru_maxrss / 1048576L << " GiB";
		return o;
	}
	if( r.ru_maxrss > 1024L ){
		o << r.ru_maxrss / 1024L << " MiB";
		return o;
	}

	o << r.ru_maxrss << " KiB";
	return o;
}

ostream& operator << ( ostream& o, Clock& c ) {
	double e = double( times( &c.t )) / c.ticks;
	double u = double( c.t.tms_utime ) / c.ticks;
	double s = double( c.t.tms_stime ) / c.ticks;

	o.setf( ::ios::fixed );
	int pr = o.precision( 3 );

	o << ( u-c.u1 ) << '\t' << s-c.s1 << '\t' << ( e-c.e1 ) << '\t';

	o.precision( pr );
	o.unsetf( ::ios::fixed );

	o << c.r << '\t';

	c.e1 = e;
	c.u1 = u;
	c.s1 = s;

	return o;
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
