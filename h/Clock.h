#ifndef __Clock_h__
#define __Clock_h__

#include <iostream>

#include <unistd.h>
#include <sys/times.h>

#include <sys/resource.h>

using namespace std;

/**
 * Print "maximum resident set size"
 * 
 * Reference: man getrusage
 * WARNING: may not display correct information on all architectures
 */
inline ostream& operator<< ( ostream& o, struct rusage& r );

class Clock
{
public:
	/**
	 * Constructor sets the start time
	 */
	Clock( ostream* o = NULL ) {
		setOutput( o );
		start();
	};

	inline void setOutput( ostream* o ){
		out = o;
	}

	/**
	 * Sets the start time (e.g. if not at initialization)
	 */
	inline void start(){
		e0 = double( times( &t )) / ticks;
		e1 = e0;

		u0 = double( t.tms_utime ) / ticks;
		u1 = u0;

		s0 = double( t.tms_stime ) / ticks;
		s1 = s0;
	};

	/**
	 * Displays elapsed time since last check or start
	 */
	inline void check( string message ) {
		if( !out )
			return;

		*out << *this << message << endl;
	};
	/**
	 * Displays elapsed time since since start
	 */
	inline void stop( string message = "TOTAL" ) {
		if( !out )
			return;

		e1=e0; u1 = u0; s1=s0;
		*out << "-----------------------------------------------" << endl;
		check( message );
	};

	friend ostream& operator << ( ostream& o, Clock& c );
private:
	struct rusage r;
	struct tms t;

	/**
	 * Elapsed time in seconds
	 * 
	 * Start time of clock; start time of current operation
	 */
	double e0, e1;

	/**
	 * "User" time in seconds
	 * 
	 * Start time of clock; start time of current operation
	 */
	double u0, u1;

	/**
	 * "System" time in seconds
	 * 
	 * Start time of clock; start time of current operation
	 */
	double s0, s1;

	ostream* out;
	static const long ticks;
};

ostream& operator << ( ostream& o, Clock& c );

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
