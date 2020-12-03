#ifndef __Error_h__
#define __Error_h__

#include <iostream>	
#include <sstream>
#include <mutex>

#include <cstring>

using namespace std;

#define error( ... ) throw Error( __FILE__, __LINE__, __VA_ARGS__ )

/**
 * Main exception class; will dump a variable arguments message to what()
 * and include the file/line of the throw (look at macro definition above)
 */
class Error: public exception {
public:
	Error( const char* file, int line, const char* message );

	template <typename... Args>
	Error( const char* file, int line, Args... args ) {
		lock.lock();
		output.str( "" );

		output << "*** error : ";
		print( args... );
		output << endl << "*** source: " << file << " (" << line << ")";
		lock.unlock();
	};

	virtual const char* what() const noexcept;
private:
	template <typename T>
	void print( T t ){
		output << t;
	};

	template <typename T, typename... Args>
	void print( T t, Args... args ){
		output << t << " ";
		print( args... );
	};

// 	int x;
	static stringstream output;
	static mutex lock;
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
