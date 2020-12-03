#include "Error.h"

stringstream Error :: output;
mutex Error :: lock;

Error::Error( const char* file, int line, const char* message ) {
	lock.lock();
	output.str( "" );
	output << "*** error : " << message;
	output << endl << "*** source: " << file << " (" << line << ")";
	lock.unlock();
};

const char* Error::what() const noexcept {
// 	WARNING: allocates storage for the error string; this will never be released
	return strdup( output.str().c_str());
};

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
