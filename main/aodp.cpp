#define __aodp_cpp__

#include "Application.h"

int main( int argc, char* argv[] )
{
	try{
		Application a( argc, argv );
		return a.run();
	}catch(Error e){
		cerr << e.what() << endl;
		return 1;
	}
}

// This file is part of aodp (the Automated Oligonucleotide Design Pipeline)
// 
// (C)	HER MAJESTY THE QUEEN IN RIGHT OF CANADA (2014)
// (C)	Manuel Zahariev mz@alumni.sfu.ca (2000-2008, 2014)
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
