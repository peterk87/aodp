#define __ThermoStructure_cpp__

#include "ThermoStructure.h"

ThermoStructure::ThermoStructure(
	const string& _s, const Position _lo, const Position _le, deque<Length>& _o, const int _mi, const int _ma, const Thermo& _th, ostream& _out ) :
	th( _th ), s( _s ), lo( _lo ), le( _le ), o( _o ), mi( _mi ), ma( _ma ), out( _out )
{
	assert( mi <= ma );
	assert( s.size() == o.size());
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
