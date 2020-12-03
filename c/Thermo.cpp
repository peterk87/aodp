#define __Thermo_cpp__

#include "Thermo.h"

using namespace std;

const double Thermo::R = 1.9872; // gas constant in cal/K/mol
const double Thermo::K = 273.15; // 0C in Kelvin
const double Thermo::T37C = Thermo::K + 37.0; // 37C in Kelvin
const double Thermo::x = 1.0;    // strand concentration divisor (SantaLucia and Hicks 2004; eq.3); always 1 for oligonucleotides

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
