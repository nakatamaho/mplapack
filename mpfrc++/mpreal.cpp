/*
    Multi-precision real number class. C++ wrapper fo MPFR library.
    Project homepage: http://www.holoborodko.com/pavel/
    Contact e-mail:   pavel@holoborodko.com

    Copyright (c) 2008-2010 Pavel Holoborodko

        Forked and adapted by Nakata Maho
        https://github.com/nakatamaho/mplapack
        Contact e-mail: maho.nakata@gmail.com

    Copyright (c) 2010-2021 Nakata Maho

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

    Contributors:
    Brian Gladman, Helmut Jarausch, Fokko Beekhof, Ulrich Mutze,
    Heinz van Saanen, Pere Constans, Dmitriy Gubanov
*/

#include "mpreal.h"
#include <cstring>

using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::string;
using std::ws;

namespace mpfr {

mp_rnd_t mpreal::default_rnd = mpfr_get_default_rounding_mode();
mp_prec_t mpreal::default_prec = mpfr_get_default_prec();
int mpreal::default_base = 10;
int mpreal::double_bits = -1;

} // namespace mpfr
