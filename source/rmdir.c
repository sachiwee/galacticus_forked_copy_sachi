// Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
//           2019
//    Andrew Benson <abenson@carnegiescience.edu>
//
// This file is part of Galacticus.
//
//    Galacticus is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Galacticus is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

//% Implements Fortran-callable wrappers around the Linux rmdir function.

#include <unistd.h>
#include <errno.h>

int rmdir_C(const char *name) {
  //% Fortran-callable wrapper around the rmdir() function to remove a directory.
  int status;
  status = rmdir(name);
  if ( status == -1 ) {
    int err = errno;
    if ( err == ENOENT ) {
      /* Directory does not exist - this is acceptable */
      return 0;
    } else {
      return err;
    }
  } else {
    return status;
  }
}
