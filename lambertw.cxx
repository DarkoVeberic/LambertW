/*
  Command-line interface to the Lambert W function

  Copyright (C) 2011 Darko Veberic, darko.veberic@ijs.si

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*

  This source implements an utility "lambertw" so that the
  numerical implementation of the Lambert W function can be
  easily obtained from the command line or be used in shell
  scripts.

  Usage:

  lambertw [branch] x

  - "branch" value is optional (default 0 is assumed) and can
    be only 0 or -1
  - x value for W(x); note that the definition range for
    branch 0 is [-1/e, infinity] and for branch -1 it is
    [-1/e, 0) where 1/e is approx. 0.367879

  Examples:

  ./lambertw 3.14     -->   1.073395661239825
  ./lambertw -0.2     -->  -0.2591711018190738
  ./lambertw 0 1      -->   0.5671432904097838
  ./lambertw 0 0      -->   0
  ./lambertw -1 -0.2  -->  -2.542641357773526
  ./lambertw -1 -0.3  -->  -1.781337023421628

*/

#include <LambertW.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

using namespace std;


void
Usage(const int argc, const char* const name)
{
  cout << "Usage: " << name << " [branch] x\n"
          "You gave " << argc << " parameters." << endl;
  exit(1);
}


int
main(int argc, char* argv[])
{
  int branch = 0;
  double x = 0;

  switch (argc) {
  case 3:
    {
      stringstream ss;
      ss << argv[argc-2];
      if (!(ss >> branch) || !(branch == -1 || branch == 0))
        Usage(argc, argv[0]);
    }
  case 2:
    {
      stringstream ss;
      ss << argv[argc-1];
      if (!(ss >> x))
        Usage(argc, argv[0]);
    }
    break;
  default:
    Usage(argc, argv[0]);
    break;
  }

  cout << setprecision(16)
       << utl::LambertW(branch, x)
       << endl;

  return 0;
}
