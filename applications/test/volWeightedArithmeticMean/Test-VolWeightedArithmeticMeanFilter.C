/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 LEMOS, University Rostock
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-VolWeightedArithmeticMeanFilter

Description
    Test app for volume weighted arithmetic mean filter.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "OFstream.H"
#include "meshTools.H"
#include "volWeightedArithmeticMeanFilter.H"
#include "simpleFilter.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    volVectorField U
    (
        IOobject
	(
	    "U",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );

    // volWeightedArithmeticMeanFilter using all point connected
    // cells or boundary faces to given cell
    volWeightedArithmeticMeanFilter filter(mesh);

    runTime++;

    // Applying the filter and write field
    filter(U)().write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
