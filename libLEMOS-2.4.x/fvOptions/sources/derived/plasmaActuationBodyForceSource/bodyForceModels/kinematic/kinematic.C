/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "kinematic.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "DimensionedField.H"
#include "fvMatrices.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kinematic, 0);

    addToRunTimeSelectionTable(bodyForceModel, kinematic, dictionary);
}


Foam::scalar Foam::kinematic::calcAxialLength() const
{
  scalar mCD=mag(CD_);

  if (mCD < SMALL)
  {
      FatalErrorIn("Foam::fv::plasmaActuatorBodyForce::kinematic::calcAxialLength() const")
          <<"CD vector has zero length!"<<endl
          <<abort(FatalError);
  }

  vector el=CD_/mCD;

  scalar l = mag(el & boundBox(UIndirectList<vector>(mesh_.C(), cells_)()).span());

  return l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kinematic::kinematic
(
    const fv::plasmaActuationBodyForce& dbd,
    const dictionary& dict,
    const fvMesh& mesh,
    const labelList& cells
)
:
    bodyForceModel(dbd, dict, typeName, mesh, cells)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::kinematic::~kinematic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::kinematic::read(const dictionary& dict)
{
  if (bodyForceModel::read(dict))
  {
    coeffs_.lookup("CD") >> CD_;

    bodyForceModel::read(dict);
    
    return true;
  }
  else
  { 
    return false;
  }
}


Foam::tmp<Foam::vectorField> Foam::kinematic::computeSup(fvMatrix<vector>& eqn)
{
    tmp<vectorField> Su(new vectorField(cells_.size(), vector::zero));

    const vectorField& U = eqn.psi();

    Su() = 0.5*CD_* magSqr(UIndirectList<vector>(U, cells_)()) / calcAxialLength();

    return Su;
}


/*
void Foam::kinematic::correct
(
    const vectorField& U,
    vectorField& force
)
{
    // do nothing
}


void Foam::kinematic::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force)
{
    // do nothing
}
*/

// ************************************************************************* //
