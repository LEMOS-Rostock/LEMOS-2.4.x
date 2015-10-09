/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "plasmaActuationBodyForce.H"
#include "bodyForceModel.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(plasmaActuationBodyForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        plasmaActuationBodyForce,
        dictionary
    );
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::plasmaActuationBodyForce::plasmaActuationBodyForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    bodyForce_(bodyForceModel::New(*this, coeffs_, mesh, cells_))
{
    fieldNames_.setSize(1);
    fieldNames_[0]="U";
    applied_.setSize(fieldNames_.size(), false);
}


void Foam::fv::plasmaActuationBodyForce::writeData(Ostream&) const
{
}

bool Foam::fv::plasmaActuationBodyForce::read(const Foam::dictionary& dict)
{
  if (option::read(dict))
  {
    return true;
  }
  else
  {
    return false;
  }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::plasmaActuationBodyForce::correct(volVectorField& U)
{
    // Nothing to do here
}



void Foam::fv::plasmaActuationBodyForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    UIndirectList<vector>(Su, cells_) = bodyForce_->computeSup(eqn);

    eqn += Su;
}


void Foam::fv::plasmaActuationBodyForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    UIndirectList<vector>(Su, cells_) = rho*bodyForce_->computeSup(eqn);

    eqn += Su;
}


// ************************************************************************* //
