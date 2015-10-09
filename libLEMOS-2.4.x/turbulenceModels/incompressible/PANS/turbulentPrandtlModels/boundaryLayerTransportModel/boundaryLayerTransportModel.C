/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University Rostock
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

\*---------------------------------------------------------------------------*/

#include "boundaryLayerTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace PANSModels
{
namespace turbulentPrandtlModels
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(boundaryLayerTransportModel, 0);
addToRunTimeSelectionTable(turbulentPrandtlModel, boundaryLayerTransportModel, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

boundaryLayerTransportModel::boundaryLayerTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& param,
    const word& turbulenceModelName
)
: 
    turbulentPrandtlModel(U, phi, param, turbulenceModelName)
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    //mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //


void boundaryLayerTransportModel::correct()
{
}


bool boundaryLayerTransportModel::read()
{
    turbulentPrandtlModel::read();
    sigmaEps_.readIfPresent(coeffDict_);
    sigmaK_.readIfPresent(coeffDict_);
    
    return true;
}

tmp<volScalarField> boundaryLayerTransportModel::sigmaK() const 
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ones",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "one",
                sigmaK_.dimensions(),
                1.0
            )
        )
    );
}

tmp<volScalarField> boundaryLayerTransportModel::sigmaEps() const 
{
    const volScalarField& fk = mesh_.lookupObject<volScalarField>("fk");
    const volScalarField& feps = mesh_.lookupObject<volScalarField>("feps");
    
    return  sigmaEps_*feps/fk;
};

tmp<volScalarField> boundaryLayerTransportModel::sigmaOmega() const 
{};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentPrandtlModels
} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
