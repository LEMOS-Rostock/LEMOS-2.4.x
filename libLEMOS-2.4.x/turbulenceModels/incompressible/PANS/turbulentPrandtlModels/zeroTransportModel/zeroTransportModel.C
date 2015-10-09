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

#include "zeroTransportModel.H"
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

defineTypeNameAndDebug(zeroTransportModel, 0);
addToRunTimeSelectionTable(turbulentPrandtlModel, zeroTransportModel, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zeroTransportModel::zeroTransportModel
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


void zeroTransportModel::correct()
{
}


bool zeroTransportModel::read()
{
    turbulentPrandtlModel::read();
    sigmaEps_.readIfPresent(coeffDict_);
    sigmaK_.readIfPresent(coeffDict_);
    
    return true;
}

tmp<volScalarField> zeroTransportModel::sigmaK() const 
{

   
const volScalarField& fk = mesh_.lookupObject<volScalarField>("fk");
const volScalarField& feps = mesh_.lookupObject<volScalarField>("feps");

return  sigmaK_*sqr(fk)/feps;
};

tmp<volScalarField> zeroTransportModel::sigmaEps() const 
{

const volScalarField& fk = mesh_.lookupObject<volScalarField>("fk");
const volScalarField& feps = mesh_.lookupObject<volScalarField>("feps");

return  sigmaEps_*sqr(fk)/feps;
};

tmp<volScalarField> zeroTransportModel::sigmaOmega() const 
{};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentPrandtlModels
} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
