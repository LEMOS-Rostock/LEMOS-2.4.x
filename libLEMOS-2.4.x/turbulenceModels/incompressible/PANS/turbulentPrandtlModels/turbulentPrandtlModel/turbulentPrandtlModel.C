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

#include "turbulentPrandtlModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace PANSModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPrandtlModel, 0);
defineRunTimeSelectionTable(turbulentPrandtlModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPrandtlModel::turbulentPrandtlModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& param,
    const word& turbulenceModelName
)
:
    IOdictionary
    (
        IOobject
        (
            "PANSProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    mesh_(U.mesh()),
    param_(param),
    coeffDict_(subOrEmptyDict(turbulenceModelName + "Coeffs")),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    )

{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    //mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<turbulentPrandtlModel> turbulentPrandtlModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& param,
    const word& turbulenceModelName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database multiple times
    
    IOdictionary PANSPropertiesDict
    (
        IOobject
        (
            "PANSProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    dictionary turbulenceModelCoeffsDict(PANSPropertiesDict.subOrEmptyDict(turbulenceModelName + "Coeffs"));


    const word modelType
    (
        turbulenceModelCoeffsDict.subOrEmptyDict("turbulentPrandtlModels").lookupOrAddDefault(param, word("zeroTransportModel"))
    );

    Info<< "Selecting turbulent Prandtl model " << modelType << " for parameter " << param << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "turbulentPrandtlModel::New"
            "("
                "const volVectorField&, "
                "const surfaceScalarField&, "
                "transportModel&, "
                "const word&"
            ")"
        )   << "Unknown turbulentPrandtlModel type "
            << modelType << nl << nl
            << "Valid turbulentPrandtlModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<turbulentPrandtlModel>
    (
        cstrIter()(U, phi, param, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentPrandtlModel::correct()
{
}


bool turbulentPrandtlModel::read()
{
  
    bool ok = IOdictionary::readData
    (
        IOdictionary::readStream
        (
            coeffDict_.dictName()
        )
    );
    IOdictionary::close();

    if (ok)
    {
        if (const dictionary* dictPtr = subDictPtr(coeffDict_.dictName()))
        {
            coeffDict_ <<= *dictPtr;
        }

        return true;
    }
    else
    {
        return false;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
