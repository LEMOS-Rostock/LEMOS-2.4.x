/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University of Rostock
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

#include "energyDissipationRate.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "dictionaryEntry.H"
#include "List.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(energyDissipationRate, 0);

    const word energyDissipationRate::turbulenceModelName = "turbulenceModel";
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::energyDissipationRate::compressible()
{
    if (obr_.foundObject<compressible::turbulenceModel>(turbulenceModelName))
    {
        return true;
    }
    else if (obr_.foundObject<incompressible::turbulenceModel>(turbulenceModelName))
    {
        return false;
    }
    else
    {
        WarningIn("Foam::word& Foam::energyDissipationRate::compressible() const")
            << "Turbulence model not found in database, deactivating";
        active_ = false;
    }

    return false;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyDissipationRate::energyDissipationRate
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    fieldName_(dict.lookup("fieldName")),
    obr_(obr),
    active_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "energyDissipationRate::energyDissipationRate"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);


    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        
        volTensorField* energyDissipationRateTensorPtr
        (
            new volTensorField
            (
                IOobject
                (
                    fieldName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedTensor("0", sqr(dimLength)/pow(dimTime,3), tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
            )
        );

        mesh.objectRegistry::store(energyDissipationRateTensorPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::energyDissipationRate::~energyDissipationRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyDissipationRate::read(const dictionary& dict)
{
    // Do nothing
}


void Foam::energyDissipationRate::execute()
{
    bool comp = compressible();

    if (!active_)
    {
        return;
    }

    const fvMesh& mesh = refCast<const fvMesh>(obr_);
        
    const volVectorField& UMean = mesh.lookupObject<volVectorField>("UMean");

    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    const volVectorField u = U - UMean;

    const volTensorField gradU(fvc::grad(U));
    const volTensorField gradu(fvc::grad(u));
        
    volTensorField& energyDissipationRateTensor =
            const_cast<volTensorField&>
            (
                mesh.lookupObject<volTensorField>(fieldName_)
            );
     
    if (comp)
    {
        const compressible::turbulenceModel& model =
            obr_.lookupObject<compressible::turbulenceModel>(turbulenceModelName);
        
        // not implemented yet
    }
    else
    {
        const incompressible::turbulenceModel& model =
            obr_.lookupObject<incompressible::turbulenceModel>(turbulenceModelName);

        energyDissipationRateTensor = 2*model.nuEff()*(gradu.T()&gradu);
    }
}

void Foam::energyDissipationRate::end()
{
    if (active_)
    {
        execute();
    }
}

void Foam::energyDissipationRate::timeSet()
{
    // Do nothing
}

void Foam::energyDissipationRate::write()
{
    if (active_)
    {
        regIOobject& obj =
            const_cast<regIOobject&>
            (
                obr_.lookupObject<regIOobject>(fieldName_)
            );

        obj.write();
    }
}

// ************************************************************************* //
