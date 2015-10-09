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

#include "extendedForces.H"
#include "wallFvPatch.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(extendedForces, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedForces::extendedForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
: forces(name, obr, dict, loadFromFiles, readFields)
{}


Foam::extendedForces::extendedForces
(
    const word& name,
    const objectRegistry& obr,
    const labelHashSet& patchSet,
    const word& pName,
    const word& UName,
    const word& rhoName,
    const scalar rhoInf,
    const scalar pRef,
    const coordinateSystem& coordSys
)
: forces(name, obr, patchSet, pName, UName, rhoName, rhoInf, pRef, coordSys)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedForces::~extendedForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::extendedForces::execute()
{

    if (!active_)
    {
        return;
    }

    forces::execute();

    initialise();


    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    if (!pressureForce_.valid())
    {
        pressureForce_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "pressureForce",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("pressureForce", dimPressure, vector::zero),
                calculatedFvPatchField<vector>::typeName
            )
        );
    }

    if (!viscousForce_.valid())
    {
        viscousForce_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "viscousForce",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("viscousForce", dimPressure, vector::zero),
                calculatedFvPatchField<vector>::typeName
            )
        );
    }

    if (directForceDensity_)
    {
        WarningIn("extendedForces::execute") << "direct force density unsupported!" << endl;
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const fvMesh& mesh = U.mesh();


        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAll(mesh.boundaryMesh(), patchI)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchI]))
            {
                const vectorField nfb = mesh.Sf().boundaryField()[patchI] / mesh.magSf().boundaryField()[patchI];

                const symmTensorField& devRhoReffb = tdevRhoReff().boundaryField()[patchI];

                pressureForce_().boundaryField()[patchI]==
                (
                    rho(p)*nfb*(p.boundaryField()[patchI] - pRef)
                );

                viscousForce_().boundaryField()[patchI]==
                (
                    nfb & devRhoReffb
                );
            }   
        }
    }        
}

// ************************************************************************* //
