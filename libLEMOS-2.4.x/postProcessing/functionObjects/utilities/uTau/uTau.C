/*---------------------------------------------------------------------------*\ 
| File modified by LEMOS (University of Rostock) 2013                         |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "uTau.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "wallPolyPatch.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
defineTypeNameAndDebug(uTau, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::uTau::writeFileHeader(const label i)
{
    // Add headers to output data
    file() <<   "#Column 1: time" << nl 
           <<   "#Column 2: uTau based on sqrt(mag(avg(nf&Reff)))" << nl
	   <<   "#Column 3: uTau based on sqrt(mag(avg(nf&Reff)-(nf&(nf&Reff))*nf))" << nl
	   <<   "#Column 4: uTau based on sqrt(mag(avg(snGrad(U)*nu)))" << nl
	   << endl;
}


void Foam::uTau::calcFrictionVelocity
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff
)
{
    vector uTauAvg(0.0, 0.0, 0.0);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
        const symmTensorField& Reffp = Reff.boundaryField()[patchI];

        vectorField ssp = (-Sfp/magSfp) & Reffp;

        vector uTauAvgp = gAverage(ssp);
        uTauAvg += uTauAvgp;

        if (log_ > 2)
        {
            Info<< "  tauWall  avg(" << pp.name() << "): " 
                << uTauAvgp << endl;
        }
    }


    scalar numWallPatches = patchSet_.size();

    scalar magUTauAvg = sqrt(mag(uTauAvg)/numWallPatches);
        
    if (Pstream::master())
    {
            file() << mesh.time().timeName() 
                   << token::TAB
                   << magUTauAvg;
    }


     if (log_ > 0)
     {
	Info<< "uTau: " << mag(uTauAvg) << endl;
     }

}

void Foam::uTau::calcFrictionVelocity2
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff
)
{
    vector uTauAvg(0.0, 0.0, 0.0);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
        const vectorField nfp = -Sfp/magSfp;
        const symmTensorField& Reffp = Reff.boundaryField()[patchI];

        vectorField tsp = (nfp & Reffp); 
        vectorField ssp = tsp - (nfp & tsp) * nfp; 

        vector uTauAvgp = gAverage(ssp);
        uTauAvg += uTauAvgp;

        if (log_ > 2)
        {
            Info<< "  tauWall  avg(" << pp.name() << "): " 
                << uTauAvgp << endl;
        }
    }


    scalar numWallPatches = patchSet_.size();

    scalar magUTauAvg = sqrt(mag(uTauAvg)/numWallPatches);
        
    if (Pstream::master())
    {
            file() << token::TAB 
                   << magUTauAvg;
    }


     if (log_ > 0)
     {
	Info<< "uTau: " << uTauAvg << endl;
     }

}


void Foam::uTau::calcFrictionVelocity
(
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& nu
)
{
    
    vector uTauAvg(0.0, 0.0, 0.0);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];


        const vectorField wallGradU = -U.boundaryField()[patchI].snGrad() * nu.internalField();
        
        vector uTauAvgp = gAverage(wallGradU);
        uTauAvg += uTauAvgp;

        if (log_ > 2)
        {
            Info<< " tauWall avg(" << pp.name() << "): "
                <<  uTauAvgp << endl;
        }
    }

    scalar numWallPatches = patchSet_.size();

    scalar magUTauAvg = sqrt(mag(uTauAvg)/numWallPatches);

    if (Pstream::master())
    {
            file() << token::TAB
                   << magUTauAvg
                   << endl;
    }

   if (log_ > 0)
   {
   	Info<< "uTau: " << uTauAvg << endl;
   }

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uTau::uTau
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(0),
    patchSet_()
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "uTau::uTau"
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uTau::~uTau()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uTau::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<label>("log", 0);


        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        patchSet_ =
            mesh.boundaryMesh().patchSet
            (
                wordReList(dict.lookupOrDefault("patches", wordReList()))
            );

        Info<< type() << " output:" << nl;

        if (patchSet_.empty())
        {
            forAll(pbm, patchI)
            {
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    patchSet_.insert(patchI);
                }
            }

            Info<< "    processing all wall patches" << nl << endl;
        }
        else
        {
            Info<< "    processing wall patches: " << nl;
            labelHashSet filteredPatchSet;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                label patchI = iter.key();
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    filteredPatchSet.insert(patchI);
                    Info<< "        " << pbm[patchI].name() << endl;
                }
                else
                {
                    WarningIn("void uTau::read(const dictionary&)")
                        << "Requested friction velocity on non-wall boundary "
                        << "type patch: " << pbm[patchI].name() << endl;
                }
            }

            Info<< endl;

            patchSet_ = filteredPatchSet;
        }

    }
}


void Foam::uTau::execute()
{
    // Do nothing - only valid on write
}


void Foam::uTau::end()
{
    // Do nothing - only valid on write
}

void Foam::uTau::timeSet()
{
    // Do nothing
}

void Foam::uTau::write()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (active_)
    {
        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (log_ > 0)
        {
            Info<< type() << " output:" << nl;
        }


        tmp<volVectorField> U = mesh.lookupObject<volVectorField>("U");


        if (mesh.foundObject<cmpModel>("turbulenceModel"))
        {
            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceModel");
        
	    calcFrictionVelocity(mesh, model.devRhoReff());
            calcFrictionVelocity(mesh, U(), model.mu()/model.rho());
        }
        else if (mesh.foundObject<icoModel>("turbulenceModel"))
        {
            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceModel");

	    calcFrictionVelocity(mesh, model.devReff());
	    calcFrictionVelocity2(mesh, model.devReff());
            calcFrictionVelocity(mesh, U(), model.nu());


        }
        else
        {
            FatalErrorIn("void Foam::uTau::write()")
                << "Unable to find turbulence model in the "
                << "database" << exit(FatalError);
        }

        if (log_ > 0)
        {
            Info<< endl;
        }

    }
}


// ************************************************************************* //
