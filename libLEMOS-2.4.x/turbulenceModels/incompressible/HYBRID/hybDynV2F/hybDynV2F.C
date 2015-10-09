/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "hybDynV2F.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hybDynV2F, 0);
addToRunTimeSelectionTable(RASModel, hybDynV2F, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> hybDynV2F::T_() const
{
    volScalarField yStar_=pow(CmuKE,0.25)*sqrt(k_)*yw_/nu();
    return max(k_/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL)),
               pos(yStarLim-yStar_)*6.0*sqrt(nu()/(epsilon_  +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL))));
}

tmp<volScalarField> hybDynV2F::L_() const
{
    volScalarField yStar_=pow(CmuKE,0.25)*sqrt(k_)*yw_/nu();
    return
        CL*max(pow(k_,1.5)/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL)),
               pos(yStarLim-yStar_)
               *CEta*pow(pow(nu(),3.0)/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL)),0.25));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
hybDynV2F::hybDynV2F
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const word& turbulenceModelName,
    const word& modelName
)
:   
    RASModel(modelName, U, phi, lamTransportModel, turbulenceModelName),
    //RASModel(typeName, U, phi, lamTransportModel),


    Cmu
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    
    
    CmuKE
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKE",
            coeffDict_,
            0.09
        )
    ),
    Ceps10
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps10",
            coeffDict_,
            1.40
        )
    ),
    
    Ceps11
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps11",
            coeffDict_,
            0.05
        )
    ),
    
    Ceps2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.90
        )
    ),
    
    C1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.40
        )
    ),
    C2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.30
        )
    ),
    CL
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            1.92
        )
    ),
    CEta
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            70.0
        )
    ),
    
    oneOnSigmaK
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaK",
            coeffDict_,
            1.00
        )
    ),
    
    oneOnSigmaEps
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "oneOnSigmaEps",
            coeffDict_,
            0.76923
        )
    ),
    
    yStarLim
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "yStarLim",
            coeffDict_,
            30.0
        )
    ),

    Cint
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            1.0
        )
    ),
    u_tau_averaged
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "u_tau_averaged",
            coeffDict_,
            0.2
        )
    ),
    
    

    f0_("f0small", dimless/dimTime, SMALL),

    yw_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    dMax(yw_),

    // Calculate viscosity (with Davidson correction - 2003)
    nut_(min(CmuKE*sqr(k_)/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL)), Cmu*v2_*T_())),

    rans
    (
        IOobject
        (
            "rans",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("rans", dimless, 0.0)
    ),
    cD_field
    (
	IOobject
        (
           "cD",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("cD", dimless, 0.0)
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_()),
    delta_(LESdelta::New("delta", U.mesh(), *this))


{
    printCoeffs();

    const pointField& ps = mesh_.points();
    const edgeList& es = mesh_.edges();
    const labelListList& eCells = mesh_.cellEdges();

    scalar xxx =0.0;
    
    forAll(dMax, I)
    {
      xxx = 0.0;
      forAll(eCells[I], J)
      {
	if ((es[eCells[I][J]].mag(ps)) > xxx) xxx = es[eCells[I][J]].mag(ps);
      }
      dMax[I] = xxx;
    }



}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> hybDynV2F::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*2*symm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

tmp<volSymmTensorField> hybDynV2F::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

tmp<fvVectorMatrix> hybDynV2F::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<fvVectorMatrix> hybDynV2F::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


dimensionedScalar hybDynV2F::cD(const volSymmTensorField& D) const
{

    volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    volSymmTensorField MM =
      sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}



bool hybDynV2F::read()
{
    if (RASModel::read())
    {
        Cmu.readIfPresent(coeffDict_);
        CmuKE.readIfPresent(coeffDict_);
        Ceps10.readIfPresent(coeffDict_);
        Ceps11.readIfPresent(coeffDict_);
        Ceps2.readIfPresent(coeffDict_);
        C1.readIfPresent(coeffDict_);
        C2.readIfPresent(coeffDict_);
        CL.readIfPresent(coeffDict_);
        CEta.readIfPresent(coeffDict_);
        oneOnSigmaK.readIfPresent(coeffDict_);
        oneOnSigmaEps.readIfPresent(coeffDict_);
        yStarLim.readIfPresent(coeffDict_);
        Cint.readIfPresent(coeffDict_);
        filter_.read(coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


void hybDynV2F::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    RASModel::correct();

    if (mesh_.moving())
    {
        yw_.correct();
    }

    volScalarField S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G = nut_*S2;

    volScalarField T_ = this->T_(); // DOESNT WORK AS T_() ?!
    volScalarField Ceps1 = Ceps10*(scalar(1.0)+Ceps11*min(sqrt(k_/v2_),scalar(10.0)));



    // Dissipation rate equation

//#   include "epsilonV2FWallI.H"
/////////////////////////////////
{
    labelList cellBoundaryFaceCount(epsilon_.size(), 0);

    const fvPatchList& patches = mesh_.boundary();

    //- Initialise the near-wall epsilon fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] = 0.0;
            }
        }
    }
    
    volScalarField nu = this->nu(); // NOT SURE ABOUT THIS!!
 
    //- Accumulate the wall face contributions to epsilon
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
#           include "checkPatchFieldTypesEpsilonV2F.H"

            const scalarField& nuw = nu.boundaryField()[patchi];

            scalarField magFaceGradU = mag(U_.boundaryField()[patchi].snGrad());

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                // For corner cells (with two boundary or more faces),
                // epsilon in the near-wall cell is calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

    		    epsilon_[faceCelli] +=
    		       2.0*nuw[facei]*k_[faceCelli]/(sqr(RASModel::y_[patchi][facei]));
            }
        }
    }


    // Perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}
///////////////////////////////////////////////////
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1*G/T_
      - fvm::Sp(Ceps2/T_, epsilon_)
    );

    epsEqn().relax();
#   include "wallDissipationV2FI.H"
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);





    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp(1.0/T_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);




    // f equation
    volScalarField L_ = this->L_();

    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/sqr(L_),f_)
      - ((C1-scalar(6.0))*v2_/k_ - 2.0/3.0*(C1-scalar(1.0)))/(sqr(L_)*T_)
      + C2*G/(k_*sqr(L_))
    );

    fEqn().relax();
    solve(fEqn);
    bound(f_, f0_);




    // v2 equation

    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(DkEff(), v2_)
     ==
    //    k_*f_ 
    // Davidson correction - 2003
        min(k_*f_, 
          -((C1-scalar(6.0))*v2_ - 2.0/3.0*k_*(C1-scalar(1.0)))/T_ + C2*G)
      - fvm::Sp(6.0*epsilon_/k_, v2_)
    );

    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, kMin_);

    // Re-calculate viscosity (with Davidson correction - 2003)

 
    volScalarField lPrandtl = 0.168*pow(k_,3/2)/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL));

    nut_ = min(CmuKE*sqr(k_)/(epsilon_ +  dimensionedScalar("SMALL",dimensionSet(0,2,-3,0,0,0,0),SMALL)), Cmu*v2_*this->T_());

    //    Pout<<"We are in correct(2) "<<endl;

   
    volSymmTensorField D = dev(symm(fvc::grad(U_)));
    cD_field = cD(D);
    volScalarField nuSGS = mag(cD(D))*sqr(delta())*sqrt(2*magSqr(D));

//    k_sgs =  cD(D)*sqr(delta())*2*magSqr(D)/0.3;

    scalar kkk = 0;
    scalar xk = 0;
    scalar x1 = 0.95;
    scalar x2 = 1.05;
    
    // Hybridization of both methods
    forAll(dMax, I)
      {
        rans[I] = 1; // By default we are in RANS zone
        kkk = lPrandtl[I] / dMax[I] / Cint.value();
          xk = (kkk - x1) / (x2-x1);
          if (kkk > x2) // We are in LES zone
            {
              nut_[I] = nuSGS[I];
              rans[I] = 0;
            }
          else
            {
             if (kkk > x1) // Smooth interpolation between RANS and LES
               {
                 nut_[I] = (nut_[I]-nuSGS[I])/3.1415*atan(20*xk/(x2-x1)-10*(x2+x1)/(x2-x1))+0.5*(nut_[I]+nuSGS[I]);
                 rans[I] = 1-xk;
               }
            }
      }
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
