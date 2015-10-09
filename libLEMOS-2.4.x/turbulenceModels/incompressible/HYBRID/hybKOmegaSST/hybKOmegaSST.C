/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "hybKOmegaSST.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hybKOmegaSST, 0);
addToRunTimeSelectionTable(RASModel, hybKOmegaSST, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

dimensionedScalar hybKOmegaSST::cD(const volSymmTensorField& D) const
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

tmp<volScalarField> hybKOmegaSST::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> hybKOmegaSST::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybKOmegaSST::hybKOmegaSST
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            1.0
        )
    ),

    CN_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CN",
            coeffDict_,
            1.0
        )
    ),
       
    d_
    (
        IOobject
        (
            "d",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("d", dimless, 0.0)
    ),   
 
    x1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "x1",
            coeffDict_,
            0.95
        )
    ),  
     
    x2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "x2",
            coeffDict_,
            1.05
        )
    ), 
 
    omegaSmall_
    (
        "omegaSmall", 
         dimensionSet(0, 0, -1, 0, 0), 
         SMALL
    ),
    
    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),

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
        dimensionedScalar("rans", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),


    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_()),
    delta_(LESdelta::New("delta", U.mesh(), *this))
    
{
    bound(omega_, omegaMin_);

    nut_ = a1_*k_/max(a1_*omega_, F2()*mag(symm(fvc::grad(U_))));
    nut_.correctBoundaryConditions();

    printCoeffs();
    
    const pointField& ps = mesh_.points();
    const edgeList& es = mesh_.edges();
    const labelListList& eCells = mesh_.cellEdges();

    scalar edgeLength = 0.0 ;
    volScalarField maxEdgeLength = y_;
    
    forAll( maxEdgeLength , I )
    {
        edgeLength = 0.0;
        forAll( eCells[I], J )
        {
	        if ( ( es[eCells[I][J]].mag( ps ) ) > edgeLength ) 
	            edgeLength = es[ eCells[I][J] ].mag( ps );
        }
        maxEdgeLength[I] = edgeLength;
    }
      
    forAll( d_ , I )
    {
        d_[I] = sqrt( 0.5*( maxEdgeLength[I]*maxEdgeLength[I] + delta()[I]*delta()[I] ) ); 
    }
    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> hybKOmegaSST::R() const
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> hybKOmegaSST::devReff() const
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


tmp<fvVectorMatrix> hybKOmegaSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<fvVectorMatrix> hybKOmegaSST::divDevRhoReff
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

bool hybKOmegaSST::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CN_.readIfPresent(coeffDict());
        x1_.readIfPresent(coeffDict());
        x2_.readIfPresent(coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


void hybKOmegaSST::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField S2 = magSqr(symm(fvc::grad(U_)));
    volScalarField G(GName(), nut_*2*S2);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega =
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_;

    volScalarField F1 = this->F1(CDkOmega);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*2*S2
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    nut_.correctBoundaryConditions();

    volScalarField lPrandtl = CN_ *((100.0/9.0)*(pow(k_,0.5)/(omega_ + omegaSmall_))); // integral length scale using Prandtl formula
    volSymmTensorField D = dev(symm(fvc::grad(U_)));
    volScalarField nuSGS = mag(cD(D))*sqr(delta())*sqrt(2*magSqr(D)); // SGS stress using DSM
 
    scalar LtoDRatio = 0;
    
    scalar x1 = x1_.value(); 
    scalar x2 = x2_.value();
    scalar xk = 0; // value between x1 and x2 - coordinate for a smoothing function
    
 
    forAll( d_, I )
    {
	    rans[I] = 1; 
	    
	    LtoDRatio = lPrandtl[I] / d_[I] ;
	    LtoDRatio /= Cint_.value();    // how many cells to resolve a vortex?

	    xk = ( LtoDRatio - x1 ) / ( x2 - x1 ); 
	    if ( LtoDRatio > x2 ) 
	    {
	        nut_[I] = nuSGS[I];   // cell is in LES regio
	        rans[I] = 0;
	    }
	    else
	    {
	        if ( LtoDRatio > x1 ) // cell is buffer region between URANS and LES (see hybrid model desc.)
	        {
	            nut_[I] =  ( nut_[I] - nuSGS[I] ) / 3.1415*atan(-40*xk/(x2-x1) + 10*(x2+x1)/(x2-x1) )  +  0.5*(nut_[I]+nuSGS[I]);
	            rans[I] = 1 - xk;
	        }
	    }        
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
