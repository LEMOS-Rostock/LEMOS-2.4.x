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

Class
    hybDynV2F

Description
    Hybrid model: Durbin V2F turbulence model in wall region + 
    dynamic model of Smagorinsky in outer region   

SourceFiles
    hybDynV2F.C

\*---------------------------------------------------------------------------*/

#ifndef hybDynV2F_H
#define hybDynV2F_H

#include "incompressible/RAS/RASModel/RASModel.H"

#include "../../../LES/LESdeltas/LESdelta/LESdelta.H"
#include "../../../LES/LESfilters/LESfilter/LESfilter.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class DurbinV2F Declaration
\*---------------------------------------------------------------------------*/

class hybDynV2F
:
    public RASModel 

{
    // Private data

        dimensionedScalar Cmu;
        dimensionedScalar CmuKE;
        dimensionedScalar Ceps10;
        dimensionedScalar Ceps11;
        dimensionedScalar Ceps2;
        dimensionedScalar C1;
        dimensionedScalar C2;
        dimensionedScalar CL;
        dimensionedScalar CEta;
        dimensionedScalar oneOnSigmaK;
        dimensionedScalar oneOnSigmaEps;
        dimensionedScalar yStarLim;
        dimensionedScalar Cint;
        dimensionedScalar u_tau_averaged;

        dimensionedScalar f0_;

        wallDist yw_;

        volScalarField k_;
        volScalarField epsilon_;
        volScalarField v2_;
        volScalarField f_;

        volScalarField dMax;
        volScalarField nut_;
        volScalarField rans;
        volScalarField cD_field;

	    autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;
        autoPtr<LESdelta> delta_;

 
    // Private member functions
        dimensionedScalar cD(const volSymmTensorField& D) const;

        tmp<volScalarField> T_() const;
        tmp<volScalarField> L_() const;

public:

    //- Runtime type information
    TypeName("hybDynV2F");

    // Constructors

        //- from components
        hybDynV2F
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    // Destructor

        ~hybDynV2F(){}


    // Member Functions

        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DkEff", oneOnSigmaK*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField> 
            (
                new volScalarField("DepsilonEff", oneOnSigmaEps*nut_ + nu())
            );
        }

        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Note that epsilonTilda is returned as epsilon.
        //  This is the appropriate variable for most purposes.
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        tmp<volScalarField> v2() const
        {
            return v2_;
        }

        tmp<volScalarField> f() const
        {
            return f_;
        }


        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }


        tmp<volSymmTensorField> R() const;

//         tmp<fvVectorMatrix> divR(volVectorField& U) const;
    
        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;
        
        virtual tmp<fvVectorMatrix> divDevRhoReff
        ( 
            const volScalarField& rho, 
            volVectorField& U 
        ) const ;
        

        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
