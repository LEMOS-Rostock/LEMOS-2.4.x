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

Class
    Foam::incompressible::RASModels::hybKOmegaSST

Description
    Implementation of the hybrid turbulence modelling approach combining Menters
    k-Omega-SST model and Dynamic Smagorinsky SGS model for incompressible flows.

    The Hybrid approach used is described in:
    @verbatim
        Kornev, N., Taranov, A., Shchukin, E. & Kleinsorge, L. 
        Development of hybrid URANS-LES methods for flow simulations in the ship stern area. 
        Ocean Engineering, 
        Vol. 38, 1831-1838. 
    @endverbatim

    The computational domain is dynamically divided into URANS and LES regions, depending
    on the integral length scale. If the length scale is larger than the LES filter size, the
    turbulent viscosity is calculated using LES SGS model. Otherwise, the cell is in URANS
    region and turbulent viscosity is calculated using k-Omega SST model.

    At the moment the implementation does not support arbitrary switching between different
    URANS/LES models. This means, that only the combination of k-Omega-SST and Dynamic Smagorinsky
    model is available

    @endverbatim

SourceFiles
    hybKOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef hybKOmegaSST_H
#define hybKOmegaSST_H

#include "incompressible/RAS/RASModel/RASModel.H"
#include "wallDist.H"
#include "../../../LES/LESdeltas/LESdelta/LESdelta.H"
#include "../../../LES/LESfilters/LESfilter/LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmega Declaration
\*---------------------------------------------------------------------------*/

class hybKOmegaSST
:
    public RASModel
{
    // Private data

        // Model coefficients
        dimensionedScalar alphaK1_;
        dimensionedScalar alphaK2_;

        dimensionedScalar alphaOmega1_;
        dimensionedScalar alphaOmega2_;

        dimensionedScalar gamma1_;
        dimensionedScalar gamma2_;

        dimensionedScalar beta1_;
        dimensionedScalar beta2_;

        dimensionedScalar betaStar_;

        dimensionedScalar a1_;
        dimensionedScalar c1_;
        
        // control parameters for hybrid approach
 
	    dimensionedScalar Cint_;        // how many cells we need to resolve a vortex ??
	    dimensionedScalar CN_;
	    
	    volScalarField  d_;  // characteristic size of the cell = sqrt ( 0.5 * ( dMax*dMax + deltaLES*deltaLES ) )
	    
	    dimensionedScalar x1_;          // value of L to d ratio , when the smooting starts and finises
	    dimensionedScalar x2_;          //     
	             
        dimensionedScalar omegaSmall_; 

        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
        
        wallDist y_;

        // Fields

        volScalarField k_;
        volScalarField omega_;
        volScalarField nut_;
        volScalarField rans;

	    autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;
        autoPtr<LESdelta> delta_;

    // Private member functions

        dimensionedScalar cD(const volSymmTensorField& D) const;

        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, beta1_, beta2_);
        }

        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, gamma1_, gamma2_);
        }


public:

    //- Runtime type information
    TypeName("hybKOmegaSST");


    // Constructors

        //- Construct from components
        hybKOmegaSST
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
   	        const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~hybKOmegaSST()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK(F1)*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega(F1)*nut_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }
        
        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff( volVectorField& U ) const;

        virtual tmp<fvVectorMatrix> divDevRhoReff
        ( 
            const volScalarField& rho, 
            volVectorField& U 
        ) const ;
       
        
        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
