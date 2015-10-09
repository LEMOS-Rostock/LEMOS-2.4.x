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


Class
    Foam::incompressible::PANSModels::kEpsilon

Group
    grpIcoPANSTurbulence

Description
    Standard k-epsilon turbulence model for incompressible flows.

    The default model coefficients correspond to the following:
    \verbatim
        kEpsilonCoeffs
        {
            Cmu         0.09;
            Ce1         1.44;
            Ce2         1.92;
            sigmaEps    1.3;
            sigmaK      1.0;
        }
    \endverbatim

SourceFiles
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilon_H
#define kEpsilon_H

#include "PANSModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace PANSModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilon Declaration
\*---------------------------------------------------------------------------*/

class kEpsilon
:
    public PANSModel
{

protected:

    // Protected data


        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;
            volScalarField fk_;
            volScalarField feps_;

    autoPtr<turbulentPrandtlModel> sigmaKModel_;
    autoPtr<turbulentPrandtlModel> sigmaEpsModel_;

public:

    //- Runtime type information
    TypeName("kEpsilon");

    // Constructors

        //- Construct from components
        kEpsilon
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kEpsilon()
    {}


    // Member Functions

        //- Return model parameter Ceps2star
        tmp<volScalarField> Ceps2Star() const;

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/(sigmaKModel_->sigmaK()) + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/(sigmaEpsModel_->sigmaEps()) + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }


        //- Return turbulent Prandtl number
        virtual tmp<volScalarField> fk() const;

        //- Return turbulent Prandtl number;
        virtual tmp<volScalarField> feps() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read PANSProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
