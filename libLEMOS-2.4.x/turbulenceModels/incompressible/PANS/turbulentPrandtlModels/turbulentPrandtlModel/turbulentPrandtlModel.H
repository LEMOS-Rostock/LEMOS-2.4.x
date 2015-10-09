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


Namespace
    Foam::incompressible::turbulentPrandtlModels:turbulentPrandtlModels

Group
    grpIcoPANSModels

Description
    Namespace for incompressible turbulent Prandtl models.

Class
    Foam::incompressible::turbulentPrandtlModel

Description
    Abstract base class for incompressible turbulent Prandtl models.

SourceFiles
    turbulentPrandtlModel.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentPrandtlModel_H
#define turbulentPrandtlModel_H

#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace PANSModels
{

/*---------------------------------------------------------------------------*\
                           Class turbulentPrandtlModel Declaration
\*---------------------------------------------------------------------------*/

class turbulentPrandtlModel
:
    public IOdictionary
{

protected:

    // Protected data

    const fvMesh& mesh_;

    word param_;

    dictionary coeffDict_;

    dimensionedScalar sigmaEps_;
    dimensionedScalar sigmaK_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        turbulentPrandtlModel(const turbulentPrandtlModel&);

        //- Disallow default bitwise assignment
        void operator=(const turbulentPrandtlModel&);


public:

    //- Runtime type information
    TypeName("turbulentPrandtlModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            turbulentPrandtlModel,
            dictionary,
            (   
                const volVectorField& U,
                const surfaceScalarField& phi,
                const word& param,
                const word& turbulenceModelName
            ),
            (U, phi, param, turbulenceModelName)
        );


    // Constructors

        //- Construct from components
        turbulentPrandtlModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& param,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected turbulentPrandtl model
        static autoPtr<turbulentPrandtlModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& param,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~turbulentPrandtlModel()
    {}


    // Member Functions

        // Access

       


        //- 
        virtual tmp<volScalarField> sigmaEps() const = 0;
        
        //-
        virtual tmp<volScalarField> sigmaK() const = 0;
        
        //-
        virtual tmp<volScalarField> sigmaOmega() const = 0;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read turbulentPrandtlProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
