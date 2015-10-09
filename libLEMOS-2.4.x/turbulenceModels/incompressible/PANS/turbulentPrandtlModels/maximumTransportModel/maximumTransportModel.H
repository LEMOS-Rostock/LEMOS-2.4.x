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
    Foam::incompressible::maximumTransportModels:maximumTransportModels

Group
    grpIcoPANSModels

Description
    Namespace for incompressible turbulent Prandtl models.

Class
    Foam::incompressible::maximumTransportModel

Description
    Zero transport model

SourceFiles
    maximumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef maximumTransportModel_H
#define maximumTransportModel_H

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
#include "turbulentPrandtlModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace PANSModels
{
namespace turbulentPrandtlModels
{
/*---------------------------------------------------------------------------*\
                           Class maximumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class maximumTransportModel
:
    public turbulentPrandtlModel
{
private:


    // Private Member Functions

        //- Disallow default bitwise copy construct
        maximumTransportModel(const maximumTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const maximumTransportModel&);


public:

    //- Runtime type information
    TypeName("maximumTransportModel");


    // Constructors

        //- Construct from components
        maximumTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& param,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~maximumTransportModel()
    {}


    // Member Functions

        // Access

        //- 
        virtual tmp<volScalarField> sigmaEps() const ;
        
        //-
        virtual tmp<volScalarField> sigmaK() const;

        //-
        virtual tmp<volScalarField> sigmaOmega() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read turbulentPrandtlProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulentPrandtlModels
} // End namespace PANSModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
