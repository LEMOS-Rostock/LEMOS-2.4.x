/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::bodyForceModel

Description
    Trim model base class

SourceFiles
    bodyForceModel.C

\*---------------------------------------------------------------------------*/

#ifndef bodyForceModel_H
#define bodyForceModel_H

#include "plasmaActuationBodyForce.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class bodyForceModel Declaration
\*---------------------------------------------------------------------------*/

class bodyForceModel
{

protected:

    // Protected data

        //- Reference to the rotor source model
        const fv::plasmaActuationBodyForce& dbd_;

        //- Name of model
        const word name_;

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Reference to cellset
        const labelList& cells_;

        //- Coefficients dictionary
        dictionary coeffs_;


public:

    //- Run-time type information
    TypeName("bodyForceModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            bodyForceModel,
            dictionary,
            (
                const fv::plasmaActuationBodyForce& dbd,
                const dictionary& dict,
                const fvMesh& mesh,
                const labelList& cells
            ),
            (dbd, dict, mesh, cells)
        );


    // Constructors

        //- Construct from components
        bodyForceModel
        (
            const fv::plasmaActuationBodyForce& dbd,
            const dictionary& dict,
            const word& name,
            const fvMesh& mesh,
            const labelList& cells
        );


    // Selectors

        //- Return a reference to the selected trim model
        static autoPtr<bodyForceModel> New
        (
            const fv::plasmaActuationBodyForce& dbd,
            const dictionary& dict,
            const fvMesh& mesh,
            const labelList& cells
        );


    //- Destructor
    virtual ~bodyForceModel();


    // Member functions

        //- Read
        virtual bool read(const dictionary& dict);

        //- Compute body force
        virtual tmp<vectorField> computeSup(fvMatrix<vector>& eqn) = 0;
/*
        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        ) = 0;

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        ) = 0;
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
