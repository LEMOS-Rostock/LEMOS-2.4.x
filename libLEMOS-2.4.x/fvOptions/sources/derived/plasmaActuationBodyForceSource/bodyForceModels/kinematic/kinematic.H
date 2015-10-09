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
    Foam::kinematic

Description
    Fixed trim coefficients

SourceFiles
    kinematic.C

\*---------------------------------------------------------------------------*/

#ifndef kinematic_H
#define kinematic_H

#include "bodyForceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class kinematic Declaration
\*---------------------------------------------------------------------------*/

class kinematic
:
    public bodyForceModel
{

private:
    
    //- Drag coefficients along the coordinate directions
    vector CD_;

    //- Determine axial length of cellset
    scalar calcAxialLength() const;

    //- Compute body force
    tmp<vectorField> computeSup(fvMatrix<vector>& eqn);


public:

    //- Run-time type information
    TypeName("kinematic");

    //- Constructor
    kinematic(const fv::plasmaActuationBodyForce& dbd, 
                      const dictionary& dict,
                      const fvMesh& mesh,
                      const labelList& cells
                     );

    //- Destructor
    virtual ~kinematic();


    // Member functions

        //- Read
        bool read(const dictionary& dict);
/*
        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        );

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        );
*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
