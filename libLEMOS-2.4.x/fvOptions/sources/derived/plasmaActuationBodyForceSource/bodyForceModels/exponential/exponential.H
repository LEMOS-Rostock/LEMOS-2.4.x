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
    Foam::exponential

Description
    Fixed trim coefficients

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "bodyForceModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public bodyForceModel
{

private:
    
    //- Drag coefficients along the coordinate directions
    vector D_;

    //- Constant 1
    List<scalar> Cx_;

    //- Constant 2
    List<scalar> Cy_;

    //- Constant 3
    List<scalar> Cz_;

    //- Maximal force
    vector fMax_;

    //- Body force field
    tmp<vectorField> Su_;

    //- Correct body force field
    bool correct_;

    //- Determine axial length of cellset
    void calcForceField();

    //- Compute body force
    tmp<vectorField> computeSup(fvMatrix<vector>& eqn);


public:

    //- Run-time type information
    TypeName("exponential");

    //- Constructor
    exponential(const fv::plasmaActuationBodyForce& dbd, 
                      const dictionary& dict,
                      const fvMesh& mesh,
                      const labelList& cells
                     );

    //- Destructor
    virtual ~exponential();


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
