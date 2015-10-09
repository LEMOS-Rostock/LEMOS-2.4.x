/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University of Rostock
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Class
    simpleGridFilter

Description
    Simple top-hat filter used in dynamic LES models. Implemented as a surface 
    integral of the face interpolate of the field using a filter width equal to
    the grid size.

SourceFiles
    simpleGridFilter.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGridFilter_H
#define simpleGridFilter_H

#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleGridFilter Declaration
\*---------------------------------------------------------------------------*/

class simpleGridFilter
:
    public LESfilter

{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        simpleGridFilter(const simpleGridFilter&);
        void operator=(const simpleGridFilter&);


public:

    //- Runtime type information
    TypeName("simpleGridFilter");


    // Constructors

        //- from components
        simpleGridFilter(const fvMesh& mesh);

        //- from IOdictionary
        simpleGridFilter(const fvMesh& mesh, const dictionary&);


    // Destructor

        ~simpleGridFilter()
        {}


    // Member Functions

        //- Read the LESfilter dictionary
        void read(const dictionary&);


    // Member Operators

        tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
