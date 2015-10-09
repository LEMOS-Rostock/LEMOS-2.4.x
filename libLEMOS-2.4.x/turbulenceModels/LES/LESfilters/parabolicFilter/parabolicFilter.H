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
    Foam::parabolicFilter

Description
    Parabolic filter using spatiotemporal information for filtering. 
 
SourceFiles
    parabolicFilter.C

Description
    Parabolic filter for LES

    \verbatim
    Kernel                 as filter          as Test filter with ratio 2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Box filter:            g = delta2/24  ->  g = delta2/6
    Spherical box filter:  g = delta2/64  ->  g = delta2/16
    Gaussian filter:       g = delta2/24  ->  g = delta2/6
    \endverbatim

SourceFiles
    parabolicFilter.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicFilter_H
#define parabolicFilter_H

#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parabolicFilter Declaration
\*---------------------------------------------------------------------------*/

class parabolicFilter
:
    public LESfilter
{
    // Private data

        scalar widthCoeff_;
        scalar theta_;
        volScalarField coeff_;
        volScalarField thetaCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        parabolicFilter(const parabolicFilter&);
        void operator=(const parabolicFilter&);


public:

    //- Runtime type information
    TypeName("parabolic");


    // Constructors

        //- Construct from components
        parabolicFilter(const fvMesh& mesh, scalar widthCoeff, scalar theta);

        //- Construct from IOdictionary
        parabolicFilter(const fvMesh& mesh, const dictionary&);


    //- Destructor
    virtual ~parabolicFilter()
    {}


    // Member Functions

        //- Read the LESfilter dictionary
        virtual void read(const dictionary&);


    // Member Operators

        virtual tmp<volScalarField> operator()
        (
            const tmp<volScalarField>&
        ) const;

        virtual tmp<volVectorField> operator()
        (
            const tmp<volVectorField>&
        ) const;

        virtual tmp<volSymmTensorField> operator()
        (
            const tmp<volSymmTensorField>&
        ) const;

        virtual tmp<volTensorField> operator()
        (
            const tmp<volTensorField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
