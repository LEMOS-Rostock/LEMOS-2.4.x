/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 LEMOS, University Rostock
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
    Foam::weightedArithmeticMeanFilter

Description
    Simple weighted arithmetic mean filter using point connected cell stencil. 
    Central stencil element is weighted by alpha, all other stencil elements 
    by beta.
 
    \verbatim

    mean                   alpha        beta            Definition
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    arithmetic mean:       1.0          1.0             (x_1+...+x_n)/n

    weighted
    arithmetic mean:       > 1.0        > 1.0           (alpha*x_1+
                                                         beta*(x_2+...+x_n)
                                                        )/(alpha+beta*(n-1))
   
    \endverbatim

SourceFiles
    weightedArithmeticMeanFilter.C

\*---------------------------------------------------------------------------*/

#ifndef weightedArithmeticMeanFilter_H
#define weightedArithmeticMeanFilter_H

#include "LESfilter.H"
#include "centredCPCCellToCellExtStencilObject.H"
#include "volFields.H"
#include "scalarMatrices.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class legendreFilter Declaration
\*---------------------------------------------------------------------------*/

class weightedArithmeticMeanFilter
:
    public LESfilter

{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        weightedArithmeticMeanFilter(const weightedArithmeticMeanFilter&);
        void operator=(const weightedArithmeticMeanFilter&);

        const extendedCentredCellToCellExtStencil& addressing_;

	//- Weight alpha for central stencil element
	scalar alpha_;

	//- Weight beta for all other stencil elements
	scalar beta_;

	//- Stencil weights 
	List<List<scalar> > stencilWeights_;

public:

    //- Runtime type information
    TypeName("weightedArithmeticMean");


    // Constructors

        //- from components
        weightedArithmeticMeanFilter(const fvMesh& mesh, scalar alpha, scalar beta);

        //- from IOdictionary
        weightedArithmeticMeanFilter(const fvMesh& mesh, const dictionary& dict);


    // Destructor
    virtual ~weightedArithmeticMeanFilter()
    {}



    // Member Functions

        //- Read the LESfilter dictionary
        void read(const dictionary& dict);

        //- Arrange element weights in stencil order 
        void stencilWeights(scalar alpha, scalar beta);

        //- Update stencil and dependencies after mesh changing
        void updateStencil();

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
