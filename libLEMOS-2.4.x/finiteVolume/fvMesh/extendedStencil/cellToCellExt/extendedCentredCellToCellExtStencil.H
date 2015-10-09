/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
 2013-11-05 LEMOS, University of Rostock: added support for transformations
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
    Foam::extendedCentredCellToCellExtStencil

Description

SourceFiles
    extendedCentredCellToCellExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCentredCellToCellExtStencil_H
#define extendedCentredCellToCellExtStencil_H

#include "extendedCellToCellExtStencil.H"
//#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellToCellExtStencil;
class mapDistribute;

/*---------------------------------------------------------------------------*\
              Class extendedCentredCellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCentredCellToCellExtStencil
:
    public extendedCellToCellExtStencil
{
    // Private data

        //- Swap map for getting neigbouring data
        autoPtr<mapDistribute> mapPtr_;

        //- Per cell the untransformed elements
        labelListList untransformedElements_;

        //- Per cell the transformed elements
        labelListList transformedElements_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedCentredCellToCellExtStencil
        (
            const extendedCentredCellToCellExtStencil&
        );

        //- Disallow default bitwise assignment
        void operator=(const extendedCentredCellToCellExtStencil&);


public:

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit extendedCentredCellToCellExtStencil(const cellToCellExtStencil&);


    // Member Functions

        //- Return reference to the parallel distribution map
        const mapDistribute& map() const
        {
            return mapPtr_();
        }

        //- Return reference to untransformed elements in the stencil
        const labelListList& untransformedElements() const
        {
            return untransformedElements_;
        }

        //- Return reference to transformed elements in the stencil
        const labelListList& transformedElements() const
        {
            return transformedElements_;
        }

        //- After removing elements from the stencil adapt the schedule (map).
        void compact();

        //- Use map to get the data into stencil order
        template<class Type>
        void collectData
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            List<List<Type> >& stencilFld
        ) const
        {
            extendedCellToCellExtStencil::collectData
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld,
                stencilFld,
		mapDistribute::transform()
            );
        }

        //- Use map to get the positions into stencil order
        template<class Type>
        void collectPositions
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            List<List<Type> >& stencilFld
        ) const
        {
            extendedCellToCellExtStencil::collectData
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld,
                stencilFld,
		mapDistribute::transformPosition()
            );
        }

        //- Sum vol field contributions to create cell values
        template<class Type, class WeightType>
        tmp
        <
            GeometricField
            <
                typename outerProduct<WeightType, Type>::type,
                fvPatchField,
                volMesh
            >
        > weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<WeightType> >& stencilWeights
        ) const
        {
            return extendedCellToCellExtStencil::weightedSum
            (
                map(),
                untransformedElements(),
                transformedElements(),
                fld,
                stencilWeights,
                mapDistribute::transform()
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
