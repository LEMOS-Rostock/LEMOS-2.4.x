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
    Foam::extendedCellToCellExtStencil

Description
    Baseclass for cell-to-cell stencils

SourceFiles
    extendedCellToCellExtStencil.C
    extendedCellToCellExtStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCellToCellExtStencil_H
#define extendedCellToCellExtStencil_H

#include "mapDistribute.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedCellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCellToCellExtStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;

public:

    // Constructors

        //- Construct from mesh
        explicit extendedCellToCellExtStencil(const polyMesh&);

        //- Write some statistics about stencil
        static void writeStencilStats
        (
            const labelListList& untransformedElements,
            const labelListList& transformedElements,
            const mapDistribute& map
        );

    // Member Functions
    
	//- Use map to get the (non-position) data into stencil order
	template<class Type, class TransformOp>
	static void collectData
	(
		const mapDistribute& map,
		const labelListList& untransformedElements,
		const labelListList& transformedElements,
		const GeometricField<Type, fvPatchField, volMesh>& fld,
		List<List<Type> >& stencilFld,
		const TransformOp& top
	);
	
        //- Sum surface field contributions to create cell values
        template<class Type, class WeightType, class TransformOp>
        static
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
            const mapDistribute& map,
	    const labelListList& untransformedElements,
	    const labelListList& transformedElements,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<WeightType> >& stencilWeights,
            const TransformOp& top
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedCellToCellExtStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
