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
    Foam::cellToCellExtStencil

Description
    baseclass for extended cell centred addressing. Contains per cell a
    list of neighbouring cells and/or boundaryfaces in global addressing. 
    Supports also transformations.

SourceFiles
    cellToCellExtStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellToCellExtStencil_H
#define cellToCellExtStencil_H

#include "globalIndex.H"
#include "globalIndexAndTransform.H"
#include "boolList.H"
#include "HashSet.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class cellToCellExtStencil Declaration
\*---------------------------------------------------------------------------*/

class cellToCellExtStencil
:
    public labelListList
{
    // Private data

        const polyMesh& mesh_;

        //- Global numbering for cells and boundary faces
        const globalIndex globalNumbering_;

        //- Global numbering for transformed elements
        const globalIndexAndTransform globalTransforms_;

protected:

	//- Untransformed elements
	labelListList untransformedElements_;
	       
	//- Transformed elements
	List<labelPairList> transformedElements_;

        //- Merge two lists.
        static void merge
        (
            const label,
            const label,
            const labelList&,
            labelList&
        );

        //- Merge two lists.
        static void merge(const label, const labelList&, labelList&);

        //- Merge two lists.
        static void merge(const labelPairList&, labelPairList&);

        //- Valid boundary faces (not empty)
        void validBoundaryFaces(boolList& isValidBFace) const;

        //- Return patch of all coupled faces.
        autoPtr<indirectPrimitivePatch> allCoupledFacesPatch() const;

        //- Return patch of coupled faces.
        autoPtr<indirectPrimitivePatch> coupledFacesPatch(const polyPatch& pp) const;

        //- Collect cell neighbours of faces in global numbering
        void insertFaceCells
        (
            const label exclude0,
            const label exclude1,
	    const boolList& nonEmptyFace,
	    const labelList& faceLabels,
	    labelHashSet& globals
	) const;
 
 	//- Collect cell neighbours of faces in global numbering
 	labelList calcFaceCells
	(
	    const boolList& nonEmptyFace,
	    const labelList& faceLabels,
	    labelHashSet& globals
	) const;
 

public:

    // Constructors

        //- Construct from mesh
        explicit cellToCellExtStencil(const polyMesh&);


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Global numbering for cells and boundary faces
        const globalIndex& globalNumbering() const
        {
            return globalNumbering_;
        }

        //- Global numbering for transformed elements
        const globalIndexAndTransform& globalTransforms() const
        {
            return globalTransforms_;
        }
        
        //- Untransformed elements
        const labelListList& untransformedElements() const
        {
	    return untransformedElements_;
        }
        
        //- Untransformed elements
        labelListList& untransformedElements()
        {
            return untransformedElements_;
        }
        
        //- Transformed elements
        const List<labelPairList>& transformedElements() const
        {
            return transformedElements_;
        }
        
        //- Transformed elements
        List<labelPairList>& transformedElements()
        {
            return transformedElements_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
