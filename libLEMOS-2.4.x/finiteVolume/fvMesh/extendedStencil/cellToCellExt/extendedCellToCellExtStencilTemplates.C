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

\*---------------------------------------------------------------------------*/

#include "extendedCellToCellExtStencil.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type, class TransformOp>
void Foam::extendedCellToCellExtStencil::collectData
(
    const mapDistribute& map,
    const labelListList& untransformedElements,
    const labelListList& transformedElements,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    List<List<Type> >& stencilFld,
    const TransformOp& top
)
{
    // 1. Construct cell data in compact addressing
    List<Type> flatFld(map.constructSize(), pTraits<Type>::zero);
 
    // Insert my internal values
    forAll(fld, cellI)
    {
 	flatFld[cellI] = fld[cellI];
    }
    // Insert my boundary values
    forAll(fld.boundaryField(), patchI)
    {
 	const fvPatchField<Type>& pfld = fld.boundaryField()[patchI];
 	
 	label nCompact = 
	    pfld.patch().start()
            -fld.mesh().nInternalFaces()
 	    +fld.mesh().nCells();
 
 	forAll(pfld, i)
 	{
 	    flatFld[nCompact++] = pfld[i];
 	}
    }
 
    // Do all swapping
    map.distribute(fld.mesh().globalData().globalTransforms(), flatFld, top);
    
    // 2. Pull to stencil
    stencilFld.setSize(untransformedElements.size());
    
    forAll(untransformedElements, elementI)
    {
    	const labelList& untransformedCompactCells = untransformedElements[elementI];
    	const labelList& transformedCompactCells = transformedElements[elementI];
    
        stencilFld[elementI].setSize(untransformedCompactCells.size()+transformedCompactCells.size());
    
	label nCompact = 0;
    	forAll(untransformedCompactCells, i)
        {
    	    stencilFld[elementI][nCompact++] = flatFld[untransformedCompactCells[i]];
        }
    	forAll(transformedCompactCells, i)
        {
    	    stencilFld[elementI][nCompact++] = flatFld[transformedCompactCells[i]];
        }
    }
}





template<class Type, class WeightType, class TransformOp>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<WeightType, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
> Foam::extendedCellToCellExtStencil::weightedSum
(
    const mapDistribute& map,
    const labelListList& untransformedElements,
    const labelListList& transformedElements,
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const List<List<WeightType> >& stencilWeights,
    const TransformOp& top
)
{
    typedef typename outerProduct<WeightType, Type>::type WeightedType;
    typedef GeometricField<WeightedType, fvPatchField, volMesh>
        WeightedFieldType;

    // Collect internal and boundary values
    List<List<Type> > stencilFld;
    collectData(map, untransformedElements, transformedElements, fld, stencilFld, top);

    const fvMesh& mesh = fld.mesh();

    tmp<WeightedFieldType> twf
    (
        new WeightedFieldType
        (
            IOobject
            (
                fld.name(),
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<WeightedType>
            (
                fld.name(),
                fld.dimensions(),
                pTraits<WeightedType>::zero
            ),
	    zeroGradientFvPatchField<WeightedType>::typeName
        )
    );


    WeightedFieldType& wf = twf();


    forAll(wf, celli)
    {
        const List<Type>& stField = stencilFld[celli];
        const List<WeightType>& stWeight = stencilWeights[celli];

        forAll(stField, i)
        {
            wf[celli] += stWeight[i]*stField[i];
        }
    }

    // Boundaries values set to zeroGradient
    
    return twf;
}


// ************************************************************************* //
