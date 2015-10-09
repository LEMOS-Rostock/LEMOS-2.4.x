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

#include "CPCCellToCellExtStencil.H"
#include "syncTools.H"
#include "dummyTransform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculates per point the coupled neighbour data (= pointCells)
void Foam::CPCCellToCellExtStencil::calcPointBoundaryData
(
    const boolList& isValidBFace,
    Map<labelPairList>& neiGlobal
) const
{
    const labelList& own = mesh().faceOwner();
    const labelList& nei = mesh().faceNeighbour();

    neiGlobal.resize(mesh().nPoints()-mesh().nInternalPoints());

    const polyBoundaryMesh& patches = mesh().boundaryMesh();

    // Loop over all patches
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
	    labelList boundaryPoints(coupledFacesPatch(pp)().meshPoints());

	    forAll(boundaryPoints, i)
	    {
		label pointI =  boundaryPoints[i];
		
	        labelList pointFaces = mesh().pointFaces()[pointI];

    		Map<labelPair> lGlobals;

		forAll(pointFaces, j)
		{
		    label faceI = pointFaces[j];

		    lGlobals.insert
            	    (	
                        own[faceI],
                        globalIndexAndTransform::encode
                        (
                            Pstream::myProcNo(),
                            own[faceI],
                            globalTransforms().addToTransformIndex
                            (
                                globalTransforms().nullTransformIndex(),
                                patchI,
                                true
                            )
                        )
                    );

		    if (mesh().isInternalFace(faceI))
        	    {

                        lGlobals.insert
                        (
                	    nei[faceI],
                	    globalIndexAndTransform::encode
                	    (
                	        Pstream::myProcNo(),
                		nei[faceI],
                		globalTransforms().addToTransformIndex
                		(
                		    globalTransforms().nullTransformIndex(),
				    patchI,
				    true
			        )
                	    )
                  	);

		    }
		    else
		    {
		        label bFaceI = faceI-mesh().nInternalFaces();
		
		        if (isValidBFace[bFaceI])
			{
			    lGlobals.insert
                    	    (
                                bFaceI+mesh().nCells(),
                        	globalIndexAndTransform::encode
                        	(
                        	    Pstream::myProcNo(),
                        	    bFaceI+mesh().nCells(),
                        	    globalTransforms().addToTransformIndex
                        	    (
                        	        globalTransforms().nullTransformIndex(),
                        		patchI,
                        		true
                        	    )    
                        	)
                            );
			}
		    }
                }

 
    		labelPairList pGlobals;

	        forAllIter(Map<labelPair>,lGlobals, iter)
                {
                    pGlobals.append(iter());
                }
    		
		neiGlobal.insert(pointI, pGlobals);

		lGlobals.clearStorage();
		pGlobals.clear();

	    }
        }
    }

    syncTools::syncPointMap
    (
        mesh(),
        neiGlobal,
        eqOp<labelPairList>(),
        Foam::dummyTransform()      // dummy transformation
    );
}



// Calculates per cell the neighbour data (= cell or boundary in global
// numbering). First element is always cell itself!
void Foam::CPCCellToCellExtStencil::calcCellStencil
(
    labelListList& untransformedElements,
    List<labelPairList>& transformedElements
) const
{
    // Calculate points on coupled patches
    labelList boundaryPoints(allCoupledFacesPatch()().meshPoints());

    // Mark boundary faces to be included in stencil (i.e. not empty)
    boolList isValidBFace;
    validBoundaryFaces(isValidBFace);

    // Swap pointCells for coupled points
    Map<labelPairList> neiGlobal;
    calcPointBoundaryData
    (
        isValidBFace,
        neiGlobal
    );

    untransformedElements.setSize(mesh().nCells());
    transformedElements.setSize(mesh().nCells());

    // Do coupled points first
    forAll(boundaryPoints, i)
    {
        label pointI = boundaryPoints[i];

        const labelPairList& pGlobals = neiGlobal[pointI];

        // Distribute to all pointCells
        const labelList& pCells = mesh().pointCells(pointI);


        // DynamicList<label> untrafoFaces(ownFaces.size());
        DynamicList<labelPair> trafoElements(pGlobals.size());

        forAll(pCells, jj)
	{
            label cellI = pCells[jj];
	    forAll(pGlobals, j)
            {

	        const labelPair& info = pGlobals[j];
                label transform =
                globalIndexAndTransform::transformIndex
                (
                    info
                );

	        if
                (
                    transform
                 == globalTransforms().nullTransformIndex()
                )
                {
                    label procI =
                        globalIndexAndTransform::processor(info);
                    label index =
                        globalIndexAndTransform::index(info);

                    untransformedElements[cellI].append
                    (
                        globalNumbering().toGlobal(procI, index)
                    );

                }
                else
                {
                    trafoElements.append(info);
                }
            }

	    merge
	    (
	        trafoElements,
	        transformedElements[cellI]
	    );
        }
    }

    neiGlobal.clear();


    // Do remaining points cells
    labelHashSet pointGlobals;

    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    {
        labelList pGlobals
        (
            calcFaceCells
            (
                isValidBFace,
                mesh().pointFaces()[pointI],
                pointGlobals
            )
        );

        const labelList& pCells = mesh().pointCells(pointI);

        forAll(pCells, j)
        {
            label cellI = pCells[j];

            merge
            (
                globalNumbering().toGlobal(cellI),
                pGlobals,
                untransformedElements[cellI]
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPCCellToCellExtStencil::CPCCellToCellExtStencil(const polyMesh& mesh)
:
    cellToCellExtStencil(mesh)
{
    // Calculate per cell the (point) connected cells (in global numbering)
    calcCellStencil(untransformedElements_, transformedElements_);
}


// ************************************************************************* //
