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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCellToCellExtStencil::extendedCellToCellExtStencil(const polyMesh& mesh)
:
    mesh_(mesh)
{
}

void Foam::extendedCellToCellExtStencil::writeStencilStats
(
    const labelListList& untransformedElements,
    const labelListList& transformedElements,
    const mapDistribute& map
)
{
    label sumSizeUntransformedElements = 0;
    label sumSizeTransformedElements = 0;
    label nSumUntransformedElements = 0;
    label nSumTransformedElements = 0;
    label minSizeUntransformedElements = labelMax;
    label minSizeTransformedElements = labelMax;
    label maxSizeUntransformedElements = labelMin;
    label maxSizeTransformedElements = labelMin;

    label minSize = labelMax;
    label maxSize = labelMin;
    label sumSize = 0;
    label nSum = 0;

    forAll(untransformedElements, i)
    {
        const labelList& sCells = untransformedElements[i];

        if (sCells.size() > 0)
        {
            sumSizeUntransformedElements += sCells.size();
            nSumUntransformedElements++;
            minSizeUntransformedElements = min(minSizeUntransformedElements, sCells.size());
            maxSizeUntransformedElements = max(maxSizeUntransformedElements, sCells.size());
        }
    }

    forAll(transformedElements, i)
    {
        const labelList& sCells = transformedElements[i];

        if (sCells.size() > 0)
        {
            sumSizeTransformedElements += sCells.size();
            nSumTransformedElements++;
            minSizeTransformedElements = min(minSizeTransformedElements, sCells.size());
            maxSizeTransformedElements = max(minSizeTransformedElements, sCells.size());
        }
    }

    minSize = minSizeUntransformedElements + minSizeTransformedElements;
    maxSize = maxSizeUntransformedElements + maxSizeTransformedElements;

    sumSize = sumSizeUntransformedElements + sumSizeTransformedElements;
    nSum = nSumUntransformedElements + nSumTransformedElements;

    reduce(sumSize, sumOp<label>());
    reduce(nSum, sumOp<label>());

    reduce(minSize, minOp<label>());
    reduce(maxSize, maxOp<label>());

    Info<< "Stencil size :" << nl
        << "    average : " << sumSize << " " <<nSum << nl
        << "    average : " << scalar(sumSize)/nSum << nl
        << "    min     : " << minSize << nl
        << "    max     : " << maxSize << nl
        << endl;

    label nSent = 0;
    label nLocal = 0;
    forAll(map.subMap(), procI)
    {
        if (procI != Pstream::myProcNo())
        {
            nSent += map.subMap()[procI].size();
        }
        else
        {
            nLocal += map.subMap()[procI].size();
        }
    }

    Info<< "Local data size : " << returnReduce(nLocal, sumOp<label>()) << nl
        << "Sent data size  : " << returnReduce(nSent, sumOp<label>()) << nl
        << endl;
}


// ************************************************************************* //
