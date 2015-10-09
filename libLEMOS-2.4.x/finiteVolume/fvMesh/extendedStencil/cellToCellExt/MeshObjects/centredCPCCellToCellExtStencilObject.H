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
    Foam::centredCPCCellToCellExtStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCPCCellToCellExtStencilObject_H
#define centredCPCCellToCellExtStencilObject_H

#include "extendedCentredCellToCellExtStencil.H"
#include "CPCCellToCellExtStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCPCCellToCellExtStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCPCCellToCellExtStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCPCCellToCellExtStencilObject
    >,
    public extendedCentredCellToCellExtStencil
{

public:

    TypeName("centredCPCCellToCellExtStencil");

    // Constructors

        //- Construct from uncompacted cell stencil
        explicit centredCPCCellToCellExtStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                centredCPCCellToCellExtStencilObject
            >(mesh),
            extendedCentredCellToCellExtStencil(CPCCellToCellExtStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCPCCellToCellExtStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
