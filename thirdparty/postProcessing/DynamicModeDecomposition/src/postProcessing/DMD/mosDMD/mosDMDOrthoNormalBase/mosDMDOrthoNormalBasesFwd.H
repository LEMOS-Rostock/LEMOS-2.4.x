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
    DMDOrthoNormalBase

Description
    DMD base typedefs

\*---------------------------------------------------------------------------*/

#ifndef mosDMDOrthoNormalBasesFwd_H
#define mosDMDOrthoNormalBasesFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type> class mosDMDOrthoNormalBase;

typedef mosDMDOrthoNormalBase<scalar> scalarMosDMDOrthoNormalBase;
typedef mosDMDOrthoNormalBase<vector> vectorMosDMDOrthoNormalBase;
typedef mosDMDOrthoNormalBase<tensor> tensorMosDMDOrthoNormalBase;
typedef mosDMDOrthoNormalBase<symmTensor> symmTensorMosDMDOrthoNormalBase;
typedef mosDMDOrthoNormalBase<sphericalTensor> sphericalTensorMosDMDOrthoNormalBase;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
