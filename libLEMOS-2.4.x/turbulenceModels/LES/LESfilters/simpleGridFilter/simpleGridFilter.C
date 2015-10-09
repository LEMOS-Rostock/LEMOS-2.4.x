/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University of Rostock
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "simpleGridFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleGridFilter, 0);
addToRunTimeSelectionTable(LESfilter, simpleGridFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
simpleGridFilter::simpleGridFilter
(
    const fvMesh& mesh
)
:
    LESfilter(mesh)
{}


simpleGridFilter::simpleGridFilter(const fvMesh& mesh, const dictionary&)
:
    LESfilter(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simpleGridFilter::read(const dictionary&)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

tmp<volScalarField> simpleGridFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    const volScalarField& uf=unFilteredField();

    surfaceScalarField surfvals=fvc::interpolate(uf);

    tmp<volScalarField> tfilteredField
    (
        new volScalarField
        (
            IOobject
            (
                "filtered",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", uf.dimensions(), pTraits<scalar>::zero)
        )
    );

    volScalarField& filteredField=tfilteredField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        filteredField[owner[faceI]]+=
	    0.5*( surfvals[faceI] + uf[owner[faceI]] ) 
	    * mesh().magSf()[faceI];
      
        filteredField[neighbour[faceI]]+=
	    0.5*( surfvals[faceI] + uf[neighbour[faceI]] ) 
	    * mesh().magSf()[faceI];
    }
  
    forAll(mesh().boundary(), patchI)
    {
        const unallocLabelList& pFaceCells =
	    mesh().boundary()[patchI].faceCells();
      
        const fvPatchScalarField& ufsf = uf.boundaryField()[patchI];

        forAll(mesh().boundary()[patchI], faceI)
        {
	    filteredField[pFaceCells[faceI]]+=
	        0.5*(ufsf[faceI]+uf[pFaceCells[faceI]]) 
	        * mesh().boundary()[patchI].magSf()[faceI];
	}
    }

    volScalarField ssum=fvc::surfaceSum(mesh().magSf());
    
    forAll(ssum, cellI)
    {
      filteredField[cellI]/=ssum[cellI];
    }
    filteredField.correctBoundaryConditions();
   
    unFilteredField.clear();
    return tfilteredField;
}


tmp<volVectorField> simpleGridFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    const volVectorField& uf=unFilteredField();

    surfaceVectorField surfvals=fvc::interpolate(uf);

    tmp<volVectorField> tfilteredField
    (
        new volVectorField
        (
            IOobject
            (
                "filtered",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector("zero", uf.dimensions(), pTraits<vector>::zero)
        )
    );

    volVectorField& filteredField=tfilteredField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        filteredField[owner[faceI]]+=
	    0.5*( surfvals[faceI] + uf[owner[faceI]] ) 
            * mesh().magSf()[faceI];
      
        filteredField[neighbour[faceI]]+=
	    0.5*( surfvals[faceI] + uf[neighbour[faceI]] ) 
            * mesh().magSf()[faceI];
    }
  
    forAll(mesh().boundary(), patchI)
    {
        const unallocLabelList& pFaceCells =
	    mesh().boundary()[patchI].faceCells();
      
        const fvPatchVectorField& ufsf = uf.boundaryField()[patchI];

        forAll(mesh().boundary()[patchI], faceI)
        {
	    filteredField[pFaceCells[faceI]]+=
	        0.5*(ufsf[faceI]+uf[pFaceCells[faceI]]) 
                * mesh().boundary()[patchI].magSf()[faceI];
	}
    }

    volScalarField ssum=fvc::surfaceSum(mesh().magSf());
    
    forAll(ssum, cellI)
    {
      filteredField[cellI]/=ssum[cellI];
    }
    filteredField.correctBoundaryConditions();

    unFilteredField.clear();

    return tfilteredField;
}


tmp<volSymmTensorField> simpleGridFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    const volSymmTensorField& uf=unFilteredField();

    surfaceSymmTensorField surfvals=fvc::interpolate(uf);

    tmp<volSymmTensorField> tfilteredField
    (
        new volSymmTensorField
        (
            IOobject
            (
                "filtered",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedSymmTensor("zero", uf.dimensions(), pTraits<symmTensor>::zero)
        )
    );

    volSymmTensorField& filteredField=tfilteredField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        filteredField[owner[faceI]]+=
	    0.5*( surfvals[faceI] + uf[owner[faceI]] ) 
            * mesh().magSf()[faceI];
      
        filteredField[neighbour[faceI]]+=
	    0.5*( surfvals[faceI] + uf[neighbour[faceI]] ) 
            * mesh().magSf()[faceI];
    }
  
    forAll(mesh().boundary(), patchI)
    {
        const unallocLabelList pFaceCells =
	    mesh().boundary()[patchI].faceCells();
      
        const fvPatchSymmTensorField& ufsf = uf.boundaryField()[patchI];

        forAll(mesh().boundary()[patchI], faceI)
        {
	    filteredField[pFaceCells[faceI]]+=
	        0.5*(ufsf[faceI]+uf[pFaceCells[faceI]]) 
                * mesh().boundary()[patchI].magSf()[faceI];
	}
    }

    volScalarField ssum=fvc::surfaceSum(mesh().magSf());
    
    forAll(ssum, cellI)
    {
      filteredField[cellI]/=ssum[cellI];
    }
    filteredField.correctBoundaryConditions();

    unFilteredField.clear();

    return tfilteredField;
}


tmp<volTensorField> simpleGridFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    const volTensorField& uf=unFilteredField();

    surfaceTensorField surfvals=fvc::interpolate(uf);

    tmp<volTensorField> tfilteredField
    (
        new volTensorField
        (
            IOobject
            (
                "filtered",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("zero", uf.dimensions(), pTraits<tensor>::zero)
        )
    );
 
    volTensorField& filteredField=tfilteredField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll(owner, faceI)
    {
        filteredField[owner[faceI]]+=
	    0.5*( surfvals[faceI] + uf[owner[faceI]] ) 
            * mesh().magSf()[faceI];
      
        filteredField[neighbour[faceI]]+=
	    0.5*( surfvals[faceI] + uf[neighbour[faceI]] ) 
            * mesh().magSf()[faceI];
    }
  
    forAll(mesh().boundary(), patchI)
    {
        const unallocLabelList pFaceCells =
	    mesh().boundary()[patchI].faceCells();
      
        const fvPatchTensorField& ufsf = uf.boundaryField()[patchI];

        forAll(mesh().boundary()[patchI], faceI)
        {
	    filteredField[pFaceCells[faceI]]+=
	        0.5*(ufsf[faceI]+uf[pFaceCells[faceI]]) 
                * mesh().boundary()[patchI].magSf()[faceI];
	}
    }

    volScalarField ssum=fvc::surfaceSum(mesh().magSf());
    
    forAll(ssum, cellI)
    {
      filteredField[cellI]/=ssum[cellI];
    }
    filteredField.correctBoundaryConditions();

    unFilteredField.clear();

    return tfilteredField;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
