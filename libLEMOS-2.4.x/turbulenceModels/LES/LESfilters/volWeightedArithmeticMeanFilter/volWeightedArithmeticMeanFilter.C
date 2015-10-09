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

\*---------------------------------------------------------------------------*/

#include "volWeightedArithmeticMeanFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volWeightedArithmeticMeanFilter, 0);
addToRunTimeSelectionTable(LESfilter, volWeightedArithmeticMeanFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
volWeightedArithmeticMeanFilter::volWeightedArithmeticMeanFilter
(
    const fvMesh& mesh
)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh))
{
    stencilWeights();
}


volWeightedArithmeticMeanFilter::volWeightedArithmeticMeanFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh))
{
    stencilWeights();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volWeightedArithmeticMeanFilter::read(const dictionary& dict)
{
    stencilWeights();
}

void volWeightedArithmeticMeanFilter::stencilWeights()
{

    const labelListList& untransformedElements = addressing_.untransformedElements();
    const labelListList& transformedElements = addressing_.transformedElements();

    // Determine volumes of stencil elements
    // Note: boundary faces are handled by using volume 
    //       of (face connected) neighbour cells.
    //       At the moment, no extrapolation is done but
    //       boundary face value as virtual cell value
    tmp<volScalarField> tCellVolumes
    (
        new volScalarField
        (
            IOobject
            (
                "cellVolumes",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            0.0,
            zeroGradientFvPatchField<scalar>::typeName
        )
    );

    volScalarField& cellVolumes = tCellVolumes();

    cellVolumes.internalField()=mesh().V();
    cellVolumes.correctBoundaryConditions();

    List<List<scalar> > stencilVolumes;

    addressing_.collectData
    (
        cellVolumes,
	stencilVolumes
    );

    stencilWeights_.resize(mesh().nCells());

    forAll(stencilWeights_, stencilWeightsI)
    {
	label nElements = untransformedElements[stencilWeightsI].size() + transformedElements[stencilWeightsI].size();

        List<scalar>& stencilWeight = stencilWeights_[stencilWeightsI];
	stencilWeight.clear();

        List<scalar>& stencilVolume = stencilVolumes[stencilWeightsI];

	if(nElements > 0)
	{
	    scalar sumStencilVolume = 0.0;

	    forAll(stencilVolume, stencilVolumeI)
	    {
	        sumStencilVolume += stencilVolume[stencilVolumeI];
            }	

	    forAll(stencilVolume, stencilVolumeI)
	    {
	        stencilWeight.append(stencilVolume[stencilVolumeI]/sumStencilVolume);
	    }
	}
    }
}

void volWeightedArithmeticMeanFilter::updateStencil()
{
    // ToDo: recreate stencil after mesh changing
    //addressing_ = (centredCPCCellToCellExtStencilObject::New(mesh()));

    //stencilWeights(alpha_, beta_);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


Foam::tmp<Foam::volVectorField> Foam::volWeightedArithmeticMeanFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const 
{
    const tmp<volVectorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volScalarField> Foam::volWeightedArithmeticMeanFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    const tmp<volScalarField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volSymmTensorField> Foam::volWeightedArithmeticMeanFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    const tmp<volSymmTensorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volTensorField> Foam::volWeightedArithmeticMeanFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    const tmp<volTensorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();


    return filteredField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
