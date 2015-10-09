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

#include "volWeightedGaussianFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(volWeightedGaussianFilter, 0);
addToRunTimeSelectionTable(LESfilter, volWeightedGaussianFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
volWeightedGaussianFilter::volWeightedGaussianFilter
(
    const fvMesh& mesh,
    const scalar gamma
)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    gamma_(gamma)
{
    stencilWeights();
}


volWeightedGaussianFilter::volWeightedGaussianFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    gamma_(dict.subDict(type() + "Coeffs").lookupOrDefault<scalar>("gamma", 6.0))
{
    stencilWeights();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void volWeightedGaussianFilter::read(const dictionary& dict)
{
    stencilWeights();
}

void volWeightedGaussianFilter::stencilWeights()
{
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


    volScalarField sqrDeltaCoeffs = pow(cellVolumes, 2.0/3.0);
    List<List<scalar> > stencilSqrDeltaCoeffs;
    addressing_.collectData
    (
        sqrDeltaCoeffs,
        stencilSqrDeltaCoeffs
    );

    const labelListList& untransformedElements = addressing_.untransformedElements();
    const labelListList& transformedElements = addressing_.transformedElements();

    List<List<vector> > stencilCenters;

    addressing_.collectPositions
    (
        mesh().C(),
	stencilCenters
    );

    stencilWeights_.resize(mesh().nCells());

    forAll(stencilWeights_, stencilWeightsI)
    {
	label nElements = untransformedElements[stencilWeightsI].size() + transformedElements[stencilWeightsI].size();

        List<scalar>& stencilWeight = stencilWeights_[stencilWeightsI];
	stencilWeight.clear();

        List<vector>& stencilCenter = stencilCenters[stencilWeightsI];
        List<scalar>& stencilSqrDeltaCoeff = stencilSqrDeltaCoeffs[stencilWeightsI];
        List<scalar>& stencilVolume = stencilVolumes[stencilWeightsI];
 
	if(nElements > 0)
	{
            const vector& X = stencilCenter[0]; 
	    
            scalar sumGStencilVolume = 0.0;

	    forAll(stencilCenter, stencilCenterI)
	    {
	        const vector& R = stencilCenter[stencilCenterI]; 
	        const scalar sqrDeltaCoeff = stencilSqrDeltaCoeff[stencilCenterI];
	        const scalar volume = stencilVolume[stencilCenterI];

		scalar magXR = magSqr(X - R); 

		scalar G = pow((gamma_/(constant::mathematical::pi*sqrDeltaCoeff)), 1.5) * exp(-1.0*gamma_/sqrDeltaCoeff * magXR);
		stencilWeight.append(G*volume);

                sumGStencilVolume += G*volume;
            }	

            forAll(stencilWeight, stencilWeightI)
            {
	        stencilWeight[stencilWeightI] /= sumGStencilVolume;
            }
	}
    }
}

void volWeightedGaussianFilter::updateStencil()
{
    // ToDo: recreate stencil after mesh changing
    //addressing_ = (centredCPCCellToCellExtStencilObject::New(mesh()));

    //stencilWeights(alpha_, beta_);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


Foam::tmp<Foam::volVectorField> Foam::volWeightedGaussianFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const 
{
    const tmp<volVectorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volScalarField> Foam::volWeightedGaussianFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    const tmp<volScalarField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volSymmTensorField> Foam::volWeightedGaussianFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    const tmp<volSymmTensorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volTensorField> Foam::volWeightedGaussianFilter::operator()
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
