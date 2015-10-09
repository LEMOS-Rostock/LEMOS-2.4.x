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

#include "weightedArithmeticMeanFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(weightedArithmeticMeanFilter, 0);
addToRunTimeSelectionTable(LESfilter, weightedArithmeticMeanFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
weightedArithmeticMeanFilter::weightedArithmeticMeanFilter
(
    const fvMesh& mesh,
    scalar alpha,
    scalar beta
)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    alpha_(alpha),
    beta_(beta)
{
    stencilWeights(alpha_, beta_);
}


weightedArithmeticMeanFilter::weightedArithmeticMeanFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh)),
    alpha_(dict.subDict(type() + "Coeffs").lookupOrDefault<scalar>("alpha", 1.0)),
    beta_(dict.subDict(type() + "Coeffs").lookupOrDefault<scalar>("beta", 1.0))
{
    stencilWeights(alpha_, beta_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void weightedArithmeticMeanFilter::read(const dictionary& dict)
{
    dict.subDict(type() + "Coeffs").readIfPresent<scalar>("alpha", alpha_);
    dict.subDict(type() + "Coeffs").readIfPresent<scalar>("beta", beta_);

    stencilWeights(alpha_, beta_);
}

void weightedArithmeticMeanFilter::stencilWeights(scalar alpha, scalar beta)
{
    const labelListList& untransformedElements = addressing_.untransformedElements();
    const labelListList& transformedElements = addressing_.transformedElements();

    stencilWeights_.resize(mesh().nCells());

    forAll(stencilWeights_, stencilWeightsI)
    {
	label nElements = untransformedElements[stencilWeightsI].size() + transformedElements[stencilWeightsI].size();

        List<scalar>& stencilWeight = stencilWeights_[stencilWeightsI];
	stencilWeight.clear();

	if(nElements > 0)
	{	
	    scalar sumWeight = alpha + beta*(nElements-1);

	    stencilWeight.append(alpha/sumWeight);
	    nElements--;

	    while(nElements)
	    {
	        stencilWeight.append(beta/sumWeight);
	        
                nElements--;
	    }
	}
    }
}

void weightedArithmeticMeanFilter::updateStencil()
{
    // ToDo: recreate stencil after mesh changing
    //addressing_ = (centredCPCCellToCellExtStencilObject::New(mesh()));

    //stencilWeights(alpha_, beta_);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


Foam::tmp<Foam::volVectorField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const 
{
    const tmp<volVectorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volScalarField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    const tmp<volScalarField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volSymmTensorField> Foam::weightedArithmeticMeanFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    const tmp<volSymmTensorField>& filteredField = addressing_.weightedSum(unFilteredField(), stencilWeights_);

    unFilteredField.clear();

    return filteredField;
}

Foam::tmp<Foam::volTensorField> Foam::weightedArithmeticMeanFilter::operator()
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
