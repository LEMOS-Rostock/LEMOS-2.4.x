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

#include "minimumFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(minimumFilter, 0);
addToRunTimeSelectionTable(LESfilter, minimumFilter, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
minimumFilter::minimumFilter
(
    const fvMesh& mesh
)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh))
{}


minimumFilter::minimumFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    addressing_(centredCPCCellToCellExtStencilObject::New(mesh))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void minimumFilter::read(const dictionary& dict)
{}

void minimumFilter::updateStencil()
{
    // ToDo: recreate stencil after mesh changing
    //addressing_ = (centredCPCCellToCellExtStencilObject::New(mesh()));

    //stencilWeights(alpha_, beta_);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::minimumFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const 
{

    tmp<volVectorField> tfilteredField
    (
        new volVectorField
        (
            IOobject
            (
                unFilteredField().name(),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<vector>
            (
                unFilteredField().name(),
                unFilteredField().dimensions(),
                pTraits<vector>::zero
            ),
            zeroGradientFvPatchField<vector>::typeName
        )
    );

    volVectorField& filteredField = tfilteredField();
    
    List<List<vector> > stencilsData(mesh().nCells());
    addressing_.collectData(unFilteredField(), stencilsData);

    forAll(stencilsData, stencilsDataI)
    {
        List<vector>& stencilData = stencilsData[stencilsDataI];
    	
	label sStencilData = stencilData.size();

	vector median = pTraits<vector>::zero;

	for (label i=0; i<3; i++)
    	{
            // Sort stencil elements in ascending order
       	    Foam::sort(stencilData, less<vector>(i));    

            if (sStencilData  > 0)
            {
                median[i] = stencilData[0][i];
            }
    	}

        filteredField[stencilsDataI] = median;
    }

    unFilteredField.clear();

    return tfilteredField;
}

Foam::tmp<Foam::volScalarField> Foam::minimumFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    tmp<volScalarField> tfilteredField
    (
        new volScalarField
        (
            IOobject
            (
                unFilteredField().name(),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<scalar>
            (
                unFilteredField().name(),
                unFilteredField().dimensions(),
                pTraits<scalar>::zero
            ),
            zeroGradientFvPatchField<scalar>::typeName
        )
    );

    volScalarField& filteredField = tfilteredField();

    List<List<scalar> > stencilsData(mesh().nCells());
    addressing_.collectData(unFilteredField(), stencilsData);

    forAll(stencilsData, stencilsDataI)
    {
        List<scalar>& stencilData = stencilsData[stencilsDataI];

        label sStencilData = stencilData.size();

        scalar median = 0.0;

        // Sort stencil elements in ascending order
        Foam::sort(stencilData);
            
        if (sStencilData  > 0)
        {
            median = stencilData[0];
        }
    
        filteredField[stencilsDataI] = median;
    }

    unFilteredField.clear();

    return tfilteredField;

}

Foam::tmp<Foam::volSymmTensorField> Foam::minimumFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    tmp<volSymmTensorField> tfilteredField
    (
        new volSymmTensorField
        (
            IOobject
            (
                unFilteredField().name(),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<symmTensor>
            (
                unFilteredField().name(),
                unFilteredField().dimensions(),
                pTraits<symmTensor>::zero
            ),
            zeroGradientFvPatchField<symmTensor>::typeName
        )
    );

    volSymmTensorField& filteredField = tfilteredField();

    List<List<symmTensor> > stencilsData(mesh().nCells());
    addressing_.collectData(unFilteredField(), stencilsData);

    forAll(stencilsData, stencilsDataI)
    {
        List<symmTensor>& stencilData = stencilsData[stencilsDataI];

        label sStencilData = stencilData.size();

        symmTensor median = pTraits<symmTensor>::zero;

        for (label i=0; i<6; i++)
        {
            // Sort stencil elements in ascending order
            Foam::sort(stencilData, less<symmTensor>(i));
            
            if (sStencilData  > 0)
            {
                median[i] = stencilData[0][i];
            }
        }
            
        filteredField[stencilsDataI] = median;
    }
    
    unFilteredField.clear();

    return tfilteredField;
}

Foam::tmp<Foam::volTensorField> Foam::minimumFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    tmp<volTensorField> tfilteredField
    (
        new volTensorField
        (
            IOobject
            (
                unFilteredField().name(),
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<tensor>
            (
                unFilteredField().name(),
                unFilteredField().dimensions(),
                pTraits<tensor>::zero
            ),
            zeroGradientFvPatchField<tensor>::typeName
        )
    );

    volTensorField& filteredField = tfilteredField();

    List<List<tensor> > stencilsData(mesh().nCells());
    addressing_.collectData(unFilteredField(), stencilsData);

    forAll(stencilsData, stencilsDataI)
    {
        List<tensor>& stencilData = stencilsData[stencilsDataI];

        label sStencilData = stencilData.size();

        tensor median = pTraits<tensor>::zero;

        for (label i=0; i<9; i++)
        {
            // Sort stencil elements in ascending order
            Foam::sort(stencilData, less<tensor>(i));
            
            if (sStencilData  > 0)
            {
                median[i] = stencilData[0][i];
            }
        }
        filteredField[stencilsDataI] = median;
    }

    unFilteredField.clear();

    return tfilteredField;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
