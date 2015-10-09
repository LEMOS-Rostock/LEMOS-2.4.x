/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University of Rostock
     \\/     M anipulation  |
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

#include "oscillatingNormalVelocityFvPatchVectorField.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void oscillatingNormalVelocityFvPatchVectorField::calcFunctionValues() 
{
    // get patch dimensions
    const vectorField& patchPoints = this->patch().patch().localPoints();

    vector minPoint = gMin(patchPoints);
    vector maxPoint = gMax(patchPoints);

    // transform coordinates into parameter space (0 < t < 1)
    // and calculate function values, must be done after every
    // mesh change otherwise only one time
    const vectorField& faceCentres = this->patch().Cf(); 
    
    forAll(faceCentres, fcI)
    {
        const vector& fc = faceCentres[fcI];

        vector tfc;
        tfc.x() =  (minPoint.x() - fc.x())/ (minPoint.x() - maxPoint.x());
        tfc.y() =  -1;
        tfc.z() =  (minPoint.z() - fc.z())/ (minPoint.z() - maxPoint.z());
        
        funcXValues_[fcI] = funcX_->value(tfc.x());
        funcZValues_[fcI] = funcZ_->value(tfc.z());
    }
}

void oscillatingNormalVelocityFvPatchVectorField::currentScale() 
{
    const scalar t = this->db().time().timeOutputValue();
    const scalar a1 = amplitude1_->value(t);
    const scalar a2 = amplitude2_->value(t);
    const scalar f = frequency_->value(t);
    const scalar p = phase_->value(t);

    scalarField yValues = a1*funcXValues_*sin(constant::mathematical::twoPi*f*t) + a2*funcXValues_*funcZValues_*sin(constant::mathematical::pi*f*t + p);

    vectorField uncorrFixedValue(this->size(), pTraits<vector>::zero);
    uncorrFixedValue.replace(vector::Y, yValues);

    const vectorField nHat(this->patch().nf());
    fixedValue_ = nHat*(nHat & uncorrFixedValue);

    if(fluxCorr_)
    {
        fixedValue_ = fixedValue_ - this->patch().Sf()* gSum(fixedValue_&this->patch().Sf())/gSum(this->patch().magSf());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingNormalVelocityFvPatchVectorField::oscillatingNormalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    fixedValue_(p.size(), pTraits<vector>::zero),
    amplitude1_(),
    amplitude2_(),
    frequency_(),
    phase_(),
    funcX_(),
    funcZ_(),
    funcXValues_(p.size(), 0.0),
    funcZValues_(p.size(), 0.0),
    curTimeIndex_(-1),
    fluxCorr_(true)
{}


oscillatingNormalVelocityFvPatchVectorField::oscillatingNormalVelocityFvPatchVectorField
(
    const oscillatingNormalVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    fixedValue_(ptf.fixedValue_),
    amplitude1_(ptf.amplitude1_().clone().ptr()),
    amplitude2_(ptf.amplitude2_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    phase_(ptf.phase_().clone().ptr()),
    funcX_(ptf.funcX_().clone().ptr()),
    funcZ_(ptf.funcZ_().clone().ptr()),
    funcXValues_(ptf.funcXValues_),
    funcZValues_(ptf.funcZValues_),
    curTimeIndex_(-1),
    fluxCorr_(true)
{}


oscillatingNormalVelocityFvPatchVectorField::oscillatingNormalVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    fixedValue_(p.size(), pTraits<vector>::zero),
    amplitude1_(DataEntry<scalar>::New("amplitude1", dict)),
    amplitude2_(DataEntry<scalar>::New("amplitude2", dict)),
    frequency_(DataEntry<scalar>::New("frequency", dict)),
    phase_(DataEntry<scalar>::New("phase", dict)),
    funcX_(DataEntry<scalar>::New("funcX", dict)),
    funcZ_(DataEntry<scalar>::New("funcZ", dict)),
    funcXValues_(p.size(), 0.0),
    funcZValues_(p.size(), 0.0),
    curTimeIndex_(-1),
    fluxCorr_(dict.lookupOrDefault("fluxCorr", true))
{
    calcFunctionValues();

    if (dict.found("value"))
    {
        fixedValueFvPatchField<vector>::operator==
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        currentScale();
        
        fixedValueFvPatchField<vector>::operator==
        (
            fixedValue_
        );
    }
}


oscillatingNormalVelocityFvPatchVectorField::oscillatingNormalVelocityFvPatchVectorField
(
    const oscillatingNormalVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    fixedValue_(ptf.fixedValue_),
    amplitude1_(ptf.amplitude1_().clone().ptr()),
    amplitude2_(ptf.amplitude2_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    phase_(ptf.phase_().clone().ptr()),
    funcX_(ptf.funcX_().clone().ptr()),
    funcZ_(ptf.funcZ_().clone().ptr()),
    funcXValues_(ptf.funcXValues_),
    funcZValues_(ptf.funcZValues_),
    curTimeIndex_(-1),
    fluxCorr_(true)
{}


oscillatingNormalVelocityFvPatchVectorField::oscillatingNormalVelocityFvPatchVectorField
(
    const oscillatingNormalVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    fixedValue_(ptf.fixedValue_),
    amplitude1_(ptf.amplitude1_().clone().ptr()),
    amplitude2_(ptf.amplitude2_().clone().ptr()),
    frequency_(ptf.frequency_().clone().ptr()),
    phase_(ptf.phase_().clone().ptr()),
    funcX_(ptf.funcX_().clone().ptr()),
    funcZ_(ptf.funcZ_().clone().ptr()),
    funcXValues_(ptf.funcXValues_),
    funcZValues_(ptf.funcZValues_),
    curTimeIndex_(-1),
    fluxCorr_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oscillatingNormalVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
}


void oscillatingNormalVelocityFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);
}


void oscillatingNormalVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        currentScale();

        fixedValueFvPatchField<vector>::operator==
        (
            fixedValue_         
        );

        curTimeIndex_ = this->db().time().timeIndex();
    }


    fixedValueFvPatchField<vector>::updateCoeffs();
}


void oscillatingNormalVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchField<vector>::write(os);
    amplitude1_->writeData(os);
    amplitude2_->writeData(os);
    frequency_->writeData(os);
    phase_->writeData(os);
    funcX_->writeData(os);
    funcZ_->writeData(os);
    os.writeKeyword("fluxCorr") << fluxCorr_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    oscillatingNormalVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
