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

Class
    Foam::oscillatingNormalVelocityFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an oscillating velocity condition in terms of
    amplitude and frequency. The boundary condition can be used to simulate 
    suction and blowing.

        \f[ 
            x_p = a1 f(x) sin(2 \pi f t) + a2 f(x) g(z) sin(\pi f t + p)
        \f]

    where

    \vartable
        x_p     | patch values
        f(x)    | modulation function in x direction
        g(z)    | modulation function in z direction
        p       | phase
        a1,a2   | amplitudes
        f       | frequency [1/s]
        t       | time [s]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        modulation   | functions values        | yes         |
        functions    |                         | yes         |
        phase        | phase value             | yes         | 
        amplitudes   | oscillation amplitudes  | yes         |
        frequency    | oscillation frequency   | yes         |
        fluxCorr     | flux adjustment          | no          | true
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            oscillatingNormalVelocity;
        refValue        uniform 5.0;
        phase           0.0;
        amplitude1      constant 1;
        amplitude2      constant 2;
        frequency       constant 10;
        funcX		csvFile; // Optional: OpenFOAM tableFile
        funcXCoeffs
        {
            nHeaderLine     0;
            refColumn       0;
            componentColumns 1(1);
            separator       " ";
            mergeSeparators 1;
            fileName        "$FOAM_CASE/constant/funcX";
        }
        funcZ           csvFile; // Optional: OpenFOAM tableFile
        funcZCoeffs
        {
            nHeaderLine     0;
            refColumn       0;
            componentColumns 1(1);
            separator       " ";
            mergeSeparators 1;
            fileName        "$FOAM_CASE/constant/funcZ";
        }
        fluxCorr        true;

    }
    \endverbatim

Note
    The amplitudes and frequency entries are DataEntry types, able to describe
    time varying functions.  The example above gives the usage for supplying
    constant values. The functions f(x) and g(z) are DataEntry supplied as CSV 
    or tableFile, first column is a parameter between 0 and 1, second column is
    function value. The provided functions must be scaled to 0 and 1, the function
    values should be bounded by 1 and -1. The boundary conditon transforms the
    spatial coordinates automatically to the interval [0,1].

    The boundary condition assumes the normal component in y direction.

SeeAlso
    Foam::DataEntry

SourceFiles
    oscillatingNormalVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingNormalVelocityFvPatchVectorField_H
#define oscillatingNormalVelocityFvPatchVectorField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "DataEntry.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class oscillatingNormalVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class oscillatingNormalVelocityFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
    // Private data

        //- Offset
        vectorField fixedValue_;

        //- Amplitude1
        autoPtr<DataEntry<scalar> > amplitude1_;

        //- Amplitude2
        autoPtr<DataEntry<scalar> > amplitude2_;

        //- Frequency
        autoPtr<DataEntry<scalar> > frequency_;

        //- Phase
        autoPtr<DataEntry<scalar> > phase_;

        //- Function1
        autoPtr<DataEntry<scalar> > funcX_;

        //- Function2
        autoPtr<DataEntry<scalar> > funcZ_;

        //- Function values of function1
        scalarField funcXValues_;

        //- Function values of function2
        scalarField funcZValues_;

        //- Current time index
        label curTimeIndex_;

        //- Correct patch flux 
        bool fluxCorr_;

    // Private Member Functions

        //- Return current patch dimensions
        void calcFunctionValues();

        //- Return current scale
        void currentScale();


public:

    //- Runtime type information
    TypeName("oscillatingNormalVelocity");


    // Constructors

        //- Construct from patch and internal field
        oscillatingNormalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillatingNormalVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oscillatingNormalVelocityFvPatchVectorField
        //  onto a new patch
        oscillatingNormalVelocityFvPatchVectorField
        (
            const oscillatingNormalVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oscillatingNormalVelocityFvPatchVectorField
        (
            const oscillatingNormalVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new oscillatingNormalVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oscillatingNormalVelocityFvPatchVectorField
        (
            const oscillatingNormalVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new oscillatingNormalVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
