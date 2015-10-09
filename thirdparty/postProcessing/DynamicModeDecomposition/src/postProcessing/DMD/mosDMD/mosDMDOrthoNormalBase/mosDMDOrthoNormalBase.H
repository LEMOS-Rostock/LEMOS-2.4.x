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


Description
    Class which assembles scalar eigen-base given a set of scalar fields.
    This is performed by establishing a matrix of snapshots, calculating and
    sorting them and providing corresponding eigen-vectors.  Eigen-values are
    sorted in increasing order.  Snapshots are added one at a time and
    when the base is completed, the calculation is triggered


Class
    DMDOrthoNormalBase


Description
    Establish DMD ortho-normal base and interpolation coefficients give a list
    of fields. 


SourceFiles
    DMDOrthoNormalBase.C

\*---------------------------------------------------------------------------*/

#ifndef DMDOrthoNormalBase_H
#define DMDOrthoNormalBase_H

#include "volFields.H"
#include "RectangularMatrix.H"
#include "DMDModel.H"
#include "mosDMDEigenBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class DMDOrthoNormalBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mosDMDOrthoNormalBase
:
    public DMDModel<Type>
{
public:

    // Public typedefs

        typedef GeometricField<Type, fvPatchField, volMesh> GeoTypeField;

        typedef Field<Type> TypeField;



private:

    // Private data

       // Number of modes to compute
       // 0 - all modes
       label baseSize_;

       // Sampling time step, must be uniform
       scalar deltaT_;

       // Sampling time step, must be uniform
       word modeSelection_;

       scalarField energyModes_;
       scalarField modeEnergies_;
       scalarField modeCmptEnergies_;

       //- Ritz values 
       FieldField<Field, complex> ritzValues_;
         
    // Private Member Functions

        //- Disallow default bitwise copy construct
        mosDMDOrthoNormalBase(const mosDMDOrthoNormalBase&);

        //- Disallow default bitwise assignment
        void operator=(const mosDMDOrthoNormalBase&);

        //- Calculates energy contribution of modes
        void calcEnergyModes();

        //- Calculate ortho-normal base
        void calcOrthoBase
        (
            const PtrList<GeoTypeField>& snapshots
        );


public:

    //- Runtime type information
        TypeName("mosDMD");

    // Constructors

        //- Construct given list of snapshots 
        mosDMDOrthoNormalBase
        (
            const fvMesh& mesh,
            const word& DMDModelName = typeName
        );


    // Destructor

        ~mosDMDOrthoNormalBase();


    // Member Functions

        //- Return base size
        label baseSize() 
        {
            return baseSize_;
        }

        
        //- Calculate ortho-normal base
        void calcBase
        (
            PtrList<GeoTypeField>& snapshots
        )
        {
            calcOrthoBase(snapshots);   
        };

        //- Return ritz values
        const FieldField<Field, complex>& ritzValues()
        {
            return ritzValues_;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mosDMDOrthoNormalBase.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
