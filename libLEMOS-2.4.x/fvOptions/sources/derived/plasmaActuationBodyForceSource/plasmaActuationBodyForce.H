/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::fv::volumeDrag

Description
    Creates a pressure gradient source

    Note: Currently only handles kinematic pressure

    \heading Source usage
    Example usage:
    \verbatim
    volumeDragCoeffs
    {
        fieldNames  (U);                        // name of velocity field
        Ubar        (10.0 0 0);                 // desired average velocity
        gradPini    gradPini [0 2 -2 0 0] 0;    // initial pressure gradient
        flowDir     (1 0 0);                    // flow direction
    }
    \endverbatim


SourceFiles
    volumeDrag.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaActuationBodyForce_H
#define plasmaActuationBodyForce_H

//#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Friend declaration
class bodyForceModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class volumeDrag Declaration
\*---------------------------------------------------------------------------*/

class plasmaActuationBodyForce
:
    public option
{
    // Private data

    autoPtr<bodyForceModel> bodyForce_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        plasmaActuationBodyForce(const plasmaActuationBodyForce&);

        //- Disallow default bitwise assignment
        void operator=(const plasmaActuationBodyForce&);


public:

    //- Runtime type information
    TypeName("plasmaActuationBodyForce");


    // Constructors

        //- Construct from explicit source name and mesh
        plasmaActuationBodyForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
