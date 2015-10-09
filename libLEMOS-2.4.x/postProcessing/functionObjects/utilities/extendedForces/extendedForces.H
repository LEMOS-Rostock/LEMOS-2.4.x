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
    Foam::extendedForces

Group
    LEMOS UtilitiesFunctionObjects

Description
    This function object calculates additionally viscous and pressure forces
    for wall boundary patches.

    Reference:  Insight CAE, a workbench for Computer-Aided Engineering 
                Hannes Kroeger <hannes@kroegeronline.net>

    Example of function object specification:
    \verbatim
    extendedForces
    {
        type        extendedForces;
        functionObjectLibs ("libLEMOS-2.4.x.so");
        ...
    }
    \endverbatim

    \heading Function object usage

    See OpenFOAM forces function object

SourceFiles
    extendedForces.C

\*---------------------------------------------------------------------------*/

#ifndef extendedForces_H
#define extendedForces_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class extendedForces Declaration
\*---------------------------------------------------------------------------*/

class extendedForces
: public forces
{

protected:

    // Protected data
    
        //- Additional force fields
        autoPtr<volVectorField> pressureForce_, viscousForce_;

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        extendedForces(const extendedForces&);

        //- Disallow default bitwise assignment
        void operator=(const extendedForces&);
 

public:

    //- Runtime type information
    TypeName("extendedForces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        extendedForces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

        //- Construct from components
        extendedForces
        (
            const word& name,
            const objectRegistry&,
            const labelHashSet& patchSet,
            const word& pName,
            const word& UName,
            const word& rhoName,
            const scalar rhoInf,
            const scalar pRef,
            const coordinateSystem& coordSys
        );


    //- Destructor
    virtual ~extendedForces();


    // Member Functions

        //- Execute, currently does nothing
        virtual void execute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
