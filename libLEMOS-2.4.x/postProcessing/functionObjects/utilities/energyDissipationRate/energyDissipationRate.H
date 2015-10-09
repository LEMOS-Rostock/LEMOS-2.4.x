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
    Foam::energyDissipationRate

Group
    LEMOS UtilitiesFunctionObjects

Description
    This function object calculates effective energy dissipation rate and performs
    averaging in time. The mean field is written into time folder according to
    output time settings. 

    At the moment energy dissipation rate function object can be used for 
    incompressible flows only. The determination of disispation rate for 
    compressible flows has to be implemented yet. 

    Note the code base for compressible flows is already implemented.

    The energy dissipation rate is calculated as follows:

    epsilon = <-tau_ij && S_ij>
  

    Example of function object specification:
    \verbatim
    energyDissipationRate
    {
        type        energyDissipationRate;
        functionObjectLibs ("libLEMOS-2.4.x.so");
        outputControl   outputTime;  // timeStep or outputTime
        outputInterval  10; 

        fieldName       epsilon;     // name of the energy dissipation field to be written
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: energyDissipationRate | yes        |
        outputControl  | timeStep or outputTime     | optional   | 
        outputInterval |                            | optional   |
        fieldName      | name of the energy dissipation field    | yes        |
    \endtable

SourceFiles
    energyDissipationRate.C

\*---------------------------------------------------------------------------*/

#ifndef energyDissipationRate_H
#define energyDissipationRate_H

#include "fieldAverage.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "OFstream.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "compressible/turbulenceModel/turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;
class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class energyDissipationRate Declaration
\*---------------------------------------------------------------------------*/

class energyDissipationRate
{
public:

        static const word turbulenceModelName;

protected:

    // Protected data
        //- 

        //- Name of this set of energyDissipationRate object
        word name_;

        //- Name of this set of energyDissipationRate object
        word fieldName_;
 
         //- Database this class is registered to
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        energyDissipationRate(const energyDissipationRate&);

        //- Disallow default bitwise assignment
        void operator=(const energyDissipationRate&);
 
        //- Return true if compressible turbulence model is identified
        bool compressible();


public:

    //- Runtime type information
    TypeName("energyDissipationRate");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        energyDissipationRate
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~energyDissipationRate();


    // Member Functions

        //- Return name of the set of energyDissipationRate
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the energyDissipationRate data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the energyDissipationRate and write
        virtual void write();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();
  
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
