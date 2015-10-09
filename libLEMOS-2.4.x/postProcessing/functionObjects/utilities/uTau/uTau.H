/*---------------------------------------------------------------------------*\ 
| File modified by LEMOS (University of Rostock) 2013                         |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::uTau

Group
    LEMOS UtilitiesFunctionObjects

Description
    This function object evaluates and outputs the friction velocity at wall
    patches.  The result is written into the file 
    ..\postprocessing\<FunctionObjectName>\<timeDir>\uTau.dat
    every time step or at userdefined timeStep/outputTime. 

    The friction velocity is calculated using two methods:

    1.) sqrt(|n*Reff|) 
    2.) sqrt(nu*|snGrad(U)|)

    The function object works for compressible and incompressible flows,
    provided that the case setup relies on a turbulence model, which provides
    the "laminar" or a "turbulent" turbulence model.

    All wall patches are included by default (global); to restrict the calculation to 
    certain wall patches, use the optional 'patches' entry. Different log levels 
    provide detailed information about min/max/avg friction velocity values. 
    Default setting is no messages during solver run.

    Example of function object specification:
    \verbatim
    uTau1
    {
        type        uTau;
        functionObjectLibs ("libLEMOS-2.4.x.so");
        log 	    0; // 0 - no output (default)
                       // 1 - averaged friction velocity only 
                       // 2 - min/max/avg friction velocity (all wall patches)
                       // 3 - min/max/avg friction velocity (single wall patches)
        patches     (".*Wall"); // optional
        outputControl   timeStep;  // timeStep or outputTime
        outputInterval  10; 
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: uTau | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

SourceFiles
    uTau.C
    IOuTau.H

\*---------------------------------------------------------------------------*/

#ifndef uTau_H
#define uTau_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "OFstream.H"

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
                       Class uTau Declaration
\*---------------------------------------------------------------------------*/

class uTau
:
    public functionObjectFile
{

protected:

    // Protected data

        //- Name of this set of uTau object
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        label log_;

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Total time
        scalar totalTime_;

        //- Total time
        label totalIter_;

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the friction velocity using wall normal
        //  and stress tensor Reff from turbulence model (traction)
        void calcFrictionVelocity
        (
            const fvMesh& mesh,
            const volSymmTensorField& Reff
        );

        //- Calculate the friction velocity using wall normal
        //  and stress tensor Reff from turbulence model
        void calcFrictionVelocity2
        (
            const fvMesh& mesh,
            const volSymmTensorField& Reff
        );

        //- Calculate the friction velocity using wall normal 
        //  and resolved velocity U at the wall
        void calcFrictionVelocity
        (
            const fvMesh& mesh,
            const volVectorField& U,
            const volScalarField& nu
        );

        //- Disallow default bitwise copy construct
        uTau(const uTau&);

        //- Disallow default bitwise assignment
        void operator=(const uTau&);


public:

    //- Runtime type information
    TypeName("uTau");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        uTau
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~uTau();


    // Member Functions

        //- Return name of the set of uTau
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the uTau data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the uTau and write
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
