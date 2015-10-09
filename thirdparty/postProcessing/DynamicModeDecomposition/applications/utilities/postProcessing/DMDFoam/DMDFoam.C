/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University of Rostock
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Application
    DMDFoam

Author
    LEMOS, University of Rostock.  All rights reserved.

Description
    Calculates dynamic mode decomposition of a given field set

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DMDModel.H"
#include "VectorSpace.H"
#include "argList.H"
#include "timeSelector.H"
#include <complex>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void addSnapshotsToPtrList
(
    fvMesh& mesh,
    Time& runTime,
    labelList& timeIndices, 
    const instantList& Times,
    const label& startTime,
    const label& endTime,
    const label& nSnapshots,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& fields, 
    const word& fieldName
)
{
    Info << "Number of snapshots: " << nSnapshots << endl;

    label snapshotI = 0;

    for (label i = startTime; i < endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        Info<< "    Reading " << fieldName << endl;
        fields.set
        (
            snapshotI,
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                ),
                mesh
            )
        );

        // Rename the field
        fields[snapshotI].rename(fieldName + name(i));
        timeIndices[snapshotI] = i;
        snapshotI++;
        Info<< endl;
    }

    timeIndices.setSize(snapshotI);
    runTime.setTime(Times[timeIndices[0]], timeIndices[0]);
}


template<class Type>
void reconstructSnapshots
(
    DMDModel<Type>& DMD,
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots,
    Time& runTime,
    const instantList& Times,
    labelList& timeIndices, 
    const word& fieldName
)
{
    const FieldField<Field, complex>& ritzValues = DMD.ritzValues();

    for(label snapI=0; snapI<snapshots.size(); snapI++)
    {
        GeometricField<Type, fvPatchField, volMesh> reconstruct
        (
            IOobject
            (
                fieldName+"DMDreconstruct",
                runTime.timeName(),
                snapshots[snapI].mesh(),
                IOobject::NO_READ
            ),
            snapshots[snapI].mesh(),
            dimensioned<Type>
            (
                "zero",
                snapshots[snapI].dimensions(),
                pTraits<Type>::zero
            )
        );

        runTime.setTime(Times[timeIndices[0]], timeIndices[0]);


        for (label baseI = 0; baseI < DMD.baseSize(); baseI++)
        {
            GeometricField<Type, fvPatchField, volMesh> DMDmodeRe
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + name(baseI) + "Real",
                    runTime.timeName(), 
                    snapshots[0].mesh(),
                    IOobject::MUST_READ
                ),
                snapshots[0].mesh()
            );

            GeometricField<Type, fvPatchField, volMesh> DMDmodeIm
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + name(baseI) + "Imag",
                    runTime.timeName(),
                    snapshots[0].mesh(),
                    IOobject::MUST_READ
                ),
                snapshots[0].mesh()
            );

            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                volScalarField reconstructCmpt = reconstruct.component(cmpt);

                const complexField& ritzValueCmpt = ritzValues[cmpt];

                std::complex<scalar> tmp;
                tmp.real() = ritzValueCmpt[baseI].Re();
                tmp.imag() = ritzValueCmpt[baseI].Im();

                tmp = std::pow(tmp, snapI);
                    
                reconstructCmpt +=
                    (tmp.real() * DMDmodeRe.component(cmpt) - tmp.imag()*DMDmodeIm.component(cmpt));

                reconstruct.replace(cmpt, reconstructCmpt);
            }
        }

        runTime.setTime(Times[timeIndices[snapI]], timeIndices[snapI]);
        reconstruct.write();

        scalar sumFieldError =
            gSum(mag
            (   
                reconstruct.internalField()
              - snapshots[snapI].internalField()
            ));

        scalar measure = gSum(mag(snapshots[snapI].internalField())) + SMALL;

        Info << "Time step: " << runTime.timeName() << nl 
             << "Field error: absolute = " << sumFieldError << " relative = "
             << sumFieldError/measure << " measure = " << measure
             << endl;

    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, true);
    argList::addOption
    ( 
        "reconstruct",
        "",
        "reconstruct original fields from dynamic modes"
    ); 

    argList args(argc, argv);

    bool reconstruct(args.optionFound("reconstruct"));

#   include "createTime.H"

    instantList Times = timeSelector::select
    (
        runTime.times(),
        args
    );


    // Get times list
    //instantList Times = runTime.times();
    if(Times.empty()) Times = runTime.times();
    

    const label startTime = 0;
    const label endTime = Times.size();
    const label nSnapshots = Times.size();

    labelList timeIndices(nSnapshots);

#   include "createMesh.H"

    IOdictionary DMDProperties
    (
        IOobject
        (
            "DMDProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word fieldName
    (
        DMDProperties.lookup("field")
    );

    runTime.setTime(Times[startTime], startTime);
    // Determine type of field
    IOobject fieldHeader
    (
        fieldName,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (fieldHeader.headerOk())
    {
        const word fieldType = fieldHeader.headerClassName();
        
        
        if (fieldType == volScalarField::typeName)
        {	    
            autoPtr<DMDModel<scalar> > DMDPtr
            (
                DMDModel<scalar>::New(mesh)
            );
            DMDModel<scalar>& DMD(DMDPtr());

            // Create a list of snapshots
            PtrList<volScalarField> fields(nSnapshots);

            addSnapshotsToPtrList(mesh, runTime, timeIndices, Times, startTime, endTime, nSnapshots, fields, fieldName);
            
            // Create DMD base
            DMD.calcBase(fields);

            if(reconstruct)
            { 
                reconstructSnapshots(DMD , fields, runTime, Times, timeIndices, fieldName);
            }
        }
        else if (fieldType == volVectorField::typeName)
        {    
	    // Create vectorDMDModel
            autoPtr<DMDModel<vector> > DMDPtr
            (
                DMDModel<vector>::New(mesh)
            );
            
            DMDModel<vector>& DMD(DMDPtr());
            
            // Create PtrList for snapshots
            PtrList<volVectorField> fields(nSnapshots);
            
            // Add snapshots to PtrList
            addSnapshotsToPtrList(mesh, runTime, timeIndices, Times, startTime, endTime, nSnapshots, fields, fieldName);
            
            // Create DMD base
            DMD.calcBase(fields);
 
            if(reconstruct)
            { 
                reconstructSnapshots(DMD , fields, runTime, Times, timeIndices, fieldName);
            }
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
        }
	    // Create vectorDMDModel
            autoPtr<DMDModel<symmTensor> > DMDPtr
            (
                DMDModel<symmTensor>::New(mesh)
            );
            
            DMDModel<symmTensor>& DMD(DMDPtr());
            
            // Create PtrList for snapshots
            PtrList<volSymmTensorField> fields(nSnapshots);
            
            // Add snapshots to PtrList
            addSnapshotsToPtrList(mesh, runTime, timeIndices, Times, startTime, endTime, nSnapshots, fields, fieldName);
            
            // Create DMD base
            DMD.calcBase(fields);
 
            if(reconstruct)
            { 
                reconstructSnapshots(DMD , fields, runTime, Times, timeIndices, fieldName);
            }
        else if (fieldType == volTensorField::typeName)
        {
	    // Create vectorDMDModel
            autoPtr<DMDModel<tensor> > DMDPtr
            (
                DMDModel<tensor>::New(mesh)
            );
            
            DMDModel<tensor>& DMD(DMDPtr());
            
            // Create PtrList for snapshots
            PtrList<volTensorField> fields(nSnapshots);
            
            // Add snapshots to PtrList
            addSnapshotsToPtrList(mesh, runTime, timeIndices, Times, startTime, endTime, nSnapshots, fields, fieldName);
            
            // Create DMD base
            DMD.calcBase(fields);
 
            if(reconstruct)
            { 
                reconstructSnapshots(DMD , fields, runTime, Times, timeIndices, fieldName);
            }
        }
    }

    Info << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
