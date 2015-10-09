/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University Rostock
     \\/     M anipulation  |
------------------------------------------------------------------------------
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

Application
    postChannelExt

Description
    Post-processes data from channel flow calculations.

    Assuming that the mesh is periodic in the x and z directions, 
    collapse given fields to a line and print them as graph output

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "channelIndex.H"
#include "makeGraph.H"
#include "Triple.H"
#include "volFieldsFwd.H"

#include "OSspecific.H"




void processField
(
    const volScalarField& field,
    const channelIndex& channelIndexing,
    const word& identifierName,
    const fileName& path,
    const word& gFormat,
    const word& moment,
    const bool asymmetric=false

)
{
    scalarField fieldValues
    (
        channelIndexing.collapse(field, asymmetric)
    );

    if(moment == "rms")
    {
        fieldValues = sqrt(mag(fieldValues));
    }

    const scalarField& y = channelIndexing.y();
    
    Info << "    Creating graph " << identifierName << endl;
    makeGraph(y, fieldValues, identifierName, path, gFormat);
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"
#   include "readTransportProperties.H"

    const word& gFormat = runTime.graphFormat();

    // Setup channel indexing for averaging over channel down to a line

    IOdictionary channelDict
    (
        IOobject
        (
            "postChannelExtDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    channelIndex channelIndexing(mesh, channelDict);

    List<Triple<word> > fields(channelDict.lookup("fieldsAndIdentifiers"));

    bool backwardsCompatibility(channelDict.lookupOrDefault("backwardsCompatibility", true));

    // For each time step read all fields
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "Processing fields for time " << runTime.timeName() << endl;
        Info << endl;

        fileName path(runTime.rootPath()/runTime.caseName()/"graphs"/runTime.timeName());
        mkDir(path);

        forAll(fields, I)
        {
           const word& fieldName = fields[I].first();
           const word& identifierName = fields[I].second();
           const word& moment = fields[I].third();

           IOobject fieldHeader
           (
               fieldName,
               runTime.timeName(),
               mesh,
               IOobject::MUST_READ
           );

           if (!fieldHeader.headerOk())
           {
               Info << endl;
               Info<< "No " << fieldName <<" field" << endl;
               Info << endl;
               continue;
           }

           Info << endl;
           Info << fieldName << endl;

           if(fieldHeader.headerClassName() == volScalarField::typeName)
           {
               Info << endl;
               Info << "    Reading field" << endl;
              
               volScalarField field(fieldHeader, mesh);

               Info << "    Collapsing field"  << endl;
               Info << endl;
               processField(field, channelIndexing, identifierName, path, gFormat, moment);
           }

           if(fieldHeader.headerClassName() == volVectorField::typeName)
           {
               Info << endl;
               Info << "    Reading field" << endl;
               
               volVectorField field(fieldHeader, mesh);
               
               Info << "    Collapsing field"  << endl;
               Info << endl;
               
               
               if((fieldHeader.name() == "U" || fieldHeader.name() == "UMean") && backwardsCompatibility)
               {
                   processField(field.component(vector::X)(), channelIndexing, identifierName, path, gFormat, "mean");
               }
               else
               {
                   processField(field.component(vector::X)(), channelIndexing, identifierName+"_x", path, gFormat, moment);
                   processField(field.component(vector::Y)(), channelIndexing, identifierName+"_y", path, gFormat, moment);
                   processField(field.component(vector::Z)(), channelIndexing, identifierName+"_z", path, gFormat, moment);
                    
               }
           }

           if(fieldHeader.headerClassName() == volSymmTensorField::typeName)
           {
               Info << endl;
               Info << "    Reading field" << endl;

               volSymmTensorField field(fieldHeader, mesh);

               Info << "    Collapsing field"  << endl;
               Info << endl;
               
               if((fieldHeader.name() == "UPrime2Mean") && backwardsCompatibility)
               {
                   processField(field.component(symmTensor::XX)(), channelIndexing, "u", path, gFormat, "rms");
                   processField(field.component(symmTensor::YY)(), channelIndexing, "v", path, gFormat, "rms");
                   processField(field.component(symmTensor::ZZ)(), channelIndexing, "w", path, gFormat, "rms");
                   processField(field.component(symmTensor::XY)(), channelIndexing, "uv", path, gFormat, "mean", true);
               }
               else
               {
                   processField(field.component(symmTensor::XX)(), channelIndexing, identifierName+"_xx", path, gFormat, moment);
                   processField(field.component(symmTensor::XY)(), channelIndexing, identifierName+"_xy", path, gFormat, moment, true);
                   processField(field.component(symmTensor::XZ)(), channelIndexing, identifierName+"_xz", path, gFormat, moment, true);
                   processField(field.component(symmTensor::YY)(), channelIndexing, identifierName+"_yy", path, gFormat, moment);
                   processField(field.component(symmTensor::YZ)(), channelIndexing, identifierName+"_yz", path, gFormat, moment, true);
                   processField(field.component(symmTensor::ZZ)(), channelIndexing, identifierName+"_zz", path, gFormat, moment);
               }
           }

           if(fieldHeader.headerClassName() == volTensorField::typeName)
           {
               Info << endl;
               Info << "    Reading field" << endl;

               volTensorField field(fieldHeader, mesh);

               Info << "    Collapsing field"  << endl;
               Info << endl;

               processField(field.component(tensor::XX)(), channelIndexing, identifierName+"_xx", path, gFormat, moment);
               processField(field.component(tensor::XY)(), channelIndexing, identifierName+"_xy", path, gFormat, moment, true);
               processField(field.component(tensor::XZ)(), channelIndexing, identifierName+"_xz", path, gFormat, moment, true);
               processField(field.component(tensor::YX)(), channelIndexing, identifierName+"_yx", path, gFormat, moment, true);
               processField(field.component(tensor::YY)(), channelIndexing, identifierName+"_yy", path, gFormat, moment);
               processField(field.component(tensor::YZ)(), channelIndexing, identifierName+"_yz", path, gFormat, moment, true);
               processField(field.component(tensor::XZ)(), channelIndexing, identifierName+"_xz", path, gFormat, moment, true);
               processField(field.component(tensor::YZ)(), channelIndexing, identifierName+"_yz", path, gFormat, moment, true);
               processField(field.component(tensor::ZZ)(), channelIndexing, identifierName+"_zz", path, gFormat, moment);
           }


        }

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
