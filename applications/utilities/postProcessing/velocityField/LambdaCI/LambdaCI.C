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

Application
    LambdaCI

Description
    Calculates and writes the LambdaCI criterion.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (Uheader.headerOk())
    {
        Info<< "    Reading U" << endl;
        volVectorField U(Uheader, mesh);

        const volTensorField gradU(fvc::grad(U));
        
	volScalarField LambdaCI
        (
            IOobject
            (
                "LambdaCI",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            gradU.component(tensor::XX)
        );

        //Calculating second Invariant of D
        volScalarField b = gradU.component(tensor::XX)*gradU.component(tensor::YY) + gradU.component(tensor::XX)*gradU.component(tensor::ZZ) + gradU.component(tensor::YY)*gradU.component(tensor::ZZ)
                 - gradU.component(tensor::XY)*gradU.component(tensor::YX) - gradU.component(tensor::XZ)*gradU.component(tensor::ZX) - gradU.component(tensor::YZ)*gradU.component(tensor::ZY);
                         
        volScalarField Q = - b/3;
        volScalarField Q3 = pow3(Q);
    
        //Calculating Third Invariant of D
        //c=-Det(D) Third Invariant
        volScalarField c = - gradU.component(tensor::XX)*gradU.component(tensor::YY)*gradU.component(tensor::ZZ) - gradU.component(tensor::XY)*gradU.component(tensor::YZ)*gradU.component(tensor::ZX)
                   - gradU.component(tensor::XZ)*gradU.component(tensor::YX)*gradU.component(tensor::ZY) + gradU.component(tensor::XZ)*gradU.component(tensor::YY)*gradU.component(tensor::ZX)
                   + gradU.component(tensor::XY)*gradU.component(tensor::YX)*gradU.component(tensor::ZZ) + gradU.component(tensor::XX)*gradU.component(tensor::YZ)*gradU.component(tensor::ZY);
                                                              
        volScalarField R = c/2;
        volScalarField R2 = sqr(R);
                
        forAll(LambdaCI, I)
	{                                                               
            //Diskriminante der charakteristischen Gleichung Delta=R2+Q3
            if (R2[I] < Q3[I])
            {
                LambdaCI[I] = 0;
                //Imaginärer Eigenwert ist 0
            }
            else
            {   
                //Berechnung der Reelen Eigenwerte und der Imaginären Anteile der komplexen Eigenwerte.
                scalar R2mQ3 = R2[I]-Q3[I];
	        scalar sqrtR2mQ3 =Foam::sqrt(R2mQ3);
                                                                                                                                 
    	        if (Q[I] < 0)
                {
		    LambdaCI[I] = Foam::pow(sqrtR2mQ3 - R[I], 1.0/3.0) + Foam::pow(sqrtR2mQ3 + R[I], 1.0/3.0);
	        }
    	        else
    	        {
    		    if (R[I] < 0)
		    {
		        LambdaCI[I] = Foam::pow(-R[I] + sqrtR2mQ3, 1.0/3.0) - Foam::pow(-R[I] - sqrtR2mQ3, 1.0/3.0);
		    }
		    else
		    {
		        LambdaCI[I] = -Foam::pow(R[I] - sqrtR2mQ3, 1.0/3.0) + Foam::pow(R[I] + sqrtR2mQ3, 1.0/3.0);
		    }
 	        }
            }
        }

        Info<< "    Writing LambdaCI" << endl;
        LambdaCI.write();
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
