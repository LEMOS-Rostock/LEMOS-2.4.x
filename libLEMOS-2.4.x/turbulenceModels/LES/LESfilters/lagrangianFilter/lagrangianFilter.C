
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

\*---------------------------------------------------------------------------*/

#include "lagrangianFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lagrangianFilter, 0);
    addToRunTimeSelectionTable(LESfilter, lagrangianFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianFilter::lagrangianFilter(const fvMesh& mesh, scalar widthCoeff, scalar theta)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    theta_(theta),
    coeff_
    (
        IOobject
        (
            "lagrangianFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    ),
    thetaCoeff_
    (
        IOobject
        (
            "lagrangianFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimTime, theta_),
        calculatedFvPatchScalarField::typeName
    )
{
    coeff_.dimensionedInternalField() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;
}


Foam::lagrangianFilter::lagrangianFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilter(mesh),
    widthCoeff_(readScalar(bd.subDict(type() + "Coeffs").lookup("widthCoeff"))),
    theta_(readScalar(bd.subDict(type() + "Coeffs").lookup("theta"))),
    coeff_
    (
        IOobject
        (
            "lagrangianFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    ),
    thetaCoeff_
    (
        IOobject
        (
            "lagrangianFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimTime, theta_),
        calculatedFvPatchScalarField::typeName
    )
{
    coeff_.dimensionedInternalField() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lagrangianFilter::read(const dictionary& bd)
{
    bd.subDict(type() + "Coeffs").lookup("widthCoeff") >> widthCoeff_;
    bd.subDict(type() + "Coeffs").lookup("theta") >> theta_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::lagrangianFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    const surfaceScalarField& phi =
    mesh().
    lookupObject<surfaceScalarField>
    (
    "phi"
    );

    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volScalarField> filteredField =
        unFilteredField() + deltaT*fvc::ddt(unFilteredField()) + deltaT*fvc::div(phi, unFilteredField()) - fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::lagrangianFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    const surfaceScalarField& phi =
    mesh().
    lookupObject<surfaceScalarField>
    (
    "phi"
    );
    
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volVectorField> filteredField =
        //unFilteredField() + deltaT*fvc::ddt(unFilteredField()) + deltaT*fvc::div(phi, unFilteredField()) - fvc::laplacian(coeff_, unFilteredField());
        unFilteredField() + thetaCoeff_*fvc::ddt(unFilteredField()) + thetaCoeff_*fvc::div(phi, unFilteredField()) - fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::lagrangianFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    const surfaceScalarField& phi =
    mesh().
    lookupObject<surfaceScalarField>
    (
    "phi"
    );

    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volSymmTensorField> filteredField =
        unFilteredField() + deltaT*fvc::ddt(unFilteredField()) + deltaT*fvc::div(phi, unFilteredField()) - fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::lagrangianFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{ 
    const surfaceScalarField& phi =
    mesh().
    lookupObject<surfaceScalarField>
    (
    "phi"
    );

    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volTensorField> filteredField =
        unFilteredField() + deltaT*fvc::ddt(unFilteredField()) + deltaT*fvc::div(phi, unFilteredField()) - fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //
