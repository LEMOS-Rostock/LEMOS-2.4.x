Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading field f\n" << endl;
volScalarField f
(
    IOobject
    (
        "f",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField scalarDiss
(
    IOobject
    (
        "scalarDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField scalarVar
(
    IOobject
    (
        "scalarVar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    0.0
);


Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//scalar D = transportProperties.lookupOrDefault<scalar>("D_" +  f.name(), 0.0);
//Info<< "Reading molecular diffusivity coefficient: "<< D << endl;

scalar Sc = transportProperties.lookupOrDefault<scalar>("Sc" +  f.name(), 1000);
scalar Sct = transportProperties.lookupOrDefault<scalar>("Sc_t" +  f.name(), 0.9);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
