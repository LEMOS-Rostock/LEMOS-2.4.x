#include "MFM.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "volFields.H"
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
 {
namespace incompressible
 {
namespace LESModels
 {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(MFM, 0);
addToRunTimeSelectionTable(LESModel, MFM, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

MFM::MFM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
  LESModel(modelName, U, phi, transport, turbulenceModelName),
  thermoDict
       (
           IOobject
           (
               "thermophysicalProperties",
               U.mesh().time().constant(),
               U.mesh(),
               IOobject::MUST_READ,
               IOobject::NO_WRITE
           )
       ),
  smoothingFilterPtr(LESfilter::New(U.mesh(), coeffDict(), "smoothingFilter")),
  smoothingFilter(smoothingFilterPtr()),
  testFilterPtr(LESfilter::New(U.mesh(), coeffDict(), "testFilter")),
  testFilter(testFilterPtr()),
  viscLengthScale_
  (
     IOobject
     (
      "viscLengthScale",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      0.0,
      zeroGradientFvPatchScalarField::typeName
  ),
  Ureynolds_
  (
     IOobject
     (
      "UsgsUsgs",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedSymmTensor("UsgsUsgs", dimensionSet(0, 2, -2, 0, 0, 0, 0), pTraits<symmTensor>::zero)
  ),
  uDelta_
  (
     IOobject
     (
      "uDelta",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     U-testFilter(U)
   ),
  N_
  (
     IOobject
     (
      "N",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("N", dimless, 10.0),
      zeroGradientFvPatchScalarField::typeName
  ),
 Re_
  (
     IOobject
     (
      "Re",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("Re", dimless, 10.0),
      zeroGradientFvPatchScalarField::typeName
  ),
  Cai_
  (
     IOobject
     (
      "Cai",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
      mesh_,
      dimensionedScalar("Cai_", dimless, 0.0),
      zeroGradientFvPatchScalarField::typeName
   )
{
  updateSubGridScaleFields(dev(symm(fvc::grad(U)))); 

  printCoeffs();
  
  Info << "MFM created" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void MFM::updateSubGridScaleFields(const volSymmTensorField& S)
{
   updateCai(S);
   updateN(S);
}


// Reynolds number based on the norm of strain rate tensor
void MFM::ReS(const volSymmTensorField& S) 
{
   Re_ = mag(S)*sqr(delta())/nu();
   max(Re_,Re_,1.0);
   Info << "updating Reynolds number ReS: " << average(Re_) << endl;

   Re_.correctBoundaryConditions();
}


// Reynolds number based on the local velocity
void MFM::ReU(const volSymmTensorField& S)
{
   Re_ = max( mag(U())*(delta())/nu() , 1.0);
//   Re_ = max(sqrt(Re_* magSqr(fvc::curl(U()))*sqr(sqr(delta())/nu()))  , 1.0);
   Info << "updating Reynolds number ReU: " << average(Re_) << endl;
   Re_.correctBoundaryConditions();
}
       

// Anisotropy factor, based on Reynolds number 
void MFM::updateCai(const volSymmTensorField& S) 
{
  volSymmTensorField ss = dev(symm(fvc::grad(uDelta_)));
  
  //volScalarField Re = (max( mag(ss)*sqr(delta())/nu(), 1.0));
  volScalarField Re = (max( mag(S)*sqr(delta())/nu(), 1.0));
  //Re = Re/max(Re);

  ////Cai_= 1.0-pow(Re,-3.0/16.0);
  //Cai_ = 2.0/constant::mathematical::pi * atan(Bsgs()*(Re-1));
  //
  Cai_ = 1.0-(1.0/(1+0.00204*pow(Re, 3))+(0.00143*pow(Re, 3))/(1+0.0005863*pow(Re, 4)));
 
   Info << "updating anisotropy factor: " << average(Cai_) << endl;
   Cai_.correctBoundaryConditions();
}


scalar MFM::Bsgs() const
{
   return readScalar(thermoDict.lookup("Bsgs"));
}

tmp<volSymmTensorField> MFM::B() const
{
  return  dev(symm(-(F1()*U()*uDelta_+ F1()*uDelta_*U()) - (sqr(F1())*uDelta_*uDelta_))) ;
}


tmp<volScalarField> MFM::F1() const
{
  return  Bsgs()*Cai_*pow((1-pow(alpha,-4.0/3.0)),-0.5)*pow(2.0,-2.0/3.0*N_)*sqrt(pow(2.0,4.0/3.0*N_)-1.0);
 //return  Bsgs()*pow((1-pow(alpha,-4.0/3.0)),-0.5)*pow(2.0,-2.0/3.0*N_)*sqrt(pow(2.0,4.0/3.0*N_)-1.0);
}


tmp<volSymmTensorField> MFM::devReff() const
{
  return  - nu()*dev(twoSymm(fvc::grad(U()))) - dev(symm((F1()*U()*uDelta_+ F1()*uDelta_*U() + sqr(F1())*uDelta_*uDelta_)));
}



tmp<volScalarField> MFM::epsilon() const
{
  return
  (
     -dev(symm(F1()*U()*uDelta_+ F1()*uDelta_*U() + sqr(F1())*uDelta_*uDelta_))
     &&
     dev(symm(fvc::grad(U())))
  );
}


tmp<volScalarField> MFM::k() const
{
   return tmp<volScalarField>
   (
      new volScalarField
      (
         IOobject
         (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         mesh_,
         dimensionedScalar("k", nu()().dimensions(), 0.0)
      )
   );
}



tmp<fvVectorMatrix> MFM::divDevReff(volVectorField& U) const
{
  surfaceScalarField phiDelta = fvc::interpolate(F1()*uDelta_)&U.mesh().Sf();

  return
  (
     - fvm::laplacian(nu(), U, "laplacian(nu,U)")
     + fvm::div(phiDelta, U, "div(MFM1)")
     - fvc::div(nu()*dev(T(fvc::grad(U))))
     + fvc::div((F1()*U*uDelta_ + sqr(F1()*uDelta_)),  "div(MFM2)")
     //+ fvc::div(smoothingFilter(F1()*uDelta_*U + F1()*U*uDelta_ + sqr(F1()*uDelta_)),  "div(MFM2)")
     //+ fvc::div((F1()*uDelta_*U + F1()*U*uDelta_ + sqr(F1()*uDelta_)),  "div(MFM2)")
  );
}

tmp<fvVectorMatrix> MFM::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

tmp<volScalarField> MFM::D() const
{
  return  Dsgs()*pow((1-pow(alpha,-4.0/3.0)),-0.5)*pow(2.0,-2.0/3.0*N_)*sqrt(pow(2.0,4.0/3.0*N_)-1.0);
}

scalar MFM::Dsgs() const
{
   return gAverage(Cai_)*Bsgs();
}


tmp<volScalarField> MFM::molecularDiffusivityCoeff(word name) const
{
  return
    (
     laminarDiffusivity_[name] 
     );
}


tmp<volVectorField> MFM::Feff(const volScalarField &f) const
{
  return
    (
      - ( laminarDiffusivity_[f.name()] ) * fvc::grad(f)
     );
}


tmp<fvScalarMatrix> MFM::divFeff(volScalarField &f) const
{
    surfaceScalarField phiDelta = fvc::interpolate(F1()*uDelta_)&U().mesh().Sf();

  return
     (
      - fvm::laplacian(laminarDiffusivity_[f.name()], f, "laplacian(Deff,F)")
      + fvm::div(phiDelta, f, "div(MFM1)")
      + fvc::div(D()*U() * (*deltaScalarFields_[f.name()]) + F1()*D()*uDelta_*(*deltaScalarFields_[f.name()]))
     );
}


void MFM::correct(const tmp<volTensorField>& gradU)
{
  LESModel::correct(gradU);

  uDelta_ = U() - testFilter(U());
   
  for ( HashTable<volScalarField&,word>::iterator iter=registeredScalarFields_.begin(); iter!=registeredScalarFields_.end(); iter++ )
  {
      const volScalarField& registeredScalarField = iter();

      (*deltaScalarFields_[registeredScalarField.name()]) = registeredScalarField-testFilter(registeredScalarField);
  } 

  viscLengthScale_ = F1();
  Ureynolds_ = B();

  updateSubGridScaleFields(dev(symm(gradU)));
}


void MFM::updateN(const volSymmTensorField& S)
{
    Info << "calculating number of cascade steps" << endl;

    ReU(S);
    // ReS(S);

    scalar cv = readScalar(thermoDict.lookup("cv"));

    scalarField tmp = cv*pow(Re_,0.75);
    max(tmp,tmp,1.0);

    N_.internalField() = log(tmp)/log(2.0);
    N_.correctBoundaryConditions();


    for ( HashTable<volScalarField&,word>::iterator iter=registeredScalarFields_.begin(); iter!=registeredScalarFields_.end(); iter++ )
    {
        const volScalarField& registeredScalarField = iter();
      
        (*NScalarFields_[registeredScalarField.name()]).internalField() = log(tmp * pow(nu()/laminarDiffusivity_[registeredScalarField.name()], 0.5))/log(2.0) ;
    }

    Info << "updating number of cascade steps: " << average(N_) << endl;
}


void MFM::registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff) 
{

	word name = f.name();

	Info << "register ScalarField " <<  name << endl;

	registeredScalarFields_.insert
	(
		name,
		f
	);

        laminarDiffusivity_.insert
        (
          name,
          volScalarField(
          	IOobject
                (
                        "D_"+name,
                        f.time().timeName(),
                        f.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                ),
                f.mesh(),
                dimensionedScalar("D_"+name,dimensionSet(0,2,-1,0,0,0,0), molecularDiffusivityCoeff)
          )
        );

        deltaScalarFields_.insert
        ( 
            name,
            new volScalarField
            (
                IOobject
                (
                    name+"Delta",
                    f.time().timeName(),
                    f.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                f.mesh(),
                dimensionedScalar(name+"Delta",f.dimensions(), 0.0)
            )
        );

        (*deltaScalarFields_[name]) = f-testFilter(f);

        NScalarFields_.insert
        (
            name,
            new volScalarField
            (
                IOobject
                (
                    "N_"+name,
                    f.time().timeName(),
                    f.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                f.mesh(),
                dimensionedScalar("N_"+name,dimless, 0.0)
            )
        );
}

//Abschluss
bool MFM::read()
{
    if (LESModel::read())
    {
        smoothingFilter.read(coeffDict());
        testFilter.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 } // End namespace LESmodels
 }
 } // End namespace Foam

// ************************************************************************* //
