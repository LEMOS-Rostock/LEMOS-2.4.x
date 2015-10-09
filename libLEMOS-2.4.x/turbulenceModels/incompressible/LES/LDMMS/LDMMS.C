#include "LDMMS.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
 {
namespace incompressible
 {
namespace LESModels
 {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LDMMS, 0);
addToRunTimeSelectionTable(LESModel, LDMMS, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LDMMS::LDMMS
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
  LESModel(modelName, U, phi, transport, turbulenceModelName),
  GenEddyVisc(U, phi, transport),
  
  filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
  filter_(filterPtr_()),
  gridFilterPtr_(LESfilter::New(U.mesh(), coeffDict())),
  gridFilter_(gridFilterPtr_()), 
  cD_
  (
    IOobject
    (
      "cD_LDMMS",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, 0.0)
    ),
        
  cDsim_
  (
    IOobject
    (
      "cDsim_LDMMS",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("zero", dimless, 0.0)
    ),
  cI_
  (
     IOobject
     (
      "cI_LDMMS",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedScalar("zero", dimless, 0.0)
     ),


  nuSgs_
  (
     IOobject
     (
      "nuSgs",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedScalar("zero", nu()().dimensions(), 0.0)
     ),

  Leo_(
     IOobject
     (
      "leonardStress",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedSymmTensor("0", U.dimensions() * U.dimensions(),
       symmTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
     ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
	dimensionedScalar("zero", U.dimensions() & U.dimensions(), 0.0)
    ),
     LLMM
     (
      IOobject
      (
       "LLMM",
       runTime_.timeName(),
       mesh_,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
       ),
      mesh_
      ),
     NNMM
     (
      IOobject
      (
       "NNMM",
       runTime_.timeName(),
       mesh_,
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
       ),
      mesh_
      ),
      MMMM
      (
       IOobject
       (
        "MMMM",
        runTime_.timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
        ),
       mesh_
       ),
     LL
    (
     IOobject
     (
      "LL",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedSymmTensor("0", U.dimensions() * U.dimensions(),  symmTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     NN
    (
     IOobject
     (
      "NN",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedSymmTensor("0", U.dimensions() * U.dimensions(),  symmTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     MM
    (
     IOobject
     (
      "MM",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh_,
     dimensionedSymmTensor("0", U.dimensions() * U.dimensions(),  symmTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     T
    (
     IOobject
     (
      "T",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
	1/(1.5*delta())*pow(mag(LLMM-NNMM)*MMMM,0.125)
     )
{
  updateSubGridScaleFields(dev(symm(fvc::grad(U)))); 

  printCoeffs();
  
  Info << "LDMMS created" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void LDMMS::updateSubGridScaleFields(const volSymmTensorField& S)
{
  nuSgs_ = cD(S)*sqr(delta())*sqrt(2.0)*mag(S);
  nuSgs_.correctBoundaryConditions();
}

tmp<volScalarField> LDMMS::cD(const volSymmTensorField& S)
{


   //T = 1/(1.5*delta())*pow(8*max((LLMM-NNMM)*MMMM,dimensionedScalar("zero", MMMM.dimensions()*MMMM.dimensions(), 0.0)),0.125);

  //volScalarField LLMM_NNMM = -8*(LLMM-NNMM)*MMMM;
  volScalarField LLMM_NNMM = (LLMM-NNMM)*MMMM;
  LLMM_NNMM.max(0.0);
  LLMM.correctBoundaryConditions(); 
  NNMM.correctBoundaryConditions(); 
  MMMM.correctBoundaryConditions(); 
   T = 1/(1.5*delta())*pow(LLMM_NNMM,0.125);
   //T = pow(-8.0*LLMM*MMMM,0.125)/(1.5*delta());
   //scalarField T1 = (pow(8*(LLMM.internalField()-NNMM.internalField())*MMMM.internalField(),0.125))/(1.5*delta());
   //T = (pow(8*(LLMM-NNMM)*MMMM,0.125))/(1.5*delta());
   //T = 1/(1.5*delta())*(pow(-8*LLMM*MMMM,0.125));
   
   T.correctBoundaryConditions();

   LL = (filter_(sqr(U())) - (sqr(filter_(U()))));

   MM =  2*sqr(delta())*(filter_((sqrt(2.0)*mag(S))*(S)) - 4*(sqrt(2.0)*mag(filter_(S)))*filter_(S)); // Zang
        //sqr((2*delta()))*sqrt(2.0)*mag(filter_(S))*filter_(S) - filter_(sqr(delta())*sqrt(2.0)*mag(S)*S);

   NN = (filter_(gridFilter_(sqr(filter_(U())))) - sqr(filter_(gridFilter_(filter_(U())))))-(filter_(gridFilter_(sqr(U())))-filter_(sqr(gridFilter_(U()))));
   //NN = (filter_(sqr(gridFilter_(U())))-sqr(filter_(gridFilter_((U())))));

  LL.correctBoundaryConditions(); 
  NN.correctBoundaryConditions(); 
  MM.correctBoundaryConditions(); 

   //NN = (filter_(sqr(gridFilter_(U())))-sqr(filter_(gridFilter_(U()))));
  
  
  // transport equation for LL*MM
  fvScalarMatrix LLMMEqn
  (
    fvm::ddt(LLMM)
    +fvm::div(phi() ,LLMM)
    == 
    T*((LL&&MM)-LLMM)
  );

  LLMMEqn.solve();


  // transport equation for MM*MM
  fvScalarMatrix MMMMEqn
  (
    fvm::ddt(MMMM)
    +fvm::div(phi(),MMMM)
    ==
    T*((MM&&MM)-MMMM)
  );

  MMMMEqn.solve();

  // transport equation for NN*MM
  fvScalarMatrix NNMMEqn
  (
    fvm::ddt(NNMM)
    +fvm::div(phi() ,NNMM)
    == 
    T*((NN&&MM)-NNMM)
  );

  NNMMEqn.solve();



  forAll(MMMM, cellid)
  {
    scalar denom = (MMMM[cellid]);

        if (denom  <= VSMALL )
		cD_[cellid]= 0.0;
        else
		cD_[cellid]= ((LLMM[cellid]-NNMM[cellid]))/denom;

        if (cD_[cellid] <= VSMALL )
                LLMM[cellid]=NNMM[cellid];
 } 
  // bounding
  //cD_.min(1.0);
  cD_.max(0.0);
  
  return cD_;
}


tmp<volScalarField> LDMMS::cI(const volSymmTensorField& S) 
{
return cI_;
}

tmp<volScalarField> LDMMS::k() const
{
  return k_;
}

tmp<volSymmTensorField> LDMMS::B() const
{
  return  ((2.0/3.0)*I)*k() - 2*nuSgs_*(symm(fvc::grad(U())));
}

tmp<volSymmTensorField> LDMMS::devReff() const
{
  return  - 2*nuEff()*dev(symm(fvc::grad(U())));
}


tmp<volScalarField> LDMMS::epsilon() const
{
  return
    (
     (  - 2.0*nuSgs_*dev(symm(fvc::grad(U()))))
     &&
     dev(symm(fvc::grad(U())))
     );
}






tmp<fvVectorMatrix> LDMMS::divDevReff(volVectorField& U) const
{
  return
    (
     - fvm::laplacian(nuEff(), U, "laplacian(nuEff,U)")
     - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
     + fvc::div(dev(Leo_))
     );
}

tmp<volVectorField> LDMMS::F(const volScalarField &f) const
{
  return
    (
      *LeoPhi_[f.name()] - (*turbulentDiffusivity_[f.name()]) * fvc::grad(f)
     );
}

tmp<volScalarField> LDMMS::molecularDiffusivityCoeff(word name) const
{
  return
    (
     laminarDiffusivity_[name] 
     );
}

tmp<volScalarField> LDMMS::turbulentDiffusivityCoeff(word name) const
{
  return
    (
     *turbulentDiffusivity_[name]
     );
}

tmp<volVectorField> LDMMS::Feff(const volScalarField &f) const
{
  return
    (
      *LeoPhi_[f.name()] - ( laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()] ) * fvc::grad(f)
     );
}


tmp<fvScalarMatrix> LDMMS::divFeff(volScalarField &f) const
{
  return
     (
      - fvm::laplacian((laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()]), f, "laplacian(Deff,F)")
      + fvc::div((*LeoPhi_[f.name()]))
      );
}

tmp<fvScalarMatrix> LDMMS::divFsgs(volScalarField &f) const
{
  return
     (
      fvm::laplacian((*turbulentDiffusivity_[f.name()]), f, "laplacian(Deff,F)")
      - fvc::div((*LeoPhi_[f.name()]))
      );
}
void LDMMS::correct(const tmp<volTensorField>& gradU)
{

  LESModel::correct(gradU);

  /// strain rate tensor
  volSymmTensorField S = dev(symm(gradU));

  //volSymmTensorField LeoFiltered_ =  (filter_(gridFilter_(sqr(filter_(U()))))-sqr(filter_(gridFilter_(filter_(U())))));
  volSymmTensorField LeoFiltered_ =  ((gridFilter_(sqr((U()))))-sqr((gridFilter_((U())))));

  Leo_ = LeoFiltered_;

  updateSubGridScaleFields(S);

  

  //nuSgs_ = cD_ * sqr(delta()) * sqrt( magSqr(S) );
  //nuSgs_.correctBoundaryConditions();

// 
//   // Calculate dynamic Coefficient for turbulent kinetic energy k, cI_:
 //volScalarField denomCI = 2*sqr(delta())*(filter_((sqrt(2.0)*mag(S))*(S)) - 4*(sqrt(2.0)*mag(filter_(S)))*filter_(S)); //sqr(delta())*(4*sqr(mag(testFilter_(S))) - testFilter_(sqr(mag(S))));
 
 //volScalarField KK = tr(NN);
 
  // forAll(cI_, cellid)
    // {
    //   if(( denomCI[cellid] < SMALL)||(KK[cellid]<SMALL))
// 	{
 //	  cI_[cellid]=0.0;
 //	}
   //    else
 //	{
//	  cI_[cellid]=( KK[cellid] )/denomCI[cellid];
 //	}
 
   //    if(cI_[cellid] < SMALL) cI_[cellid]=0.0;
  //   }
 
 //  k_ = cI_*sqr(delta())*sqr(mag(S));
 //  k_.correctBoundaryConditions();
 
   for ( HashTable<volScalarField&,word>::iterator iter=registeredScalarFields_.begin(); iter!=registeredScalarFields_.end(); iter++ )
   {	
 
       const volScalarField& F = iter();
 
       volVectorField& LeoPhi = *LeoPhi_[F.name()];
       volScalarField& cF = *cF_[F.name()];
       volScalarField& TF = *TF_[F.name()];
 
       volScalarField& LFMF = *LFMF_[F.name()];
       volScalarField& NFMF = *NFMF_[F.name()];
       volScalarField& MFMF = *MFMF_[F.name()];

       volVectorField gradF = fvc::grad(F);
 	
 
	volScalarField LFMF_NFMF = (LFMF-NFMF)*MFMF;
  	LFMF_NFMF.max(0.0);
  	LFMF.correctBoundaryConditions();
  	NFMF.correctBoundaryConditions();
  	MFMF.correctBoundaryConditions();
   	TF = 1/(1.5*delta())*pow(LFMF_NFMF,0.25);


   LeoPhi =  ( gridFilter_(U() * F) ) - (gridFilter_(U()) * ( gridFilter_(F)));
 
   // Calculate directly dynamic Coefficient cF for turbulent viscosity ScSgs by tripple filtering
   volVectorField MF = 2*sqr(delta())*(filter_((sqrt(2.0)*mag(S))*(gradF)) - 4*(sqrt(2.0)*mag(filter_(S)))*filter_(gradF)); // Zang

  volVectorField LF =  filter_(U()*F) - (filter_(U()) * filter_(F));
        
  volVectorField NF =  filter_ ( gridFilter_ ( filter_( U() ) * filter_( F ) ) ) - ( filter_ ( gridFilter_ ( filter_( U() ) ) ) * filter_ ( gridFilter_ ( filter_( F ) ) ) ) -( filter_(gridFilter_(U()*F))-filter_( gridFilter_(U())*gridFilter_(F)));

   // transport equation for LL*MM
  fvScalarMatrix LFMFEqn
  (
    fvm::ddt(LFMF)
    +fvm::div(phi() ,LFMF)
    ==
    TF*((LF&&MF)-LFMF)
  );

  LFMFEqn.solve();

  // transport equation for MM*MM
  fvScalarMatrix MFMFEqn
  (
    fvm::ddt(MFMF)
    +fvm::div(phi(),MFMF)
    ==
    TF*((MF&&MF)-MFMF)
  );

  MFMFEqn.solve();

  // transport equation for NN*MM
  fvScalarMatrix NFMFEqn
  (
    fvm::ddt(NFMF)
    +fvm::div(phi() ,NFMF)
    ==
    TF*((NF&&MF)-NFMF)
  );

  NFMFEqn.solve();


  forAll(MFMF, cellid)
  {
    scalar denom = (MFMF[cellid]);
    scalar num   = (LFMF[cellid]-NFMF[cellid]);

        if (denom  <= 1e-10 )
                cF[cellid]= 0.0;
        else
                cF[cellid]= num/(denom);

        if (cF[cellid] <= SMALL )
                LFMF[cellid]=NFMF[cellid];
 }
  // bounding
  //cF_.min(1.0);
  cF.max(0.0);

      volScalarField& Dt = *turbulentDiffusivity_[F.name()];
 
      forAll(cF, cellid)
          {
 	      scalar dt=0.0;
 
               if(cD_[cellid] > VSMALL)  dt = cF[cellid]*nuSgs_[cellid]/(cD_[cellid]);
 
             if( dt<SMALL)
                 {
                         Dt[cellid] = 0.0;
                 }
                 else
                 {
                         Dt[cellid] = dt;
                 }
 
          }
 
       Dt.correctBoundaryConditions();
     }
}


void LDMMS::registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff) 
{

	Info << "register ScalarField " << endl;
	word name = f.name();


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

        turbulentDiffusivity_.insert
        (
          name,
          new volScalarField (
                IOobject
                (
                        "Dt_"+name,
                        f.time().timeName(),
                        f.mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                ),
                f.mesh(),
                dimensionedScalar("Dt_"+name,dimensionSet(0,2,-1,0,0,0,0), 0.0)
          )
        );


        LeoPhi_.insert
        (
                name,
                new volVectorField
                (
                        IOobject
                        (
                                "leonardScalarFlux_" + name,
                                runTime_.timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedVector("zero",  dimensionSet(0,1,-1,0,0,0,0) * f.dimensions() , pTraits<vector>::zero)
                )
        );

        cF_.insert
        (
                name,
                new volScalarField
                (
                        IOobject
                        (
                                "c" + name + "_LDMMS",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar("zero",  dimless, 0.0)
                )
        );


       LFMF_.insert
        (
                name,
                new volScalarField
                (
                        IOobject
                        (
                                "LFMF" + name + "_LDMMS",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh_
                        //dimensionedScalar("zero",  sqr(U().dimensions())*f.dimensions(), 0.0)
                )
        );

       NFMF_.insert
        (
                name,
                new volScalarField
                (
                        IOobject
                        (
                                "NFMF" + name + "_LDMMS",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh_
                        //dimensionedScalar("zero",  sqr(U().dimensions())*f.dimensions(), 0.0)
                )
        );
	MFMF_.insert
        (
                name,
                new volScalarField
                (
                        IOobject
                        (
                                "MFMF" + name + "_LDMMS",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh_
                        //dimensionedScalar("zero",  sqr(U().dimensions())*f.dimensions(), 0.0)
                )
        );

       volScalarField LFMF = *LFMF_[f.name()];
       volScalarField NFMF = *NFMF_[f.name()];
       volScalarField MFMF = *MFMF_[f.name()];


        TF_.insert
        (
                name,
                new volScalarField
                (
                        IOobject
                        (
                                "TF" + name + "_LDMMS",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                        //mesh_,
			pow(mag(*LFMF_[f.name()]-(*NFMF_[f.name()]))*(*MFMF_[f.name()]),0.25)/(1.5*delta())
                )
        );


}

//Abschluss
bool LDMMS::read()
{
    if (LESModel::read())
    {
        filter_.read(coeffDict());

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
