//#include <stdio.h>
#include "dynamicMixedModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
 {
namespace incompressible
 {
namespace LESModels
 {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicMixedModel, 0);
addToRunTimeSelectionTable(LESModel, dynamicMixedModel, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicMixedModel::dynamicMixedModel
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

  testFilterPtr_(LESfilter::New(U.mesh(), coeffDict())),
  testFilter_(testFilterPtr_()),

  gridFilterPtr_(LESfilter::New(U.mesh(), coeffDict())),
  gridFilter_(gridFilterPtr_()),

   LESfilterDict(coeffDict()),
  cD_
  (
     IOobject
     (
      "cD_dynamicMixedModel",
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
      "cI_dynamicMixedModel",
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
     dimensionedTensor("0", U.dimensions() * U.dimensions(),
       tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
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
    H
    (
     IOobject
     (
      "H",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedTensor("0", U.dimensions() * U.dimensions(),  tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     HPhi
    (
     IOobject
     (
      "HPhi",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedVector("0", U.dimensions() ,  pTraits<vector>::zero)
     ),
    L
    (
     IOobject
     (
      "L",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedTensor("0", U.dimensions() * U.dimensions(),  tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     L_m
    (
     IOobject
     (
      "L_m",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedTensor("0", U.dimensions() * U.dimensions(),  tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
     L_t
    (
     IOobject
     (
      "L_t",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedTensor("0", U.dimensions() * U.dimensions(),  tensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
     ),
    LPhi
    (
     IOobject
     (
      "LPhi",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedVector("0", U.dimensions() ,  pTraits<vector>::zero)
     ),
     LPhi_m
    (
     IOobject
     (
      "LPhi_m",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedVector("0", U.dimensions() ,  pTraits<vector>::zero)
     ),
     LPhi_t
    (
     IOobject
     (
      "LPhi_t",
      runTime_.timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh_,
     dimensionedVector("0", U.dimensions() ,  pTraits<vector>::zero)
     )
{
Info << "DMM created" << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> dynamicMixedModel::k() const
{
  return k_;
}

tmp<volSymmTensorField> dynamicMixedModel::B() const
{
  return symm(Leo_) + ((2.0/3.0)*I)*k() - 2*nuSgs_*(symm(fvc::grad(U())));
}

tmp<volSymmTensorField> dynamicMixedModel::devReff() const
{
  return dev(symm(Leo_)) - 2*nuEff()*dev(symm(fvc::grad(U())));
}


tmp<volScalarField> dynamicMixedModel::epsilon() const
{
  return
    (
     ( Leo_ - 2.0*nuSgs_*dev(symm(fvc::grad(U())))) 
     && 
     dev(symm(fvc::grad(U())))
     );
}


tmp<fvVectorMatrix> dynamicMixedModel::divDevReff(volVectorField& U) const
{
  return
    (
     - fvm::laplacian(nuEff(), U, "laplacian(nuEff,U)")
     - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
     + fvc::div(dev(Leo_))
     );
}

tmp<volVectorField> dynamicMixedModel::F(const volScalarField &f) const
{
  return
    (
     *LeoPhi_[f.name()] - (*turbulentDiffusivity_[f.name()]) * fvc::grad(f)
     );
}

tmp<volScalarField> dynamicMixedModel::molecularDiffusivityCoeff(word name) const
{
  return
    (
     laminarDiffusivity_[name] 
     );
}

tmp<volScalarField> dynamicMixedModel::turbulentDiffusivityCoeff(word name) const
{
  return
    (
     *turbulentDiffusivity_[name]
     );
}

tmp<volVectorField> dynamicMixedModel::Feff(const volScalarField &f) const
{
  return
    (
     *LeoPhi_[f.name()] - ( laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()] ) * fvc::grad(f)
     );
}


tmp<fvScalarMatrix> dynamicMixedModel::divFeff(volScalarField &f) const
{
  return
     (
      fvm::laplacian((laminarDiffusivity_[f.name()] + *turbulentDiffusivity_[f.name()]), f, "laplacian(Deff,F)")
      - fvc::div((*LeoPhi_[f.name()]))
      );
}

void dynamicMixedModel::correct(const tmp<volTensorField>& gradU)
{

  LESModel::correct(gradU);

  
  word calcLeonardStress_(LESfilterDict.lookup("calcLeonardStress"));

  word calcDynamicConstant_(LESfilterDict.lookup("calcDynamicConstant"));

  volSymmTensorField D = dev( symm( gradU() ) );

  //U field filtering operations
  volVectorField test_filter_U = testFilter_(U());
  volVectorField grid_filter_U = gridFilter_(U());
  volVectorField grid_test_filter_U = gridFilter_(test_filter_U);
  volVectorField test_grid_test_filter_U = testFilter_(grid_test_filter_U);

  volTensorField test_filter_sqrU = testFilter_(U()*U());
  volTensorField test_grid_filter_sqrU = testFilter_(gridFilter_(U()*U()));
  volTensorField test_sqr_grid_filter_U = testFilter_(grid_filter_U*grid_filter_U);
  volTensorField test_grid_filter_sqr_test_filter_U = testFilter_(gridFilter_(test_filter_U * test_filter_U));

  volTensorField LeoFiltered_ =  test_grid_filter_sqr_test_filter_U - ( test_grid_test_filter_U  *  test_grid_test_filter_U) ;

  if (calcLeonardStress_ == "direct")
  {
	// Leonard stress tensor, directly calculated by tripple filtering
	Leo_ =  LeoFiltered_;
  }
  else{

	// Leonard stress tensor, calculated by clark approximation
	L = test_filter_sqrU - test_filter_U * test_filter_U;
	L_t =  LeoFiltered_;
	//}
	//
	// apply clipping
	forAll(L, cellI)
	{
		for (int j=0;j<9;j++)
	 	{
	 		if ( fabs( L_t[cellI].component(j) ) > 1.25 * fabs( L[cellI].component(j) ) )
	 			L_t[cellI].component(j) = 1.25 * L[cellI].component(j);
	 	}
	}

	// assign values
	Leo_ =  L_t ;
  }
  
  if (calcDynamicConstant_ == "direct")
  {

        // Calculate directly dynamic Coefficient cD_ for turbulent viscosity nusgs by tripple filtering

	volSymmTensorField M = -2.0 * testFilter_( sqr( delta() ) ) * mag( testFilter_( D ) ) * testFilter_( D ) + 2.0 * sqr( delta() ) * testFilter_( mag( D ) * D );

	volTensorField L = test_filter_sqrU - test_filter_U * test_filter_U;

	H =  LeoFiltered_ - ( test_grid_filter_sqrU  - test_sqr_grid_filter_U) ;

        volTensorField devLH = dev(L) - dev(H);

        forAll(cD_, cellid)
        {
        	scalar denom=M[cellid] && M[cellid];
         	if( denom<SMALL)
         	{
         		cD_[cellid]=0.0;
         	}
		else
		{
			cD_[cellid]=( M[cellid] && devLH[cellid] )/denom;
		}
              //if( cD_[cellid]<SMALL) cD_[cellid]=0.0;
	}
	
	// bounding
	cD_.min(1.0);
	cD_.max(0.0);	
  }
  else
  {
         // Calculate dynamic Coefficient for turbulent viscosity nusgs, cD_:
	  L = test_filter_sqrU - test_filter_U * test_filter_U;
	  L_m = test_grid_filter_sqrU - test_sqr_grid_filter_U;
	  L_t = LeoFiltered_;
	// }
	 // apply clipping
	 forAll(L, cellI)
	 {
	 	for (int j=0;j<9;j++)
	 	{
	 		if ( fabs( L_t[cellI].component(j) ) > 1.25 * fabs( L[cellI].component(j) ) )
	 			L_t[cellI].component(j) = 1.25 * L[cellI].component(j);

	 		if ( fabs( L_m[cellI].component(j) ) > 0.25 * fabs( L[cellI].component(j) ) )
	 			L_m[cellI].component(j) = 0.25 * L[cellI].component(j);
	 	}
	 }
	
	 H = L_t - L_m;


        volSymmTensorField M = -2.0 * testFilter_( sqr( delta() ) ) * mag( testFilter_( D ) ) * testFilter_( D ) + 2.0 * sqr( delta() ) * testFilter_( mag( D ) * D );


	 volTensorField devLH = dev(L) - dev(H);

         forAll(cD_, cellid)
         {
                scalar denom=M[cellid] && M[cellid];

		scalar numerator = M[cellid] && devLH[cellid];

                //if(( denom<SMALL)||(numerator<SMALL))
                if( denom<SMALL)
                {
                        cD_[cellid]=0.0;
                }
                else
                {
                        cD_[cellid]=numerator/denom;
                } 

              //if( cD_[cellid]<SMALL) cD_[cellid]=0.0;

         }

        // bounding
        cD_.min(1.0);
        cD_.max(0.0); 
  }
	


  nuSgs_ = cD_ * sqr( delta() ) * sqrt( magSqr( D ) );
  nuSgs_.correctBoundaryConditions();



  // Calculate dynamic Coefficient for turbulent kinetic energy k, cI_:
  volScalarField denomCI = sqr(delta())*(4*sqr(mag(testFilter_(D))) - testFilter_(sqr(mag(D))));

  volScalarField KK = tr(H);

  forAll(cI_, cellid)
    {
      if(( denomCI[cellid] < SMALL)||(KK[cellid]<SMALL))
	{
	  cI_[cellid]=0.0;
	}
      else
	{
	  cI_[cellid]=( KK[cellid] )/denomCI[cellid];
	}

      if(cI_[cellid] < SMALL) cI_[cellid]=0.0;
    }

  k_ = cI_*sqr(delta())*sqr(mag(D))+tr(Leo_);
  k_.correctBoundaryConditions();

  for ( HashTable<volScalarField&,word>::iterator iter=registeredScalarFields_.begin(); iter!=registeredScalarFields_.end(); iter++ )
  {	

      const volScalarField& F = iter();

      volVectorField& LeoPhi = *LeoPhi_[F.name()];
      volScalarField& cF = *cF_[F.name()];

      volVectorField gradF = fvc::grad(F);
	
  //F field filtering operations
  volScalarField test_filter_F = testFilter_(F);
  volScalarField grid_filter_F = gridFilter_(F);
  volScalarField grid_test_filter_F = gridFilter_(test_filter_F);
  volScalarField test_grid_test_filter_F = testFilter_(grid_test_filter_F);

  volVectorField test_filter_U_F = testFilter_(U()*F);
  volVectorField test_grid_filter_U_F = testFilter_(gridFilter_(U()*F));
  volVectorField test_sqr_grid_filter_U_F = testFilter_(grid_filter_U*grid_filter_F);
  volVectorField test_grid_filter_test_filter_U_F = testFilter_(gridFilter_(test_filter_U * test_filter_F));

  volVectorField LeoPhiFiltered_ =  test_grid_filter_test_filter_U_F - ( test_grid_test_filter_U  *  test_grid_test_filter_F) ;

	if (calcLeonardStress_ == "direct")
  	{
        	// Leonard stress tensor, directly calculated by tripple filtering
        	// LeoPhi =  testFilter_ ( gridFilter_ ( testFilter_( U() ) * testFilter_( F ) ) ) - ( testFilter_ ( gridFilter_ ( testFilter_( U() ) ) ) * testFilter_ ( gridFilter_ ( testFilter_( F ) ) ) ) ;
        	LeoPhi = LeoPhiFiltered_;
        }
        else
        {
        	// Leonard stress tensor, calculated by clark approximation
	        LPhi = test_filter_U_F - test_filter_U * test_filter_F;
         	LPhi_t = LeoPhiFiltered_;
		//}
                // apply clipping
                forAll(LPhi, cellI)
		{
        		for (int j=0;j<3;j++)
        		{
        			if ( fabs( LPhi_t[cellI].component(j) ) > 1.25 * fabs( LPhi[cellI].component(j) ) )
        				LPhi_t[cellI].component(j) = 1.25 * LPhi[cellI].component(j);
        		}
        	}
        
        	// assign values
        	LeoPhi =  LPhi_t ;
       	}                                        
        
 if (calcDynamicConstant_ == "direct")
  {

        // Calculate directly dynamic Coefficient cF for turbulent viscosity ScSgs by tripple filtering
        volVectorField M = -2.0 * testFilter_( sqr( delta() ) ) * mag( testFilter_( D ) ) * testFilter_( gradF ) + 2.0 * sqr( delta() ) * testFilter_( mag( D ) * gradF );
        
        volVectorField L =  test_filter_U_F - test_filter_U * test_filter_F;
        
        HPhi =  LeoPhiFiltered_ - ( test_grid_filter_U_F  - test_sqr_grid_filter_U_F );
        
        volVectorField LHPhi = L - HPhi;
        
        forAll(cF, cellid)
        {
        	scalar denom=M[cellid] && M[cellid];

	       scalar numerator=M[cellid] && LHPhi[cellid];

                if(( denom<SMALL)||(numerator<SMALL))
        	{
        		cF[cellid]=0.0;
        	}
        	else
        	{
        		cF[cellid]= numerator/denom;
                }
		if(cF[cellid]<SMALL) cF[cellid]=0.0;	
        }
        
	// bounding
        cF.min(1.0);
        cF.max(0.0);
  }
else
  {
	// Calculate dynamic Coefficient for turbulent viscosity ScSgs, cF by clark approximation
	 LPhi = test_filter_U_F - test_filter_U * test_filter_F;
	 LPhi_m = test_grid_filter_U_F -test_sqr_grid_filter_U_F;
         LPhi_t = LeoPhiFiltered_;
	//}	
         // apply clipping
         forAll(LPhi, cellI)
         {
         	for (int j=0;j<3;j++)
         	{
			if ( fabs( LPhi_t[cellI].component(j) ) > 1.25 * fabs( LPhi[cellI].component(j) ) )
         			LPhi_t[cellI].component(j) = 1.25 * LPhi[cellI].component(j);
         
         		if ( fabs( LPhi_m[cellI].component(j) ) > 0.25 * fabs( LPhi[cellI].component(j) ) )
         			LPhi_m[cellI].component(j) = 0.25 * LPhi[cellI].component(j);
                }
         }
         
         HPhi = LPhi_t - LPhi_m;
         
         
         volVectorField M = -2.0 * testFilter_( sqr( delta() ) ) * mag( testFilter_( D ) ) * testFilter_( gradF ) + 2.0 * sqr( delta() ) * testFilter_( mag( D ) * gradF );
        
 
         volVectorField LHPhi = LPhi - HPhi;
         
         forAll(cF, cellid)
         {
         	scalar denom=M[cellid] && M[cellid];

		scalar numerator=M[cellid] && LHPhi[cellid];

         	if(( denom<SMALL)||(numerator<SMALL))
         	{
         		cF[cellid]=0.0;
         	}
         	else
         	{
         		cF[cellid]=numerator/denom;
         	}

		if(cF[cellid]<SMALL) cF[cellid]=0.0;	


         }
         
         // bounding
         cF.min(1.0);
         cF.max(0.0);
  }



     volScalarField& Dt = *turbulentDiffusivity_[F.name()];


     forAll(cF, cellid)
         {
	      scalar dt=0.0;

              if(cD_[cellid] > SMALL)  dt = cF[cellid]*nuSgs_[cellid]/(cD_[cellid]);

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


void dynamicMixedModel::registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff) 
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
                                "c" + name + "_dynamicMixedModel",
                                runTime_.timeName(),
                                mesh_,
                                IOobject::NO_READ,
                                IOobject::AUTO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar("zero",  dimless, 0.0)
                )
        );

}

//Abschluss
bool dynamicMixedModel::read()
{
    if (LESModel::read())
    {
      testFilter_.read(coeffDict());
      gridFilter_.read(coeffDict());
      
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
