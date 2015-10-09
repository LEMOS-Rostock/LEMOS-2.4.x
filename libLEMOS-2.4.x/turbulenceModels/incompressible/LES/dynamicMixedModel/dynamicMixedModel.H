
#ifndef dynamicMixedModel_H
#define dynamicMixedModel_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/LES/lnInclude/GenEddyVisc.H"
#include "LESfilter.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class dynMixedSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

class dynamicMixedModel
:
   public GenEddyVisc
{

private:

  // Disallow default bitwise copy construct and assignment
  dynamicMixedModel(const dynamicMixedModel&);
  dynamicMixedModel& operator=(const dynamicMixedModel&);


  autoPtr<LESfilter> testFilterPtr_;
  LESfilter& testFilter_;
  autoPtr<LESfilter> gridFilterPtr_;
  LESfilter& gridFilter_;

  const dictionary LESfilterDict;

protected:

  volScalarField cD_;
  volScalarField cI_;
  volScalarField nuSgs_;
  volTensorField Leo_;
  volScalarField k_;
  volTensorField H;
  volVectorField HPhi;
  volTensorField L;
  volTensorField L_m;
  volTensorField L_t;
  volVectorField LPhi;
  volVectorField LPhi_m;
  volVectorField LPhi_t;

  HashTable<volScalarField&, word> registeredScalarFields_;
  HashTable<volScalarField, word> laminarDiffusivity_;
  HashPtrTable<volScalarField, word> turbulentDiffusivity_;
  HashPtrTable<volVectorField,word> LeoPhi_;
  HashPtrTable<volScalarField,word> cF_;

public:

    //- Runtime type information
    TypeName("dynamicMixedModel");

    // Constructors

        //- Constructors from components
        dynamicMixedModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

        //- Constructors incl. scalar sgs
        dynamicMixedModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& F,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

    // Destructor

        ~dynamicMixedModel()
        {}


    // Member Functions

        //- Return B.
        tmp<volSymmTensorField> B() const;

        //- Return B.
        tmp<volSymmTensorField> devReff() const;

        //- Return sub-grid disipation rate
        tmp<volScalarField> epsilon() const;

        //- Return the SGS viscosity.
        tmp<volScalarField> nuSgs() const
        {
	  return nuSgs_;
        };
        // - Return the SGS kinetic energy
        tmp<volScalarField> k() const;

        //- Returns div(B).
        // This is the additional term due to the filtering of the NSE.
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Correct Eddy-Viscosity and related properties
        void correct(const tmp<volTensorField>& gradU);

	//- Return the sub-grid scalar flux vector.
	tmp<volVectorField> F(const volScalarField &f) const;
	
 	//- Return the sub-grid scalar flux vector.
	tmp<volVectorField> Feff(const volScalarField &f) const;

        tmp<fvScalarMatrix> divFeff(volScalarField &f) const;

	void registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff);

        tmp<volScalarField> molecularDiffusivityCoeff(word name) const;

        tmp<volScalarField> turbulentDiffusivityCoeff(word name) const;
	
        //- Read turbulenceProperties dictionary
        bool read();

      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESmodels
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
