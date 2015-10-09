
#ifndef LDMMS_H
#define LDMMS_H

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

class LDMMS
:
   public GenEddyVisc
{

private:

  // Disallow default bitwise copy construct and assignment
  LDMMS(const LDMMS&);
  LDMMS& operator=(const LDMMS&);


  autoPtr<LESfilter> filterPtr_;
  LESfilter& filter_;

  autoPtr<LESfilter> gridFilterPtr_;
  LESfilter& gridFilter_;
  
   const dictionary LESfilterDict;

protected:

  volScalarField cD_;
  volScalarField cDsim_;

  //volScalarField cD_;
  volScalarField cI_;
  volScalarField nuSgs_;
  volSymmTensorField Leo_;
  volScalarField k_;
  
  volScalarField LLMM;
  volScalarField NNMM;
  volScalarField MMMM;
  volSymmTensorField LL;
  volSymmTensorField NN;
  volSymmTensorField MM;
  
  
  volScalarField T;

  HashTable<volScalarField&, word> registeredScalarFields_;
  HashTable<volScalarField, word> laminarDiffusivity_;
  HashPtrTable<volScalarField, word> turbulentDiffusivity_;
  HashPtrTable<volVectorField,word> LeoPhi_;
  HashPtrTable<volScalarField,word> cF_;

  HashPtrTable<volScalarField,word> LFMF_;
  HashPtrTable<volScalarField,word> NFMF_;
  HashPtrTable<volScalarField,word> MFMF_;
  HashPtrTable<volScalarField,word> TF_;


public:

    //- Runtime type information
    TypeName("LDMMS");

    // Constructors

        //- Constructors from components
        LDMMS
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

        //- Constructors incl. scalar sgs
        LDMMS
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& F,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

    // Destructor

        ~LDMMS()
        {}

  // Private Member Functions
  
    //- Update sub-grid scale fields
    void updateSubGridScaleFields(const volSymmTensorField& S);

    //- Calculate coefficients cD, cI from filtering velocity field
    tmp<volScalarField> cD(const volSymmTensorField& S);
    tmp<volScalarField> cI(const volSymmTensorField& S);






    // Member Functions
    
	// Update 

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
        tmp<fvScalarMatrix> divFsgs(volScalarField &f) const;

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
