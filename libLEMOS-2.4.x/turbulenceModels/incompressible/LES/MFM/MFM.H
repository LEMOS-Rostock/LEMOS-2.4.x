
#ifndef MFM_H
#define MFM_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/LES/lnInclude/LESModel.H"
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

class MFM
:
  //public GenEddyVisc
  virtual public LESModel
  {

private:

  // Disallow default bitwise copy construct and assignment
  MFM(const MFM&);
  MFM& operator=(const MFM&);

  IOdictionary thermoDict;

  autoPtr<LESfilter> smoothingFilterPtr;
  LESfilter& smoothingFilter;

  autoPtr<LESfilter> testFilterPtr;
  LESfilter& testFilter;

  const dictionary LESfilterDict;

  static const label alpha=2;
  static const scalar  phi0=2.0/3.0;

protected:

  volScalarField viscLengthScale_;
  volSymmTensorField Ureynolds_;
  volVectorField uDelta_;

  volScalarField N_;
  volScalarField Re_;
  volScalarField Cai_;

  HashTable<volScalarField&, word> registeredScalarFields_;
  HashTable<volScalarField, word> laminarDiffusivity_;
  HashPtrTable<volScalarField, word> deltaScalarFields_;
  HashPtrTable<volScalarField, word> NScalarFields_;

public:

    //- Runtime type information
    TypeName("MFM");

    // Constructors

        //- Constructors from components
        MFM
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

        //- Constructors incl. scalar sgs
        MFM
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& F,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );

    // Destructor

        ~MFM()
        {}

  // Private Member Functions
  
    //- Update sub-grid scale fields
    void updateSubGridScaleFields(const volSymmTensorField& S);
    void ReS(const volSymmTensorField& S);
    void ReU(const volSymmTensorField& S);

    void updateCai(const volSymmTensorField& S);

    scalar Bsgs() const;
    scalar Dsgs() const;
    
    // Member Functions
    
	// Update 

        //- Return B.
        tmp<volSymmTensorField> B() const;

        //- Return B.
        tmp<volScalarField> D() const;

        tmp<volScalarField> F1() const;

        //- Return B.
        tmp<volSymmTensorField> devReff() const;

        //- Return sub-grid disipation rate
        tmp<volScalarField> epsilon() const;

        //- Return the SGS viscosity.
        tmp<volScalarField> nuSgs() const
        {
           return tmp<volScalarField>
           (
              new volScalarField
              (
                 IOobject
                 (
                    "nuSgs",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                 mesh_,
                 dimensionedScalar("nuSgs", nu()().dimensions(), 0.0)
              )
           );
        };
        // - Return the SGS kinetic energy
        tmp<volScalarField> k() const;

        //- Returns div(B).
        // This is the additional term due to the filtering of the NSE.
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;


        //- Correct Eddy-Viscosity and related properties
        void correct(const tmp<volTensorField>& gradU);


        //- Return number of cascade steps.
        void updateN(const volSymmTensorField& S);


	//- Return the sub-grid scalar flux vector.
	//tmp<volVectorField> F(const volScalarField &f) const;
	
 	//- Return the sub-grid scalar flux vector.
	tmp<volVectorField> Feff(const volScalarField &f) const;

        tmp<fvScalarMatrix> divFeff(volScalarField &f) const;
        tmp<fvScalarMatrix> divFsgs(volScalarField &f) const;

	void registerScalarField(volScalarField &f, scalar molecularDiffusivityCoeff);

        tmp<volScalarField> molecularDiffusivityCoeff(word name) const;

	
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
