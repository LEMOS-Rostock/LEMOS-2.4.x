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


Namespace
    Foam::DMDModels

Group
    grpDMDModels

Description
    Namespace for DMD models.

Class
    Foam::DMDModel

Description
    Abstract base class for dynamic decomposition models.

SourceFiles
    DMDModel.C

\*---------------------------------------------------------------------------*/

#ifndef DMDModel_H
#define DMDModel_H

#include "volFields.H"
#include "fvMatrices.H"
#include "Matrix.H"
#include "IOdictionary.H"
#include "complexFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


    void multiply
    (
	scalarSquareMatrix& ans,         // value changed in return
	const scalarRectangularMatrix& A,
	const scalarRectangularMatrix& B,
	const scalarRectangularMatrix& C,
	const DiagonalMatrix<scalar>& D
    );
    
    void multiply
    (
	scalarSquareMatrix& ans,         // value changed in return
	const DiagonalMatrix<scalar>& A,
	const scalarSquareMatrix& B,
	const DiagonalMatrix<scalar>& C
    );




// Forward declarations
class fvMesh;
/*---------------------------------------------------------------------------*\
                           Class PANSModel Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DMDModel
:
    public IOdictionary
{

protected:

    // Protected data


	//- 
	const fvMesh& mesh_;
	
        //- Model coefficients dictionary
        dictionary coeffDict_;

    // Protected Member Functions
    
    template<class Type2>
    void printMatrix(const List<Type2>& matrix, word name, label r = 0, label c = 0, label nRows = 0, label nColumns = 0);
    
    template<class Form, class Type2>
    void printMatrix(const Matrix<Form, Type2>& matrix, word name, label r = 0, label c = 0, label nRows = 0, label nColumns = 0);
    



private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DMDModel(const DMDModel&);

        //- Disallow default bitwise assignment
        void operator=(const DMDModel&);


public:

    //- Runtime type information
    TypeName("DMDModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            DMDModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& DMDModelName
            ),
            (mesh, DMDModelName)
        );


    // Constructors

        //- Construct from components
        DMDModel
        (
            const fvMesh& mesh,
            const word& DMDModelName
        );


    // Selectors

        //- Return a reference to the selected PANS model
        static autoPtr<DMDModel<Type> > New
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~DMDModel()
    {}


    // Member Functions

        // Access
          virtual label baseSize() = 0;

          virtual const  FieldField<Field, complex>& ritzValues() = 0;

          
          virtual void calcBase(PtrList<GeometricField<Type, fvPatchField, volMesh> >& fields) = 0;
        
      //  virtual void calcModes(PtrList<GeometricField<Type, fvPatchField, volMesh> >& fields) = 0;

        //- Read DMDProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDMDModel(Type)                                                    \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        DMDModel<Type>, 		                                      \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        DMDModel<Type>,                                                       \
        dictionary                                                            \
    );                                            


#define makeDMDModelType(SS, Type)					      \
									      \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);			      \
									      \
    DMDModel<Type>::							      \
        adddictionaryConstructorToTable<SS<Type> >                            \
             add##SS##Type##ConstructorToTable_;

//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DMDModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
