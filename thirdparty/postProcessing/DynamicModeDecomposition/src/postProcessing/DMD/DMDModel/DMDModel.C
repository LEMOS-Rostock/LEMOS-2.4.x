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

#include "DMDModel.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
template<class Type>
template<class Type2>
void DMDModel<Type>::printMatrix(const List<Type2>& matrix, word name, label r, label c, label nRows, label nColumns)
{

    Info << endl;
    Info << endl;
    Info << "Matrix Coeffs of " << name << endl;
    Info << endl;
    
    
    if (r+nRows > matrix.size() || nRows < 1) 
    {
        nRows = matrix.size() - r;
    }
    
    if (c+nColumns > matrix.size() || nColumns < 1)
    {
        nColumns = matrix.size() - c;
    }
    
    for (label i = r; i< (r+nRows); i++)
    {
	for (label j = c; j<(c+nColumns); j++)
	{
	    if (i == j)
	    {
	        Info << setw(8) << matrix[i] << " ";
	    }
	    else
	    {
	        Info << setw(8) << "0" << " ";
	    }
	}
	
	Info << endl;
    }
    
    Info << endl;
}



template<class Type>
template<class Form, class Type2>
void DMDModel<Type>::printMatrix(const Matrix<Form, Type2>& matrix, word name, label r, label c, label nRows, label nColumns)
{
    if (r >= matrix.n() || r < 0)
    {
        FatalErrorIn
        (
            "printMatrix("
            "Matrix<Form, Type2>& matrix,"
            "word name, "
            "label r, "
            "label c, "
            "label rows, "
            "label columns "
        )   << "r is greater than the number of rows in matrix! "
	    << "r must have a value between 0 and matrix.n()-1"
            << abort(FatalError);
    }
    
    if (c >= matrix.m() || c < 0)
    {
        FatalErrorIn
        (
            "printMatrix("
            "Matrix<Form, Type2>& matrix,"
            "word name, "
            "label r, "
            "label c, "
            "label rows, "
            "label columns "
        )   << "c is greater than the number of columns in matrix! "
	    << "c must have a value between 0 and matrix.m()-1"
            << abort(FatalError);
    }
    
    if (nRows < 0)
    {
        FatalErrorIn
        (
            "printMatrix("
            "Matrix<Form, Type2>& matrix,"
            "word name, "
            "label r, "
            "label c, "
            "label rows, "
            "label columns "
        )   << "nRows must have a value greater or equal 0 "
            << abort(FatalError);
    }
    
    if (nColumns < 0)
    {
        FatalErrorIn
        (
            "printMatrix("
            "Matrix<Form, Type2>& matrix,"
            "word name, "
            "label r, "
            "label c, "
            "label rows, "
            "label columns "
        )   << "nColumns must have a value greater or equal 0 "
            << abort(FatalError);
    }

    Info << endl;
    Info << endl;
    Info << "Matrix Coeffs of " << name << endl;
    Info << endl;
    
    
    if (r+nRows > matrix.n() || nRows < 1) 
    {
        nRows = matrix.n() - r;
    }
    
    if (c+nColumns > matrix.m() || nColumns < 1)
    {
        nColumns = matrix.m() - c;
    }
    
    for (label i = r; i< (r+nRows); i++)
    {
	for (label j = c; j<(c+nColumns); j++)
	{
	    Info << setw(8) << matrix[i][j] << " ";
	}
	
	Info << endl;
    }
    
    Info << endl;
}



void multiply
(
    scalarSquareMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B,
    const scalarRectangularMatrix& C,
    const scalarDiagonalMatrix& D
)
{
    if (A.m() != B.n())
    {
        FatalErrorIn
        (
            "multiply("
            "scalarRectangularMatrix& answer),"
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "const DiagonalMatrix<scalar>& D"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.n()
            << abort(FatalError);
    }
    
    if (B.m() != C.n())
    {
        FatalErrorIn
        (
              "multiply("
            "scalarRectangularMatrix& answer),"
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "const DiagonalMatrix<scalar>& D"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.m() << " and C.n = " << C.n()
            << abort(FatalError);
    }
    
    
    if (C.m() != D.size())
    {
        FatalErrorIn
        (
            "multiply("
            "scalarRectangularMatrix& answer),"
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "const DiagonalMatrix<scalar>& D"
        )   << "C and D must have identical inner dimensions but C.m = "
            << C.m() << " and D.n = " << D.size()
            << abort(FatalError);
    }


    ans = scalarSquareMatrix(D.size(), D.size(), scalar(0));

    for (register label i = 0; i < A.n(); i++)
    {
        for (register label g = 0; g < C.m(); g++)
        {
            for (register label l = 0; l < C.n(); l++)
            {
                scalar ab = 0;
                for (register label j = 0; j < A.m(); j++)
                {   
                    ab += A[i][j]*B[j][l];
                }
                ans[i][g] += C[l][g] * ab;
            }
            
            ans[i][g] = ans[i][g] * D[g];
        }
    }

}


void multiply
(
    scalarSquareMatrix& ans,         // value changed in return
    const scalarDiagonalMatrix& A,
    const scalarSquareMatrix& B,
    const scalarDiagonalMatrix& C
)
{
    if (A.size() != B.n())
    {
        FatalErrorIn
        (
            "multiply("
            "scalarRectangularMatrix& answer),"
            "const DiagonalMatrix<scalar>& A, "
            "const scalarRectangularMatrix& B, "
            "const DiagonalMatrix<scalar>& C"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.size() << " and B.n = " << B.n()
            << abort(FatalError);
    }
    
    if (B.m() != C.size())
    {
        FatalErrorIn
        (
             "multiply("
            "scalarRectangularMatrix& answer),"
            "const DiagonalMatrix<scalar>& A, "
            "const scalarRectangularMatrix& B, "
            "const DiagonalMatrix<scalar>& C"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.m() << " and C.n = " << C.size()
            << abort(FatalError);
    }
    
    
    ans = scalarSquareMatrix(B.n(), B.n(), scalar(0));

    for (register label i = 0; i < A.size(); i++)
    {
        for (register label j = 0; j < C.size(); j++)
        {               
            ans[i][j] = A[i] * B[i][j] * C[j];
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type>
DMDModel<Type>::DMDModel
(
    const fvMesh& mesh,
    const word& DMDModelName
)
:
    IOdictionary
    (
        IOobject
        (
            "DMDProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),
    mesh_(mesh),
    coeffDict_(subOrEmptyDict(DMDModelName + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
autoPtr<DMDModel<Type> > 
DMDModel<Type>::New
(
   const fvMesh& mesh
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word DMDModelName
    (
        IOdictionary
        (
            IOobject
            (
                "DMDProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("DMDModel")
    );

    Info<< "Selecting DMD model " << DMDModelName << endl;

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "DMDModel::New"
            "("
                "const fvMesh& mesh,"
                "const PtrList<GeometricField<Type, fvPatchField, volMesh> >& fields"
            ")"
        )   << "DMD model table is empty"
            << exit(FatalError);
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(DMDModelName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DMDModel::New"
            "("
                "const fvMesh& mesh,"
                "const PtrList<GeometricField<Type, fvPatchField, volMesh> >& fields"
            ")"
        )   << "Unknown DMDModel type "
            << DMDModelName << nl << nl
            << "Valid DMDModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Foam::DMDModel<Type> >
    (
        cstrIter()(mesh, DMDModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
bool DMDModel<Type>::read()
{
    if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
    {
        coeffDict_ <<= *dictPtr;
        
        return true;
    }
    else
    {
        return false;
    }
}


}
// ************************************************************************* //
