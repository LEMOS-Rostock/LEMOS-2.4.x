/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 LEMOS, University Rostock
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    DMDEigenBase

\*---------------------------------------------------------------------------*/

#include "mosDMDEigenBase.H"
#include "volFields.H"
#include "mosDMD.H"
#include "fvMatrices.H"
#include "DMDModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mosDMDEigenBase::realSymmEigenSolver(const  Eigen::MatrixXd& M, Eigen::DiagonalMatrix<scalar, Eigen::Dynamic>& S, Eigen::MatrixXd& U)
{    
    // Solve eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
    eigenSolver.compute(M);
      
    // Sort eigenvalues and corresponding eigenvectors
    // in descending order
    SortableList<scalar> sortedList(M.rows());

    forAll (sortedList, i)
    {
        sortedList[i] = eigenSolver.eigenvalues()[i];
    }
         
    // Do sort 
    sortedList.sort();
      
    label n = 0;
    forAllReverse(sortedList, i)
    {
        S.diagonal()[n] = sortedList[i];
        U.col(n) = eigenSolver.eigenvectors().col(sortedList.indices()[i]);
        
        n++;
    }
}

void Foam::mosDMDEigenBase::realNonsymmEigenSolver(const Eigen::MatrixXd& M, Eigen::DiagonalMatrix<std::complex<scalar>, Eigen::Dynamic>& D, Eigen::MatrixXcd& V)
{
    // Solve eigenvalues and eigenvectors
    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver;
    eigenSolver.compute(M);

    D = eigenSolver.eigenvalues().asDiagonal();
    V = eigenSolver.eigenvectors();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given a list of fields
Foam::mosDMDEigenBase::mosDMDEigenBase(const PtrList<volScalarField>& snapshots, const label nSnapshots)
:
    cumEigenValues_(nSnapshots-1),
    ritzValues_(nSnapshots-1),
    coeffs_(nSnapshots-1),
    H_(Eigen::MatrixXd::Zero(nSnapshots, nSnapshots)),
    modeNorms_(nSnapshots-1)
{
    // Calculate the snapshot of the field with all available fields except the last one
    for (label snapI = 0; snapI < nSnapshots; snapI++)
    {
        for (label snapJ = 0; snapJ <= snapI; snapJ++)
        {
            // Calculate the inner product and insert it into the matrix
            H_(snapI, snapJ) =
                mosDMD::projection
                (
                    snapshots[snapI],
                    snapshots[snapJ]
                );

            if (snapI != snapJ)
            {
                H_(snapJ,snapI) = H_(snapI, snapJ);
            }       
        }
    }

    // 
    Eigen::MatrixXd U(Eigen::MatrixXd::Zero(nSnapshots-1, nSnapshots-1));
    Eigen::DiagonalMatrix<scalar, Eigen::Dynamic> S(Eigen::DiagonalMatrix<scalar, Eigen::Dynamic>(nSnapshots-1));

    // Solve eigen base and retrieve eigenvalues and eigenvectors
    // in descending order
    realSymmEigenSolver(H_.block(0,0,nSnapshots-1,nSnapshots-1), S, U);    
  
    // Calculate square root and inverse of eigenvalues
    S.diagonal() = S.inverse().diagonal().cwiseMax(0.0).cwiseSqrt();
    
    // Calculate matrix M
    Eigen::MatrixXd M = S*U.adjoint()*H_.block(0,1,nSnapshots-1,nSnapshots-1)*U*S;
   
    Eigen::MatrixXcd V(Eigen::MatrixXcd::Zero(M.rows(), M.cols()));

    // Solve eigen base and retrieve eigenvalues and eigenvectors
    Eigen::DiagonalMatrix<std::complex<scalar>, Eigen::Dynamic> eigenValues(M.rows());
    realNonsymmEigenSolver(M, eigenValues, V);

    // Clean up matrix M, not needed anymore
    M.resize(0,0);
    
    // Calculate matrix D
    Eigen::DiagonalMatrix<std::complex<scalar>, Eigen::Dynamic> D(((V.adjoint()*V).inverse()*V.adjoint()*S*U.adjoint()*H_.block(0,0,nSnapshots-1,1)));

    // Calculate build coeff matrix
    Eigen::MatrixXcd buildCoeff =(U*S*V*D);

    // Calculate mode norms
    Eigen::DiagonalMatrix<std::complex<scalar>, Eigen::Dynamic> modeNorms((buildCoeff.adjoint()*H_.block(0,0,nSnapshots-1,nSnapshots-1)*buildCoeff).diagonal());
    
    // Clean up matrix H, not needed anymore
    H_.resize(0,0);
   
    for(label n=0; n<nSnapshots-1; n++)
    {       
        // Grep eigenvalues
        ritzValues_[n].Re() = eigenValues.diagonal()[n].real();
        ritzValues_[n].Im() = eigenValues.diagonal()[n].imag();

        complexField* ptrCoeffs = new complexField(buildCoeff.rows());
    
        forAll(*ptrCoeffs, I)
        {
            // Grep build coeff
            complex& coeff = (*ptrCoeffs)[I];
                                   
            coeff.Re() = buildCoeff(I,n).real();
            coeff.Im() = buildCoeff(I,n).imag();
        }
        
        coeffs_.set
        (
            n,
            ptrCoeffs
        ); 

        scalar norm =  modeNorms.diagonal()[n].real();

	if(norm < 0.0)
        {
            FatalErrorIn
            (
                "Mode norm has negative values"
                "This happens when the rank of input matrix"
                "is much less than the number of columns"
            ) << "Try using fewer input vectors"
            << abort(FatalError);
        }

        modeNorms_[n] = norm;
    } 

    // Clean up matrices 
    V.resize(0,0);
    S.resize(0);
    D.resize(0);
    modeNorms.resize(0);
}

// ************************************************************************* //
