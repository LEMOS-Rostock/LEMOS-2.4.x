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


Description
    Class which assembles scalar eigen-base given a set of scalar fields.
    This is performed by establishing a matrix of snapshots, calculating and
    sorting them and providing corresponding eigen-vectors.  Eigen-values are
    sorted in increasing order.  Snapshots are added one at a time and
    when the base is completed, the calculation is triggered

Class
    Foam::DMDEigenBase

SourceFiles
    DMDModel.C

\*---------------------------------------------------------------------------*/

#ifndef mosDMDEigenBase_H
#define mosDMDEigenBase_H

#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "complexFields.H"

#include <eigen3/Eigen/Eigenvalues>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                       Class DMDEigenBase Declaration
\*---------------------------------------------------------------------------*/

class mosDMDEigenBase
{
    // Private data

        //- Cumulative relative eigen-values.  Renormalised to sum to 1
        scalarField cumEigenValues_;

        //- Complex Ritz values
        complexField ritzValues_;

        //- Build coeffs for DMD modes
        FieldField<Field, complex> coeffs_;
        
        Eigen::MatrixXd H_;

        //- Squared L2 norm of modes (amplitudes)        
        scalarField modeNorms_;
        
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        mosDMDEigenBase(const mosDMDEigenBase&);

        //- Disallow default bitwise assignment
        void operator=(const mosDMDEigenBase&);


        //- Calculate eigen base (for symmetric real matrices)
        void realSymmEigenSolver(const Eigen::MatrixXd& M, Eigen::DiagonalMatrix<scalar, Eigen::Dynamic>& S, Eigen::MatrixXd& U);
        
         //- Calculate eigen base
        void realNonsymmEigenSolver(const Eigen::MatrixXd& M, Eigen::DiagonalMatrix<std::complex<scalar>, Eigen::Dynamic>& D, Eigen::MatrixXcd& V);


public:

    // Constructors

        //- Construct given a list of fields
        mosDMDEigenBase(const PtrList<volScalarField>& snapshots, const label nSnapshots);


    // Destructor - default


    // Member Functions

        //- Return cumulative eigen-values
        const scalarField& cumulativeEigenValues() const
        {
            return cumEigenValues_;
        }

         //- Return ritz values
         const complexField& ritzValues() const
         {
             return ritzValues_;
         }

        //- Return build coeffs
        const FieldField<Field, complex>& coeffs() const
        {
            return coeffs_;
        }
        
        //- Return mode norms
        const scalarField& modeNorms() const
        {
            return modeNorms_;
        }
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mosDMDEigenBase.C"
#endif

#endif

// ************************************************************************* //
