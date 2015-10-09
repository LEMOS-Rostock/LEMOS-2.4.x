/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "DimensionedField.H"
#include "fvMatrices.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponential, 0);

    addToRunTimeSelectionTable(bodyForceModel, exponential, dictionary);
}


void Foam::exponential::calcForceField()
{
  scalar magD=mag(D_);

  if (magD < SMALL)
  {
      WarningIn("Foam::fv::plasmaActuatorBodyForce::exponential::calcAxialLength() const")
          <<"D vector has zero length!"<<endl;
  }

  vector nD = D_/magD;

  vector min = boundBox(UIndirectList<vector>(mesh_.C(), cells_)()).min();
  vector max = boundBox(UIndirectList<vector>(mesh_.C(), cells_)()).max();

  vectorField cellCenters(UIndirectList<vector>(mesh_.C(), cells_)());

  forAll(cellCenters, cellI)
  {
      const vector& cf = cellCenters[cellI];

      vector tcf;
      tcf.x() =  (min.x() - cf.x())/ (min.x() - max.x());
      tcf.y() =  (min.y() - cf.y())/ (min.y() - max.y());
      tcf.z() =  (min.z() - cf.z())/ (min.z() - max.z());

      scalar X = (Cx_[0]+Cx_[1]*tcf.x()+Cx_[2]*sqr(tcf.x()))*exp(-1.0*Cx_[3]*pow(tcf.x(), Cx_[4]));
      scalar Y = (Cy_[0]+Cy_[1]*tcf.y()+Cy_[2]*sqr(tcf.y()))*exp(-1.0*Cy_[3]*pow(tcf.y(), Cy_[4]));
      scalar Z = (Cz_[0]+Cz_[1]*tcf.z()+Cz_[2]*sqr(tcf.z()))*exp(-1.0*Cz_[3]*pow(tcf.z(), Cz_[4]));

       Su_()[cellI] = nD*(nD&fMax_)*X*Y*Z;
  }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exponential::exponential
(
    const fv::plasmaActuationBodyForce& dbd,
    const dictionary& dict,
    const fvMesh& mesh,
    const labelList& cells
)
:
    bodyForceModel(dbd, dict, typeName, mesh, cells),
    Su_(new vectorField(cells.size(), pTraits<vector>::zero)),
    correct_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::exponential::read(const dictionary& dict)
{
//  if (bodyForceModel::read(dict))
//  {
    coeffs_.lookup("D") >> D_;
    coeffs_.lookup("Cx") >> Cx_;
    coeffs_.lookup("Cy") >> Cy_;
    coeffs_.lookup("Cz") >> Cz_;
    coeffs_.lookup("fMax") >> fMax_;

    bodyForceModel::read(dict);
    
    return true;
//  }
//  else
//  { 
//    return false;
//  }
}


Foam::tmp<Foam::vectorField> Foam::exponential::computeSup(fvMatrix<vector>& eqn)
{
    if(correct_)
    {
        calcForceField();

        correct_ = false;
    }

    return Su_;
}


/*
void Foam::exponential::correct
(
    const vectorField& U,
    vectorField& force
)
{
    // do nothing
}


void Foam::exponential::correct
(
    const volScalarField rho,
    const vectorField& U,
    vectorField& force)
{
    // do nothing
}
*/

// ************************************************************************* //
