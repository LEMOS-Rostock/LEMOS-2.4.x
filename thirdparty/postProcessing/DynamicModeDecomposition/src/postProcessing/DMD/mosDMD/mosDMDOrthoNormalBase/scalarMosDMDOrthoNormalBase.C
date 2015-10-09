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


Class
    DMDOrthoNormalBase

Description
    Establish DMD ortho-normal base and interpolation coefficients give a list
    of fields.

\*---------------------------------------------------------------------------*/

#include "scalarMosDMDOrthoNormalBase.H"
#include "mosDMD.H"
#include "mosDMDEigenBase.H"
#include "IOmanip.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::mosDMDOrthoNormalBase<Foam::scalar>::calcOrthoBase
(
    const PtrList<volScalarField>& snapshots
)
{
    const label nSnapshots = snapshots.size();

    // Calculate eigen base for each component
    mosDMDEigenBase eigenBase(snapshots, nSnapshots);


    // Compute selected base size or entire otherwise
    if ((baseSize_ <= 0) || (baseSize_ > nSnapshots-1))
    {
        baseSize_ = nSnapshots-1;
    }

    Info << "Base size: " << baseSize_ << endl;


    // Compute mode energies
    SortableList<scalar> sortedList(nSnapshots-1);
    scalar T = (nSnapshots-1)*deltaT_;
    
    for (label baseI = 0; baseI < nSnapshots-1; baseI++)
    {
        const scalarField& amplitudes = eigenBase.modeNorms();

        if(modeSelection_ == "energy")
        {
            const complexField& ritzValues = eigenBase.ritzValues();

            std::complex<scalar> tmp;
            tmp.real() = ritzValues[baseI].Re();
            tmp.imag() = ritzValues[baseI].Im();
            tmp = std::log(tmp+SMALL)/(deltaT_*2.0*constant::mathematical::pi);
        
            sortedList[baseI] = amplitudes[baseI]*(Foam::exp(2.0*tmp.real()*T)-1)/(2.0*tmp.real()*T+SMALL);
        }
        else
        {
            sortedList[baseI] = sqr(amplitudes[baseI]);
        }
    }

    // Do sort 
    sortedList.sort();
    modeEnergies_.resize(nSnapshots-1); 

    scalarField indicesModeEnergies(nSnapshots-1, 0.0);
    label n = 0;

    forAllReverse(sortedList, i)
    {
        modeEnergies_[n] = sortedList[i];
        indicesModeEnergies[n] = sortedList.indices()[i];

        n++;
    }


    // Establish base
    for (label baseI = 0; baseI < baseSize_; baseI++)
    {
        volScalarField* onBaseRePtr
        (
            new volScalarField
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + Foam::name(baseI) + "Real",
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensionedScalar
                (
                    "zero",
                    snapshots[0].dimensions(),
                    0.0
                )
            )
        );
        volScalarField& onBaseRe = *onBaseRePtr;
        
        volScalarField* onBaseImPtr
        (
            new volScalarField
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + Foam::name(baseI) + "Imag",
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensionedScalar
                (
                    "zero",
                    snapshots[0].dimensions(),
                    0.0
                )
            )
        );
        volScalarField& onBaseIm = *onBaseImPtr;


        // Valid component, grab complex build coeffs
        const complexField& coeffs =
            eigenBase.coeffs()[indicesModeEnergies[baseI]];

        forAll (coeffs, coeffI)
        {
            onBaseRe +=
                coeffs[coeffI].Re()*snapshots[coeffI];
                        
            onBaseIm +=
                coeffs[coeffI].Im()*snapshots[coeffI];
        }

        onBaseRe.write();
        onBaseIm.write();
    }

    // Avoid multiple output of eigen values
    // only master process should write out
    if(Pstream::master())
    {
        fileName outputDir;
     
        if (Pstream::parRun())
        {
            outputDir = snapshots[0].time().path()/".."/"postProcessing/DMD"/snapshots[0].name();
        }
        else
        {
            outputDir = snapshots[0].time().path()/"postProcessing/DMD"/snapshots[0].name();
        }

        if (!isDir(outputDir))
        {
            mkDir(outputDir);
        }

        OFstream ritzValuesOf(outputDir/"ritzValues");

        ritzValuesOf << "# Ritz values obtained from dynamic mode decomposition" << nl
                     << "# " << nl
                     << "# " << nl
                     << "# First entry is the Ritz value's corresponding mode number." << nl
                     << "# Second entry is the Ritz value, e.g." << nl
                     << "# " << nl
                     << "# <mode> <real part> <imaginary part>" << nl 
                     << "#" << nl
                     << endl;

        OFstream amplitudesOf(outputDir/"amplitudes");

	amplitudesOf << "# Amplitudes of modes" << nl
                     << "# " << nl
                     << "# " << nl
                     << "# First entry is the mode number." << nl
                     << "# Second entry is the amplitude, e.g." << nl
                     << "# " << nl
                     << "# <mode> <amplitude>" << nl 
                     << "#" << nl
                     << endl;

        OFstream growthRatesOf(outputDir/"growthRates");

        growthRatesOf << "# Growth rates of modes, defined as log||RitzValues||/deltaT" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Second entry is the growth rate, e.g." << nl
                      << "# " << nl
                      << "# <mode> <growthRate>" << nl 
                      << "#" << nl
                      << endl;

        OFstream frequenciesOf(outputDir/"frequencies");

        frequenciesOf << "# Frequency of modes, defined as arg(RitzValues)/deltaT" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Second entry is the frequency" << nl
                      << "# " << nl
                      << "# <mode> <frequency>" << nl 
                      << "#" << nl
                      << endl;

        OFstream energiesOf(outputDir/"energies");

        energiesOf << "# Energies of modes, defined as amplitude*(exp(2*growthRate*T)-1)/(2*growthRate*T)" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Second entry ist the energy" << nl
                      << "# " << nl
                      << "# <mode> <energy>" << nl 
                      << "#" << nl
                      << endl;

        for (label baseI = 0; baseI < nSnapshots-1; baseI++)
        {
            // Grap mode index
            const label& modeEnergyI = indicesModeEnergies[baseI];

            ritzValuesOf << setw(5) << baseI;
            amplitudesOf << setw(5) << baseI;
            growthRatesOf << setw(5) << baseI;
            frequenciesOf << setw(5) << baseI;
            energiesOf << setw(5) << baseI;


            // Valid component, grab complex eigen values (also called Ritz values)
            const complexField& ritzValues = eigenBase.ritzValues();

            ritzValuesOf << setw(20) << ritzValues[modeEnergyI].Re() << " " 
	                 << setw(20) << ritzValues[modeEnergyI].Im() << " ";   

            // Valid component, grab mode amplitude
            const scalarField& amplitudes = eigenBase.modeNorms();

            amplitudesOf << setw(20) << amplitudes[modeEnergyI] << " "; 

            // Valid component, calculate frequency and growth rate from Ritz value
            std::complex<scalar> tmp;
            tmp.real() = ritzValues[modeEnergyI].Re();
            tmp.imag() = ritzValues[modeEnergyI].Im();
            tmp = std::log(tmp+SMALL)/(deltaT_*2.0*constant::mathematical::pi);
                    
            // Valid component, calculate growth rate from Ritz value
            //growthRatesOf << setw(15) << Foam::log(sqrt(sqr(ritzValues[baseI].Re())+sqr(ritzValues[baseI].Im())))/deltaT_ << " "; 
            growthRatesOf << setw(20) << tmp.real() << " ";   
                    
            // Valid component, calculate frequency from Ritz value
            //frequenciesOf << setw(20) << std::atan2(ritzValues[baseI].Im(), ritzValues[baseI].Re())/(deltaT_*2.0*constant::mathematical::pi) << " "; 
            frequenciesOf << setw(20) << tmp.imag() << " "; 

            // Energy contribution 
            if(modeSelection_ == "energy")
            {
                scalar T = (nSnapshots-1)*deltaT_;
                energiesOf << setw(20) << amplitudes[modeEnergyI]*(Foam::exp(2.0*tmp.real()*T)-1)/(2.0*tmp.real()*T+SMALL) << " ";
            }
            else
            {
                energiesOf << setw(20) <<  amplitudes[modeEnergyI] << " ";
            }

            ritzValuesOf << endl;
            amplitudesOf << endl;
            growthRatesOf << endl;
            frequenciesOf << endl;
            energiesOf << endl;
        }
    }

    // Assemble Ritz values to compact field for reconstruction 

    ritzValues_.setSize(1);
    const complexField& ritzValues = eigenBase.ritzValues();

    complexField* sortedRitzValues = new complexField(baseSize_);

    forAll(*sortedRitzValues, baseI)
    {
        (*sortedRitzValues)[baseI] = ritzValues[indicesModeEnergies[baseI]]; 
    }

    // Valid component, grab complex eigen values (also called Ritz values)
    ritzValues_.set
    (
        0,
        sortedRitzValues
    );

}


// ************************************************************************* //
