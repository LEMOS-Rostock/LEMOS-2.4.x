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

#include "mosDMDOrthoNormalBase.H"
#include "mosDMD.H"
#include "mosDMDEigenBase.H"
#include "IOmanip.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::mosDMDOrthoNormalBase<Type>::calcOrthoBase
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots
)
{
    // Calculate eigen base for each component
    PtrList<mosDMDEigenBase> eigenBase(pTraits<Type>::nComponents);

    const label nSnapshots = snapshots.size();

    typename
    powProduct<Vector<label>, pTraits<Type>::rank>::type validComponents
    (
        pow
        (
            snapshots[0].mesh().solutionD(),
            pTraits
            <
                typename powProduct<Vector<label>,
                pTraits<Type>::rank
            >::type>::zero
        )
    );

    label nValidCmpts = 0;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        // Component not valid, skipping
        if (validComponents[cmpt] == -1) continue;

        // Create a list of snapshots
        PtrList<volScalarField> sf(nSnapshots);

        for (label i=0; i< nSnapshots; i++)
        {
            sf.set(i, new volScalarField(snapshots[i].component(cmpt)));
        }

        // Create eigen base
        eigenBase.set(cmpt, new mosDMDEigenBase(sf, nSnapshots));
   
        nValidCmpts++;
    }

    eigenBase.setSize(nValidCmpts);

    Info << "Number of valid eigen components: " << nValidCmpts << endl;

    // Compute selected base size or entire otherwise
    if ((baseSize_ <= 0) || (baseSize_ > nSnapshots-1))
    {
        baseSize_ = nSnapshots-1;
    }

    Info << "Base size: " << baseSize_ << endl;

    // Compute mode energies
    SortableList<scalar> sortedList(nSnapshots-1,0.0);
    scalar T = (nSnapshots-1)*deltaT_;
    
    for (label baseI = 0; baseI < nSnapshots-1; baseI++)
    {
        nValidCmpts = 0;
        scalar energyCmpts = 0.0;

        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            if (validComponents[cmpt] != -1)
            {
                const scalarField& amplitudes = eigenBase[nValidCmpts].modeNorms();

                if(modeSelection_ == "energy")
                {
                    const complexField& ritzValues = eigenBase[nValidCmpts].ritzValues();

                    std::complex<scalar> tmp;
                    tmp.real() = ritzValues[baseI].Re();
                    tmp.imag() = ritzValues[baseI].Im();
                
                    tmp = std::log(tmp+SMALL)/(deltaT_*2.0*constant::mathematical::pi);

                    energyCmpts += sqr(amplitudes[baseI]*(Foam::exp(2.0*tmp.real()*T)-1)/(2.0*tmp.real()*T+SMALL));
                }
                else
                {
                    energyCmpts += sqr(amplitudes[baseI]);
                }
                
                nValidCmpts++;
            }

        }

        sortedList[baseI] = std::sqrt(energyCmpts);
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
        GeometricField<Type, fvPatchField, volMesh>* onBaseRePtr
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + name(baseI) + "Real",
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensioned<Type>
                (
                    "zero",
                    snapshots[0].dimensions(),
                    pTraits<Type>::zero
                )
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& onBaseRe = *onBaseRePtr;
        
        GeometricField<Type, fvPatchField, volMesh>* onBaseImPtr
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    snapshots[0].name() + "DMD" + name(baseI) + "Imag",
                    snapshots[0].time().timeName(),
                    snapshots[0].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                snapshots[0].mesh(),
                dimensioned<Type>
                (
                    "zero",
                    snapshots[0].dimensions(),
                    pTraits<Type>::zero
                )
            )
        );
        GeometricField<Type, fvPatchField, volMesh>& onBaseIm = *onBaseImPtr;

        nValidCmpts = 0;

        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            if (validComponents[cmpt] != -1)
            {
                // Valid component, grab complex build coeffs
                const complexField& coeffs =
                    eigenBase[nValidCmpts].coeffs()[indicesModeEnergies[baseI]];

                nValidCmpts++;

                volScalarField onBaseReCmpt = onBaseRe.component(cmpt);
                volScalarField onBaseImCmpt = onBaseIm.component(cmpt);

                forAll (coeffs, coeffI)
                {
                    onBaseReCmpt +=
                        coeffs[coeffI].Re()*snapshots[coeffI].component(cmpt);
                        
                    onBaseImCmpt +=
                        coeffs[coeffI].Im()*snapshots[coeffI].component(cmpt);
                }


                onBaseRe.replace(cmpt, onBaseReCmpt);
                onBaseIm.replace(cmpt, onBaseImCmpt);

            }
            else
            {
                // Component invalid.  Grab first snapshot.
                onBaseRe.replace
                (
                    cmpt,
                    snapshots[0].component(cmpt)
                );
                
                // Component invalid.  Grab first snapshot.
                onBaseIm.replace
                (
                    cmpt,
                    snapshots[0].component(cmpt)
                );
            }
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
                     << "# Subsequent entries are the ritz values in componentwise order, e.g." << nl
                     << "# " << nl
                     << "# <mode> <real part cmpt_1> <imaginary part cmpt_1>   <real part cmpt_2> <imaginary part cmpt_2>   ..." << nl 
                     << "#" << nl
                     << endl;

        OFstream amplitudesOf(outputDir/"amplitudes");

	amplitudesOf << "# Amplitudes of modes" << nl
                     << "# " << nl
                     << "# " << nl
                     << "# First entry is the mode number." << nl
                     << "# Subsequent entries are the amplitudes of the mode components" << nl
                     << "# " << nl
                     << "# <mode> <amplitude cmpt_1>   <amplitude cmpt_2>   <amplitude cmpt_3>   ..." << nl 
                     << "#" << nl
                     << endl;

        OFstream growthRatesOf(outputDir/"growthRates");

        growthRatesOf << "# Growth rates of modes, defined as log||RitzValues||/deltaT" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Subsequent entries are the growth rates of the mode components" << nl
                      << "# " << nl
                      << "# <mode> <growthRate cmpt_1>   <growthRate cmpt_2>   <growthRate cmpt_3>   ..." << nl 
                      << "#" << nl
                      << endl;

        OFstream frequenciesOf(outputDir/"frequencies");

        frequenciesOf << "# Frequency of modes, defined as arg(RitzValues)/deltaT" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Subsequent entries are the frequencies of the mode components" << nl
                      << "# " << nl
                      << "# <mode> <frequency cmpt_1>   <frequency cmpt_2>   <frequency cmpt_3>   ..." << nl 
                      << "#" << nl
                      << endl;

        OFstream energiesOf(outputDir/"energies");

        energiesOf << "# Energies of modes, defined as amplitude*(exp(2*growthRate*T)-1)/(2*growthRate*T)" << nl
                      << "# " << nl
                      << "# " << nl
                      << "# First entry is the mode number." << nl
                      << "# Subsequent entries are the energies of the mode components" << nl
                      << "# " << nl
                      << "# <mode> <energy cmpt_1>   <energy cmpt_2>   <energy cmpt_3>   ..." << nl 
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

            nValidCmpts = 0;
            
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                if (validComponents[cmpt] != -1)
                {

                    // Valid component, grab complex eigen values (also called Ritz values)
                    const complexField& ritzValues = eigenBase[nValidCmpts].ritzValues();

                    ritzValuesOf << setw(20) << ritzValues[modeEnergyI].Re() << " " 
		                 << setw(20) << ritzValues[modeEnergyI].Im() << " ";   

                    // Valid component, grab mode amplitude
                    const scalarField& amplitudes = eigenBase[nValidCmpts].modeNorms();

                    amplitudesOf << setw(20) << amplitudes[modeEnergyI] << " "; 

                    // Valid component, calculate growth rate from Ritz value

                    // Valid component, calculate frequency and growth rate from Ritz value
                    std::complex<scalar> tmp;
                    tmp.real() = ritzValues[modeEnergyI].Re();
                    tmp.imag() = ritzValues[modeEnergyI].Im();
                    tmp = std::log(tmp+SMALL)/(deltaT_*2.0*constant::mathematical::pi);
                    
                    // Valid component, calculate growth rate from Ritz value
                    //growthRatesOf << setw(15) << Foam::log(sqrt(sqr(ritzValues[baseI].Re())+sqr(ritzValues[baseI].Im())))/deltaT_ << " "; 
                    growthRatesOf << setw(20) << tmp.real() << " ";   
                    
                    // Valid component, calculate frequency from Ritz value
                    //frequenciesOf << setw(15) << atan2(ritzValues[baseI].Im(), ritzValues[baseI].Re())/deltaT_ << " "; 
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
                
                    nValidCmpts++;
                }
                else
                {
                    ritzValuesOf << setw(5) << snapshots[0].component(cmpt)()[0] << " " 
                                 << setw(5) << snapshots[0].component(cmpt)()[0] << " ";
                    
                    amplitudesOf << setw(5) << snapshots[0].component(cmpt)()[0] << " "; 

                    growthRatesOf << setw(5) << snapshots[0].component(cmpt)()[0] << " "; 

                    frequenciesOf << setw(5) << snapshots[0].component(cmpt)()[0] << " "; 

                    energiesOf << setw(5) << snapshots[0].component(cmpt)()[0] << " "; 
                }
              

            }

            ritzValuesOf << endl;
            amplitudesOf << endl;
            growthRatesOf << endl;
            frequenciesOf << endl;
            energiesOf << endl;
        }
    }

    // Assemble Ritz values to compact field for reconstruction 
    ritzValues_.setSize(pTraits<Type>::nComponents);
    nValidCmpts = 0;

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {

        if (validComponents[cmpt] != -1)
        {
            const complexField& ritzValues = eigenBase[nValidCmpts].ritzValues();

            complexField* sortedRitzValues = new complexField(baseSize_);

            forAll(*sortedRitzValues, baseI)
            {
                 (*sortedRitzValues)[baseI] = ritzValues[indicesModeEnergies[baseI]]; 
            }

            // Valid component, grab complex eigen values (also called Ritz values)
            ritzValues_.set
            (
                cmpt,
                sortedRitzValues
            );

            nValidCmpts++;
        }
        else
        {
            ritzValues_.set
            (
                cmpt,
                new complexField(baseSize_, complex(snapshots[0].component(cmpt)()[0], snapshots[0].component(cmpt)()[0]))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mosDMDOrthoNormalBase<Type>::mosDMDOrthoNormalBase
(
    const fvMesh& mesh,
    const word& DMDModelName
)
:
    DMDModel<Type>(mesh, DMDModelName),
    baseSize_(readLabel((*this).coeffDict_.lookup("nModes"))),
    deltaT_(readScalar((*this).coeffDict_.lookup("deltaT"))),
    modeSelection_((*this).coeffDict_.lookup("modeSelection"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::mosDMDOrthoNormalBase<Type>::~mosDMDOrthoNormalBase()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
