/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

Class
    Foam::compressibleTwoPhaseVoFMixtureC

Description
    Class to represent a mixture of two rhoFluidMulticomponentThermo-based phases

SourceFiles
    compressibleTwoPhaseVoFMixtureC.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhaseVoFMixtureC_H
#define compressibleTwoPhaseVoFMixtureC_H

#include "twoPhaseVoFMixture.H"
#include "compressibleTwoPhases.H"
#include "rhoFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressibleTwoPhaseVoFMixtureC Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhaseVoFMixtureC
:
    public twoPhaseVoFMixture,
    virtual public compressibleTwoPhases,
    public viscosity
{
    // Private Data

        //- Switch to choose between solving for internal energy
        //  or total internal energy which is the default
        Switch totalInternalEnergy_;

        //- Pressure
        volScalarField p_;

        //- Mixture temperature
        volScalarField T_;


        // --- Species definition
        word species1Name_;
        word species2Name_;

        //- Species mass-fractions per phase

        autoPtr<volScalarField> s1Phase1_;
        autoPtr<volScalarField> s2Phase1_;
        autoPtr<volScalarField> s1Phase2_;
        autoPtr<volScalarField> s2Phase2_;

        //- Indices in thermo species tables (per phase thermo)
        label s1IndexThermo1_;
        label s2IndexThermo1_;
        label s1IndexThermo2_;
        label s2IndexThermo2_;


        //- Thermo-package of phase 1
        autoPtr<rhoFluidMulticomponentThermo> thermo1_;

        //- Thermo-package of phase 2
        autoPtr<rhoFluidMulticomponentThermo> thermo2_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture viscosity
        volScalarField nu_;


public:

    //- Runtime type information
    TypeName("compressibleTwoPhaseVoFMixtureC");


    // Constructors

        //- Construct from a mesh
        compressibleTwoPhaseVoFMixtureC(const fvMesh& mesh);


    //- Destructor
    virtual ~compressibleTwoPhaseVoFMixtureC();


    // Member Functions

        //- Return true to solve for total internal energy
        //  return false to solve for internal energy
        bool totalInternalEnergy() const
        {
            return totalInternalEnergy_;
        }

        //- Return pressure [Pa]
        volScalarField& p()
        {
            return p_;
        }

        //- Return mixture temperature [K]
        volScalarField& T()
        {
            return T_;
        }

        //- Return mixture temperature [K]
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return the thermo for phase 1
        const rhoFluidMulticomponentThermo& thermo1() const
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        const rhoFluidMulticomponentThermo& thermo2() const
        {
            return thermo2_();
        }

        //- Return the thermo for phase 1
        rhoFluidMulticomponentThermo& thermo1()
        {
            return thermo1_();
        }

        //- Return the thermo for phase 2
        rhoFluidMulticomponentThermo& thermo2()
        {
            return thermo2_();
        }

        //- Return the density of phase 1
        const volScalarField& rho1() const
        {
            return thermo1_->rho();
        }

        //- Return the density of phase 2
        const volScalarField& rho2() const
        {
            return thermo2_->rho();
        }

            /////
            // Species accessors
            volScalarField& s1Phase1() { return s1Phase1_(); }
            volScalarField& s2Phase1() { return s2Phase1_(); }
            volScalarField& s1Phase2() { return s1Phase2_(); }
            volScalarField& s2Phase2() { return s2Phase2_(); }

            const volScalarField& s1Phase1() const { return s1Phase1_(); }
            const volScalarField& s2Phase1() const { return s2Phase1_(); }
            const volScalarField& s1Phase2() const { return s1Phase2_(); }
            const volScalarField& s2Phase2() const { return s2Phase2_(); }

            const word& species1Name() const { return species1Name_; }
            const word& species2Name() const { return species2Name_; }
            //////

        //- The fluid is incompressible if both phases are incompressible
        bool incompressible() const
        {
            return
                thermo1_->incompressible()
             && thermo2_->incompressible();
        }

        //- Return mixture density [kg/m^3]
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the mixture compressibility/density
        virtual tmp<volScalarField> psiByRho() const;

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();
        
        //- Correct the Composition of each phase
        virtual void correctComposition();

        //- Update mixture properties
        virtual void correct();

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }


    // IO

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
