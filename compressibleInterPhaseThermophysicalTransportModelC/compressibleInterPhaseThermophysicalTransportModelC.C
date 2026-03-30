/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "compressibleInterPhaseThermophysicalTransportModelC.H"
#include "compressibleInterPhaseTransportModelC.H"
#include "fvcSnGrad.H"
#include "compressibleTwoPhaseVoFMixtureC.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleInterPhaseThermophysicalTransportModelC::
compressibleInterPhaseThermophysicalTransportModelC
(
    const compressibleInterPhaseTransportModelC& momentumTransport
)
:
    thermophysicalTransportModel(momentumTransport.mixture_.mesh(), word::null),
    momentumTransport_(momentumTransport),
    massDiffusivity_(momentumTransport) 
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::compressibleInterPhaseThermophysicalTransportModelC::read()
{
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::kappaEff() const
{
    const compressibleTwoPhaseVoFMixtureC& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.momentumTransport1_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.momentumTransport2_->nut()
            );
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            );
    }
}


Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::kappaEff
(
    const label patchi
) const
{
    const compressibleTwoPhaseVoFMixtureC& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *momentumTransport_.momentumTransport1_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *momentumTransport_.momentumTransport2_->nut(patchi)
            );
    }
    else
    {
        return
            mixture_.alpha1().boundaryField()[patchi]
           *(
                mixture_.thermo1().kappa().boundaryField()[patchi]
              + mixture_.thermo1().rho(patchi)
               *mixture_.thermo1().Cp().boundaryField()[patchi]
               *momentumTransport_.mixtureMomentumTransport_->nut(patchi)
            )
          + mixture_.alpha2().boundaryField()[patchi]
           *(
                mixture_.thermo2().kappa().boundaryField()[patchi]
              + mixture_.thermo2().rho(patchi)
               *mixture_.thermo2().Cp().boundaryField()[patchi]
               *momentumTransport_.mixtureMomentumTransport_->nut(patchi)
            );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::alphaEff() const
{
    const compressibleTwoPhaseVoFMixtureC& mixture_ =
        momentumTransport_.mixture_;

    if (momentumTransport_.twoPhaseTransport_)
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.momentumTransport1_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.momentumTransport2_->nut()
            )/mixture_.thermo2().Cv();
    }
    else
    {
        return
            mixture_.alpha1()
           *(
                mixture_.thermo1().kappa()
              + mixture_.thermo1().rho()*mixture_.thermo1().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )/mixture_.thermo1().Cv()
          + mixture_.alpha2()
           *(
                mixture_.thermo2().kappa()
              + mixture_.thermo2().rho()*mixture_.thermo2().Cp()
               *momentumTransport_.mixtureMomentumTransport_->nut()
            )/mixture_.thermo2().Cv();
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::q() const
{
    return surfaceScalarField::New
    (
        "q",
        -fvc::interpolate(kappaEff())
        *fvc::snGrad(momentumTransport_.mixture_.T())
    );
}


Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::q
(
    const label patchi
) const
{
    return
      - kappaEff(patchi)
       *momentumTransport_.mixture_.T().boundaryField()[patchi].snGrad();
}


Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::qCorr
(
    const label patchi
) const
{
    return tmp<scalarField>(nullptr);
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::compressibleInterPhaseThermophysicalTransportModelC::divq
(
    volScalarField& he
) const
{
//    NotImplemented;
//
//    return tmp<fvScalarMatrix>(nullptr);
    const compressibleTwoPhaseVoFMixtureC& mixture = momentumTransport_.mixture_;

    // decide which phase this 'he' belongs to by name
    // (matches your fields: he.phase1 / he.phase2)
    const bool isPhase1 = (he.name().find("phase1") != string::npos);

    tmp<volScalarField> talphahe;

    if (momentumTransport_.twoPhaseTransport_)
    {
        if (isPhase1)
        {
            const rhoFluidMulticomponentThermo& thermo = mixture.thermo1();

            talphahe =
                mixture.alpha1()
               *(
                    thermo.kappa()
                  + thermo.rho()*thermo.Cp()*momentumTransport_.momentumTransport1_->nut()
                )
               /thermo.Cpv();
        }
        else
        {
            const rhoFluidMulticomponentThermo& thermo = mixture.thermo2();

            talphahe =
                mixture.alpha2()
               *(
                    thermo.kappa()
                  + thermo.rho()*thermo.Cp()*momentumTransport_.momentumTransport2_->nut()
                )
               /thermo.Cpv();
        }
    }
    else
    {
        // single mixture turbulence model (same nut for both)
        if (isPhase1)
        {
            const rhoFluidMulticomponentThermo& thermo = mixture.thermo1();

            talphahe =
                mixture.alpha1()
               *(
                    thermo.kappa()
                  + thermo.rho()*thermo.Cp()*momentumTransport_.mixtureMomentumTransport_->nut()
                )
               /thermo.Cpv();
        }
        else
        {
            const rhoFluidMulticomponentThermo& thermo = mixture.thermo2();

            talphahe =
                mixture.alpha2()
               *(
                    thermo.kappa()
                  + thermo.rho()*thermo.Cp()*momentumTransport_.mixtureMomentumTransport_->nut()
                )
               /thermo.Cpv();
        }
    }

    // unityLewisFourier-style: divq(he) = -laplacian(alphahe, he)
    // (here alphahe already includes the phase fraction alpha_i)
    return -fvm::laplacian(talphahe(), he);



}
        // * * * * * * * * * * * Mass Transfer * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::DEff() const
{
    return massDiffusivity_.DEff();
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::D1Eff() const
{
    return massDiffusivity_.D1Eff();
}

Foam::tmp<Foam::volScalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::D2Eff() const
{
    return massDiffusivity_.D2Eff();
}

Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::DEff
(
    const label patchi
) const
{
   return massDiffusivity_.DEff(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::D1Eff
(
    const label patchi
) const
{
   return massDiffusivity_.D1Eff(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::compressibleInterPhaseThermophysicalTransportModelC::D2Eff
(
    const label patchi
) const
{
   return massDiffusivity_.D2Eff(patchi);
}
void Foam::compressibleInterPhaseThermophysicalTransportModelC::predict()
{}


void Foam::compressibleInterPhaseThermophysicalTransportModelC::correct()
{}


// ************************************************************************* //
