/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "compressibleVoFC.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "upwind.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoFC::thermophysicalPredictor()
{
  compositionPredictor();

    const volScalarField& rho1(mixture.rho1());
    const volScalarField& rho2(mixture.rho2());

    const volScalarField& Cpv1(mixture.thermo1().Cpv());
    const volScalarField& Cpv2(mixture.thermo2().Cpv());


    volScalarField& e1 = mixture_.thermo1().he();
    volScalarField& e2 = mixture_.thermo2().he();


    const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
    const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));
        volScalarField& T = mixture_.T();

   // const volScalarField rhoCpv1("rhoCpv1", rho1*Cpv1);
   // const volScalarField rhoCpv2("rhoCpv2", rho2*Cpv2);
    const volScalarField rhoCpv("rhoCpv", alpha1*rho1*Cpv1+alpha2*rho2*Cpv2);

    

    
    const surfaceScalarField alphaRhoCpvPhi1
    (
        "alphaRhoCpvPhi1",
        fvc::interpolate(Cpv1)*alphaRhoPhi1
    );

    const surfaceScalarField alphaRhoCpvPhi2
    (
        "alphaRhoCpvPhi2",
        fvc::interpolate(Cpv2)*alphaRhoPhi2
    );

    const surfaceScalarField rhoCpvPhi
    (   
        "rhoCpvPhi",
        alphaRhoCpvPhi1 + alphaRhoCpvPhi2
    );

    fvScalarMatrix TEqn
    (
      fvm::ddt(rhoCpv, T)
    + fvm::div(rhoCpvPhi, T)
    - fvm::Sp(fvc::ddt(rhoCpv) + fvc::div(rhoCpvPhi), T)
    - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
    ==
      (e1Source&e1)
    + (e2Source&e2)
    );


    TEqn.relax();

    fvConstraints().constrain(TEqn);

    TEqn.solve();

    fvConstraints().constrain(T);
    // thermo update

    mixture_.correctThermo();
    mixture_.correct();
}


// ************************************************************************* //
