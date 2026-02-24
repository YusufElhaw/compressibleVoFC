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

   // const volScalarField& Cpv1(mixture.thermo1().Cpv());
   // const volScalarField& Cpv2(mixture.thermo2().Cpv());


    volScalarField& he1 = mixture_.thermo1().he();
    volScalarField& he2 = mixture_.thermo2().he();


    const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, he1));
    const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, he2));
Info<< "Defintion of EEqn1"<<endl;

    fvScalarMatrix  EEqn1
    (
        fvm::ddt(alpha1, rho1, he1)
      + fvm::div(alphaRhoPhi1, he1)
      - fvm::Sp(contErr1, he1)

      + fvc::ddt(alpha1, rho1, K) 
      + fvc::div(alphaRhoPhi1, K)
      //- contErr1*K

      - thermophysicalTransport.divq(he1) 
     ==
    e1Source
    );
Info<< "Defintion of EEqn2"<<endl;
    
    fvScalarMatrix  EEqn2
    (
        fvm::ddt(alpha2, rho2, he2)
      + fvm::div(alphaRhoPhi2, he2)
      - fvm::Sp(contErr2, he2)

      + fvc::ddt(alpha2, rho2, K) 
      + fvc::div(alphaRhoPhi2, K)
      //- contErr2*K

      - thermophysicalTransport.divq(he2) 
     ==
    e2Source
    );
    
    //  Relax/Constraints/Solve

    EEqn1.relax();
    fvConstraints().constrain(EEqn1);
    EEqn1.solve();
    fvConstraints().constrain(he1);
    Info<< "solved EEqn1"<<endl;

    EEqn2.relax();
    Info<< "realxed EEqn2"<<endl;
    fvConstraints().constrain(EEqn2);
    Info<< "fvConstraints().constrain EEqn2"<<endl;
    EEqn2.solve();
    Info<< "solved EEqn2"<<endl;
    fvConstraints().constrain(he2);
    Info<< "fvConstraints().constrain(he2); EEqn2"<<endl;

    // thermo update

    mixture_.correctThermo();
    mixture_.correct();
}


// ************************************************************************* //
