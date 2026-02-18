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

  const volScalarField& e1(mixture.thermo1().he());
  const volScalarField& e2(mixture.thermo2().he());

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture_.T();

  // Cv(T,p,Y...) from rhoFluidMulticomponentThermo -> variable field
  const volScalarField Cp1(mixture.thermo1().Cp());
  const volScalarField Cp2(mixture.thermo2().Cp());

  // Conservative coefficients (alpha*rho*Cv)
  const volScalarField aRhoCp1("aRhoCp1", alpha1*rho1*Cp1);
  const volScalarField aRhoCp2("aRhoCp2", alpha2*rho2*Cp2);

  const Foam::upwind<scalar> upCp1(mesh, alphaRhoPhi1);
  const Foam::upwind<scalar> upCp2(mesh, alphaRhoPhi2);

  const surfaceScalarField Cp1f("Cp1f", upCp1.interpolate(Cp1));
  const surfaceScalarField Cp2f("Cp2f", upCp2.interpolate(Cp2));

  const surfaceScalarField aRhoCpPhi1("aRhoCpPhi1", Cp1f*alphaRhoPhi1);
  const surfaceScalarField aRhoCpPhi2("aRhoCpPhi2", Cp2f*alphaRhoPhi2);

  // Product-rule correction:
  //  - T * [ ddt(alpha*rho*Cv) + div(alpha*rho*phi*Cv) ]
  // Implemented as -Sp( ddt(...) + div(...), T )
  const volScalarField prodCorr1("prodCorr1", fvc::ddt(aRhoCp1) + fvc::div(aRhoCpPhi1));
  const volScalarField prodCorr2("prodCorr2", fvc::ddt(aRhoCp2) + fvc::div(aRhoCpPhi2));

  fvScalarMatrix TEqn
  (
      correction
      (
          (
              fvm::ddt(aRhoCp1, T)
            + fvm::div(aRhoCpPhi1, T)
            - fvm::Sp(prodCorr1, T)              // <-- key fix
            - (
                  e1Source.hasDiag()
                ? fvm::Sp(contErr1()*Cp1, T) + fvm::Sp(e1Source.A()*Cp1, T)
                : fvm::Sp(contErr1()*Cp1, T)
              )
          )
        + (
              fvm::ddt(aRhoCp2, T)
            + fvm::div(aRhoCpPhi2, T)
            - fvm::Sp(prodCorr2, T)              // <-- key fix
            - (
                  e2Source.hasDiag()
                ? fvm::Sp(contErr2()*Cp2, T) + fvm::Sp(e2Source.A()*Cp2, T)
                : fvm::Sp(contErr2()*Cp2, T)
              )
          )
      )

    + fvc::ddt(alpha1, rho1, e1) + fvc::div(alphaRhoPhi1, e1)
    - contErr1()*e1
    + fvc::ddt(alpha2, rho2, e2) + fvc::div(alphaRhoPhi2, e2)
    - contErr2()*e2

    - fvm::laplacian(thermophysicalTransport.kappaEff(), T)

    + (
          mixture.totalInternalEnergy()
        ?
          fvc::div(fvc::absolute(phi, U), p)()()
        + (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))()()
        - (U()&(fvModels().source(rho, U)&U)()) - (contErr1() + contErr2())*K
        :
          p*fvc::div(fvc::absolute(phi, U))()()
      )
    ==
      (e1Source&e1)
    + (e2Source&e2)
  );

  TEqn.relax();

  fvConstraints().constrain(TEqn);

  TEqn.solve();

  fvConstraints().constrain(T);

  mixture_.correctThermo();
  mixture_.correct();
}


// ************************************************************************* //
