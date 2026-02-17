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

/*
  const volScalarField& rho1(mixture.rho1());
  const volScalarField& rho2(mixture.rho2());

  const volScalarField& e1(mixture.thermo1().he());
  const volScalarField& e2(mixture.thermo2().he());

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture_.T();

  // Cv(T,p,Y...) from rhoFluidMulticomponentThermo -> variable field
  const volScalarField Cv1(mixture.thermo1().Cv());
  const volScalarField Cv2(mixture.thermo2().Cv());

  // Conservative coefficients (alpha*rho*Cv)
  const volScalarField aRhoCv1("aRhoCv1", alpha1*rho1*Cv1);
  const volScalarField aRhoCv2("aRhoCv2", alpha2*rho2*Cv2);

  const Foam::upwind<scalar> upCv1(mesh, alphaRhoPhi1);
  const Foam::upwind<scalar> upCv2(mesh, alphaRhoPhi2);

  const surfaceScalarField Cv1f("Cv1f", upCv1.interpolate(Cv1));
  const surfaceScalarField Cv2f("Cv2f", upCv2.interpolate(Cv2));

  const surfaceScalarField aRhoCvPhi1("aRhoCvPhi1", Cv1f*alphaRhoPhi1);
  const surfaceScalarField aRhoCvPhi2("aRhoCvPhi2", Cv2f*alphaRhoPhi2);

  // Product-rule correction:
  //  - T * [ ddt(alpha*rho*Cv) + div(alpha*rho*phi*Cv) ]
  // Implemented as -Sp( ddt(...) + div(...), T )
  const volScalarField prodCorr1("prodCorr1", fvc::ddt(aRhoCv1) + fvc::div(aRhoCvPhi1));
  const volScalarField prodCorr2("prodCorr2", fvc::ddt(aRhoCv2) + fvc::div(aRhoCvPhi2));

  fvScalarMatrix TEqn
  (
      correction
      (
          (
              fvm::ddt(aRhoCv1, T)
            + fvm::div(aRhoCvPhi1, T)
            - fvm::Sp(prodCorr1, T)              // <-- key fix
            - (
                  e1Source.hasDiag()
                ? fvm::Sp(contErr1()*Cv1, T) + fvm::Sp(e1Source.A()*Cv1, T)
                : fvm::Sp(contErr1()*Cv1, T)
              )
          )
        + (
              fvm::ddt(aRhoCv2, T)
            + fvm::div(aRhoCvPhi2, T)
            - fvm::Sp(prodCorr2, T)              // <-- key fix
            - (
                  e2Source.hasDiag()
                ? fvm::Sp(contErr2()*Cv2, T) + fvm::Sp(e2Source.A()*Cv2, T)
                : fvm::Sp(contErr2()*Cv2, T)
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
*/

  // original Equation of OpenFoam // good for pureMixture // bad for multicomponentMixture
  const volScalarField& rho1(mixture.rho1());
  const volScalarField& rho2(mixture.rho2());

  const volScalarField& e1(mixture.thermo1().he());
  const volScalarField& e2(mixture.thermo2().he());

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture_.T();

  fvScalarMatrix TEqn
  (
      correction
      (
          mixture.thermo1().Cp()()
          *(
              fvm::ddt(alpha1, rho1, T) + fvm::div(alphaRhoPhi1, T)
            - (
                  e1Source.hasDiag()
                ? fvm::Sp(contErr1(), T) + fvm::Sp(e1Source.A(), T)
                : fvm::Sp(contErr1(), T)
              )
          )
        + mixture.thermo2().Cp()()
          *(
              fvm::ddt(alpha2, rho2, T) + fvm::div(alphaRhoPhi2, T)
            - (
                  e2Source.hasDiag()
                ? fvm::Sp(contErr2(), T) + fvm::Sp(e2Source.A(), T)
                : fvm::Sp(contErr2(), T)
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
//
  TEqn.relax();

  fvConstraints().constrain(TEqn);

  TEqn.solve();

  fvConstraints().constrain(T);

  mixture_.correctThermo();
  mixture_.correct();
  compositionPredictor();
}


// ************************************************************************* //
