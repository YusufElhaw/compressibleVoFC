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
#include "fvcLaplacian.H"
#include "upwind.H"

namespace
{
  template<class FieldType>
void writeDebugField
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime,
    const Foam::word& name,
    const FieldType& field
)
{
    FieldType debugField
    (
        Foam::IOobject
        (
            name,
            runTime.name(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        field
    );

    debugField.write();
}


template<class FieldType>
void writeDebugField
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime,
    const Foam::word& name,
    const Foam::tmp<FieldType>& tfield
)
{
    const FieldType& field = tfield();

    FieldType debugField
    (
        Foam::IOobject
        (
            name,
            runTime.name(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        field
    );

    debugField.write();

}
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
/*
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
}*/
void Foam::solvers::compressibleVoFC::thermophysicalPredictor()
{ 
  compositionPredictor();
  //T_0  298.15; // reference temperature of 298.15 K   $FOAM_ETC/controlDic
  Info <<nl<< "thermophysicalPredictor läuft"<<endl;
  
  const volScalarField& rho1(mixture.thermo1().rho()); 
  const volScalarField& rho2(mixture.thermo2().rho());
  //const volScalarField& rho(mixture.rho());

  const volScalarField& Cpv1(mixture.thermo1().Cpv());
  const volScalarField& Cpv2(mixture.thermo2().Cpv());
  const volScalarField& Cpv((alpha1*rho1*Cpv1+alpha2*rho2*Cpv2)/rho);

  // Energy of the phases (in physicalProperties Choose Enthalpy or InternalEnergy)
  const volScalarField& e1(mixture.thermo1().he());   //he Enthalpy or InternalEnergy
  const volScalarField& e2(mixture.thermo2().he());   //he Enthalpy or InternalEnergy

  const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
  const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

  volScalarField& T = mixture_.T();
const dimensionedScalar T0("T0", dimTemperature, 298.15);


    const volScalarField limAlpha1( min(max(alpha1, scalar(0)), scalar(1)) );
    const wordList TpatchTypes = T.boundaryField().types();
    surfaceScalarField alphaEffRho
    (
        IOobject
        (
            "alphaEffRho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho*thermophysicalTransport.kappaEff()/(Cpv*(rho1*alpha1+rho2*alpha2)))
      
    );
//Info<<"Cpv1 beträgt"<<Cpv1<<endl;

  volScalarField H
  (
      IOobject
      (
        "H",
        mesh.time().name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("H", dimEnergy/dimMass, 0.0),
      T.boundaryField().types()
  );

  H = ( (T - T0)*(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2) )/rho;

  H.correctBoundaryConditions();
  Info<<"thermophysical 1"<<endl;
  const scalar RelaxFac = 1.0;
  fvScalarMatrix EEqn
    (
        fvm::ddt(rho, H)
      + fvm::div(rhoPhi, H)
      ==
        fvc::laplacian(thermophysicalTransport.kappaEff(), T)
      + RelaxFac*
      (   
        fvm::laplacian(alphaEffRho, H) 
      - fvc::laplacian(alphaEffRho, H) 
      ) 
            
    );
  Info<<"thermophysical 2"<<endl;

      EEqn.relax();
      fvConstraints().constrain(EEqn);
      EEqn.solve();
  Info<<"thermophysical 3"<<endl;
  if (runTime.writeTime())
    {
  
  writeDebugField(mesh, runTime, "alphaEffRho_debug", alphaEffRho);
  writeDebugField(mesh, runTime, "fvm_laplacian_alphaEffRho_H_debug", fvc::laplacian(alphaEffRho, H));
    }
  
  
      //Now reevaluate T for the updated enthalpy fields
      T = T0 + rho*H/(limAlpha1*rho1*Cpv1 + (1.0 - limAlpha1)*rho2*Cpv2);
      mixture_.correctThermo();
      mixture_.correct();
}        

// ************************************************************************* //
