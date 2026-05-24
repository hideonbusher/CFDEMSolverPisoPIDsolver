    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//===============================
// particle interaction modelling
//===============================

    Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 1.0)
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "\nCreating density field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
        //dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), 1.0)
    );

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//===============================

//#   include "createPhi.H"
#ifndef createPhi_H
#define createPhi_H
Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
 (
     IOobject
     (
        "phi",
        runTime.timeName(),
         mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
     ),
     linearInterpolate(U*voidfraction) & mesh.Sf()
 );
#endif

Info<< "Generating interstitial face flux field phiByVoidfraction\n" << endl;
surfaceScalarField phiByVoidfraction
(
 IOobject
 (
    "phiByVoidfraction",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
 ),
 linearInterpolate(U) & mesh.Sf()
);

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );



#include "IOdictionary.H"

// 读取 controlDict 文件
IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,  // 自动检测 controlDict 变化
        IOobject::NO_WRITE
    )
);

// 读取 gravity 字典
dictionary gravityDict = controlDict.subDict("gravity");
scalar gravity_x = gravityDict.lookupOrDefault<scalar>("x", 0.0);
scalar gravity_y = gravityDict.lookupOrDefault<scalar>("y", 0.0);
scalar gravity_z = gravityDict.lookupOrDefault<scalar>("z", -9.81);
vector gravityValue(gravity_x, gravity_y, gravity_z);

// 读取 drive 字典
dictionary driveDict = controlDict.subDict("drive");
scalar drive_x = driveDict.lookupOrDefault<scalar>("x", 0.0);
scalar drive_y = driveDict.lookupOrDefault<scalar>("y", 0.0);
scalar drive_z = driveDict.lookupOrDefault<scalar>("z", 0.0);
vector driveValue(drive_x, drive_y, drive_z);

// 读取 flowControl 字典tiqu liuliang
dictionary flowControlDict = controlDict.subDict("flowControl");
scalar Q_target = flowControlDict.lookupOrDefault<scalar>("Q_target", 0.0);
scalar Kp = flowControlDict.lookupOrDefault<scalar>("Kp", 1.0);
scalar Ki = flowControlDict.lookupOrDefault<scalar>("Ki", 0.01);
scalar Kd = flowControlDict.lookupOrDefault<scalar>("Kd", 0.0);  // 微分增益
scalar max_integral = flowControlDict.lookupOrDefault<scalar>("max_integral", 1);
scalar x_plane = flowControlDict.lookupOrDefault<scalar>("x_plane", 0.1); 

// scalar Kd = flowControlDict.lookupOrDefault<scalar>("Kd", 0.0);

//*dimensionedScalar error_integral
//(
//    "error_integral",
//    dimensionSet(1, 0, -1, 0, 0, 0, 0), // 单位：m³
//    0.0  // 初始值
//);


// 定义 dimensionedVector
const dimensionedVector gravity("gravity", dimensionSet(0,1,-2,0,0,0,0), gravityValue);
const dimensionedVector drive("drive", dimensionSet(0,1,-2,0,0,0,0), driveValue);


// 输出调试信息
Info << "Gravity value from controlDict: " << gravityValue << endl;
Info << "Drive value from controlDict: " << driveValue << endl;

