#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mainheader.h"

double dt=0.05,time
      ,FinalTime=1000.0
      ,DensW=1.0//998.21//997.0
      ,DensA=1.205//1.185
      ,ViscW=1.0/5000.//0.001003//0.0008899
      ,ViscA=0.00001820755//0.00001831
      ,Ti=0.16*pow(1000.0,-1.0/8.0)   ///saman initial Kturb 0.0001 to 10
      ,Tm=10.0                  ///saman initial MUturb 0.1 to 100
      ,sigk=1.0///saman zaribe nuturb baraie K
      ,sigeps=1.3///saman zaribe nuturb baraie Eps
      ,ceps1=1.44///saman zaribe source term Eps
      ,ceps2=1.92;///saman zaribe source term Eps
double sig=0.5
      ,sigstar=0.5
      ,alf=5.0/9.0
      ,beta=3.0/40.0
      ,betastar=9.0/100.0;
struct MyVector Gravity = {0.0,0.0,0.0};
int NoOfTimeStep = 0
   ,SaveCounter = 100
   ,ResumeMode = 0;

char Dimension='2'; // 2D or 3D mesh :2D mesh should creat in Z-X Plan
double Mesh_Truncation=1.0e-10; // 2D or 3D mesh

char Epsequation='C';///saman switch hal moadele Eps
                     /// A -> Eps^2 AZ TIME GHABL
                     /// B -> Eps^2 YEKI AZ TIME GHABL
                     /// C -> Eps^2 YEKI AZ TIME GHABL VA ZARIBE E DAR SOURCE TERM
char Omequation='C';///saman switch hal moadele Om
                     /// A -> Om^2 AZ TIME GHABL
                     /// B -> Om^2 YEKI AZ TIME GHABL
                     /// C -> Om^2 YEKI AZ TIME GHABL VA ZARIBE E DAR SOURCE TERM
char WallFunction='B';///saman switch wall function
                     /// A -> wall function shomareie 1
                     /// B -> wall function shomareie 2


/// /////Laspack Library///////////
#include "./laspack/laspack.h"
/// /////Code Library//////////////
#include "general.h"
#include "mesh.h"
#include "nse.h"
#include "vof.h"
#include "kep.h"
#include "kom.h"
#include "rbe.h"
#include "exporter.h"
#include "cdfmanager.h"

void GlobalInitialization(struct Grid *Grids)
{
    int i;
    time=0.0;

    FillGrid(Grids,"cavity.neu");
//    FillGrid(Grids,"wedge.neu");
//    FillGrid(Grids,"naca0012_5degree.neu");

    InitialConditionForRigidBodyVariables(Grids);
    InitialConditionForFlowVariables(Grids);
    Write2VTK(Grids,"InitialCondition",0);
}

int main(void)
{
    int i;
    struct Grid Grids;
    double stime;
    char FileName[80];
    sprintf(FileName,"cavity");


    /// (START) export pressure vs time for following coordinates:
//    int j,k;
//    struct MyVector coordinate,*arrayVect;
//    k=20; ///no of coordinates
//    double a;
//    a=8./(double)(k+1);
//    arrayVect = (struct MyVector *)malloc(sizeof(struct MyVector)*k);
//    arrayVect[0].X=1.5;
//    arrayVect[0].Y=-0.5;
//    arrayVect[0].Z=-4+a;
//    for (j=1;j<k;j++)
//    {
//        CreateValue2CSV(FileName,j);
//        /// write the formula
//        arrayVect[j].X=1.5;
//        arrayVect[j].Y=-0.5;
//        arrayVect[j].Z=arrayVect[j-1].Z+a;
//    }
//    for (j=0;j<k;j++)
//        printf("%f, %f, %f\n",arrayVect[j].X,arrayVect[j].Y,arrayVect[j].Z);
//    printf("Press ENTER to continue");
//    getchar();
    /// (END)



    if (ResumeMode)
    {
        GlobalInitialization(&Grids);
        ReadFromCDF(&Grids,FileName);
    }
    else
    {
        GlobalInitialization(&Grids);
        WriteToCDF(&Grids,FileName);
        CreateTHF(FileName,Grids.NRBG);

    }
    stime=time;


    for(time=stime+dt;time<FinalTime;time+=dt)
    {
        NoOfTimeStep++;
        printf("\n_________________________________________________________________________\n");
        printf("Time = %f\n",time);
        double RBEError=0.0;

//        CellProp(&Grids,10);
        /// check the Boundary & inital Conditions first!
//        SolveKEEquation(&Grids,KEConvectionMethod);
//        SolveKOEquation(&Grids,KOConvectionMethod);

        SolveNSEquation(&Grids,ConvectionMethod
                              ,TemporalTermDiscretizationMethodForMomentum
                              ,TimeDiscretizationMethodForConvection
                              ,TimeDiscretizationMethodForDiffusion);

        ///saman ezafe shodane K - Epsilon
//        getchar();

//        VoFSolver(&Grids,VoFMethod
//                        ,TemporalTermDiscretizationMethodForVoF
//                        ,TimeDiscretizationMethodForVoF);

//        SolveRBEquation(&RBEError,&Grids);
//        printf("RBEError = %e\n",RBEError);

/// ///////////////////For LevelSet Method//////////////////////////////////////
//                LSSolver(&Grids
//                         ,TemporalTermDiscretizationMethodForMomentum
//                         ,TimeDiscretizationMethodForConvection);
//
//                int PCell;
//                for (PCell=0;PCell<Grids.NC;PCell++)
//                    Grids.Variables.NAlpha[PCell]=GetVolumeFraction(Grids.Variables.NPhi[PCell],Grids.Cells[PCell].Volume);
/// /////////////////////////////////////////////////////////////////////////////

//        AdvanceGrids(&Grids);
//        getchar();

        EqualSVect(Grids.Variables.PU,Grids.Variables.U ,Grids.NC);
        EqualSVect(Grids.Variables.PV,Grids.Variables.V ,Grids.NC);
        EqualSVect(Grids.Variables.PW,Grids.Variables.W ,Grids.NC);
        EqualSVect(Grids.Variables.U ,Grids.Variables.NU,Grids.NC);
        EqualSVect(Grids.Variables.V ,Grids.Variables.NV,Grids.NC);
        EqualSVect(Grids.Variables.W ,Grids.Variables.NW,Grids.NC);
        EqualSVect(Grids.Variables.P ,Grids.Variables.NP,Grids.NC);
        EqualSurfFlux(Grids.Variables.SF,Grids.Variables.NSF,&Grids);
        EqualSurfFlux(Grids.Variables.AF,Grids.Variables.NAF,&Grids);

        ///saman update maghadir K - Epsilon
        EqualSVect(Grids.Variables.PK,Grids.Variables.K ,Grids.NC);
        EqualSVect(Grids.Variables.K ,Grids.Variables.NK,Grids.NC);
        EqualSVect(Grids.Variables.PEps,Grids.Variables.Eps ,Grids.NC);
        EqualSVect(Grids.Variables.Eps ,Grids.Variables.NEps,Grids.NC);
        EqualSVect(Grids.Variables.POm,Grids.Variables.Om ,Grids.NC);
        EqualSVect(Grids.Variables.Om ,Grids.Variables.NOm,Grids.NC);

        EqualSVect(Grids.Variables.PAlpha,Grids.Variables.Alpha,Grids.NC);
        EqualSVect(Grids.Variables.Alpha,Grids.Variables.NAlpha,Grids.NC);
//        EqualSVect(Grids.Variables.PPhi,Grids.Variables.Phi,Grids.NC);
//        EqualSVect(Grids.Variables.Phi,Grids.Variables.NPhi,Grids.NC);

//        UpdateTHF(FileName,&Grids,Grids.NRBG);

        /// (START) update CSV files
//        for (j=0;j<k;j++)
//        {
//            coordinate.X=arrayVect[j].X;
//            coordinate.Y=arrayVect[j].Y;
//            coordinate.Z=arrayVect[j].Z;
//            UpdateCSV(FileName,&Grids,coordinate,j);
//        }
//        printf("%s.csv are updated Successfully.\n",FileName);
        /// (END)


        if (NoOfTimeStep%SaveCounter==0)
        {
            Write2VTK(&Grids,FileName,NoOfTimeStep/SaveCounter);
            WriteToCDF(&Grids,FileName);
        }
    }
    return(0);
}
