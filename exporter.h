double NodalValue(int NodeNumber,double *CellValue,struct Grid *Grids)
{
    int j,PCell;
    double InvVol,ValDivVol;

    InvVol=0.0;ValDivVol=0.0;
    for(j=0;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
    {
        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
        InvVol=InvVol+1.0/Grids->Cells[PCell].Volume;
        ValDivVol=ValDivVol+CellValue[PCell]/Grids->Cells[PCell].Volume;
    }
    return(ValDivVol/InvVol);
}

void FindBodyNodes(double *NodeO,struct Grid *Grids)
{
    int i,k,NumOfNodes_1,FaceNum,PCell;
    for(k=0;k<Grids->NN;k++)
        NodeO[k]=0.0;
    for(i=0;i<Grids->NRB;i++)
    {
        PCell=Grids->RigidBoundaryCells[i].CellNum;
        FaceNum=Grids->RigidBoundaryCells[i].FaceNum;
        switch (Grids->RCells[PCell].Type)
        {
            case 4:
                NumOfNodes_1=BREAK_NODES[FaceNum][0];
                for (k=0;k<NumOfNodes_1;k++)
                    NodeO[Grids->RCells[PCell].NodeList[BREAK_NODES[FaceNum][k+1]]]=1.0;
            break;
            case 6:
                NumOfNodes_1=TETRA_NODES[FaceNum][0];
                for (k=0;k<NumOfNodes_1;k++)
                    NodeO[Grids->RCells[PCell].NodeList[TETRA_NODES[FaceNum][k+1]]]=1.0;
            break;
            case 5:
                NumOfNodes_1=WEDGE_NODES[FaceNum][0];
                for (k=0;k<NumOfNodes_1;k++)
                    NodeO[Grids->RCells[PCell].NodeList[WEDGE_NODES[FaceNum][k+1]]]=1.0;
            break;
            case 7:
                NumOfNodes_1=PYRAMID_NODES[FaceNum][0];
                for (k=0;k<NumOfNodes_1;k++)
                    NodeO[Grids->RCells[PCell].NodeList[PYRAMID_NODES[FaceNum][k+1]]]=1.0;
            break;
        }
    }
}

double NodalValueAlpha(int NodeNumber,double *CellValue,struct Grid *Grids)
{
    int j,PCell;
    double AlphaPlus=0.0,AlphaMinus=0.0,Val,
           NumPlus=0.0,NumMinus=0.0;

    for(j=0;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
    {
        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
        Val=CellValue[PCell];
        if (Val>0.5)
        {
            AlphaPlus+=Val;
            NumPlus+=1.0;
        }
        else
        {
            AlphaMinus+=Val;
            NumMinus+=1.0;
        }
    }
    if (NumMinus==0)
    {
        AlphaMinus=AlphaPlus;
        NumMinus=NumPlus;
    }
    if (NumPlus==0)
    {
        AlphaPlus=AlphaMinus;
        NumPlus=NumMinus;
    }
    return (0.5*(AlphaMinus/NumMinus+AlphaPlus/NumPlus));

}

//double NodalValueAlpha(int NodeNumber,double *CellValue,struct Grid *Grids)
//{
//    int j,PCell;
//    double MaxVal,MinVal,Val;
//    double MaxVol,MinVol;
//
//    PCell=Grids->Nodes[NodeNumber].SharedCells[0];
//    MaxVal=CellValue[PCell];
//    MaxVol=Grids->Cells[PCell].Volume;
//    for(j=1;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
//    {
//        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
//        Val=CellValue[PCell];
//            if (Val>MaxVal)
//            {
//                MaxVal=Val;
//                MaxVol=Grids->Cells[PCell].Volume;
//            }
//    }
//
//    PCell=Grids->Nodes[NodeNumber].SharedCells[0];
//    MinVal=CellValue[PCell];
//    MinVol=Grids->Cells[PCell].Volume;
//    for(j=1;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
//    {
//        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
//        Val=CellValue[PCell];
//            if (Val<MinVal)
//            {
//                MinVal=Val;
//                MinVol=Grids->Cells[PCell].Volume;
//            }
//    }
//
//    return (MaxVal+MinVal)/2.0;
////    return (MaxVal/MaxVol+MinVal/MinVol)/(1./MaxVol+1./MinVol);
//}

double NodalValue_LSE(int NodeNumber,double *CellValue,struct MyVector *LSEGrad,struct Grid *Grids)
{
    int j,PCell;
    double InvDist = 0.0,
           TotalInvDist = 0.0,
           ValPerInvDist = 0.0;

    struct MyVector Cell2Node;
    for(j=0;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
    {
        PCell=Grids->Nodes[NodeNumber].SharedCells[j];
        Cell2Node=Minus(Grids->Nodes[NodeNumber].Pos,
                        Grids->Cells[PCell].CellCenter);
        InvDist = 1.0/sqrt(Dot(Cell2Node,Cell2Node));
        TotalInvDist += InvDist;
        ValPerInvDist += (CellValue[PCell]+Dot(LSEGrad[PCell],Cell2Node))*InvDist;
    }
    return(ValPerInvDist/TotalInvDist);

}

void  Write2VTK(struct Grid *Grids,char *Name,int Counter)
{
    int i;
    char path[80];
    FILE *VTK;

    sprintf(path,"./result/%s_%d.vtk",Name,Counter);
    if( (VTK=fopen(path,"w"))==NULL )
    {
        puts("vtk writing error.................Exiting ( Write2txt 000)");
        exit(-1);
    }

    int NN,NC;
    NN=Grids->NN;
    NC=Grids->NC;

//Creating Hedear
//////////////////////////////////////////
    fprintf(VTK,"# vtk DataFile Version 2.2\n");
    fprintf(VTK,"Unstructured Grid Example\n");
    fprintf(VTK,"ASCII\n");
    fprintf(VTK,"DATASET UNSTRUCTURED_GRID\n");
//////////////////////////////////////////
    fprintf(VTK,"POINTS %d double\n",NN);
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e %e %e\n",Grids->Nodes[i].Pos.X,Grids->Nodes[i].Pos.Y,Grids->Nodes[i].Pos.Z);
/////////////////////////////////////////
    int TotalCellData=0;
    for (i=0;i<NC;i++)
        TotalCellData+=Grids->RCells[i].NodePerCell;
    TotalCellData+=Grids->NC;

    fprintf(VTK,"CELLS %d %d\n",NC,TotalCellData);
    for (i=0;i<NC;i++)
    {

        fprintf(VTK,"%d ",Grids->RCells[i].NodePerCell);
        int *GAMBIT_2_VTK;
        switch (Grids->RCells[i].Type)
        {
            case 4: GAMBIT_2_VTK=GAMBIT_2_VTK_BREAK;   break;
            case 5: GAMBIT_2_VTK=GAMBIT_2_VTK_WEDGE;   break;
            case 6: GAMBIT_2_VTK=GAMBIT_2_VTK_TETRA;   break;
            case 7: GAMBIT_2_VTK=GAMBIT_2_VTK_PYRAMID; break;
        }
        int j;
        for(j=0;j<Grids->RCells[i].NodePerCell;j++)
            fprintf(VTK,"%d ",Grids->RCells[i].NodeList[GAMBIT_2_VTK[j]]);
        fprintf(VTK,"\n");

    }

/////////////////////////////////////////
    fprintf(VTK,"CELL_TYPES %d\n",NC);
    for (i=0;i<NC;i++)
    {
        switch (Grids->RCells[i].Type)
        {
            case 4: fprintf(VTK,"%d\n",12); break;
            case 5: fprintf(VTK,"%d\n",13); break;
            case 6: fprintf(VTK,"%d\n",10); break;
            case 7: fprintf(VTK,"%d\n",14); break;
        }
    }

/////////////////////////////////////////
    fprintf(VTK,"CELL_DATA %d\n",NC);

    fprintf(VTK,"SCALARS U float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NU[i]);

    fprintf(VTK,"SCALARS V float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NV[i]);

    fprintf(VTK,"SCALARS W float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NW[i]);

    fprintf(VTK,"SCALARS P float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NP[i]);

    fprintf(VTK,"SCALARS K float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NK[i]);

    fprintf(VTK,"SCALARS E float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NEps[i]);

    fprintf(VTK,"SCALARS Omega float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NOm[i]);

    fprintf(VTK,"SCALARS A float 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f\n",Grids->Variables.NAlpha[i]);

    fprintf(VTK,"VECTORS Vel_Vec float\n");
      for (i=0;i<NC;i++)
        fprintf(VTK,"%.15f %.15f %.15f\n",Grids->Variables.NU[i],Grids->Variables.NV[i],Grids->Variables.NW[i]);
/////////////////////////////////////////


/////////////////////////////////////////
    fprintf(VTK,"POINT_DATA %d\n",NN);
///////////////////////////////////////
    fprintf(VTK,"SCALARS U double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids));
///////////////////////////////////////
    fprintf(VTK,"SCALARS V double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids));
///////////////////////////////////////
    fprintf(VTK,"SCALARS W double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids));
///////////////////////////////////////
    fprintf(VTK,"SCALARS P double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NP,Grids->Variables.PressLSEGrad,Grids));
///////////////////////////////////////
    fprintf(VTK,"SCALARS K double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NK,Grids->Variables.KLSEGrad,Grids));
//////////////////////////////////////
    fprintf(VTK,"SCALARS E double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NEps,Grids->Variables.EpsLSEGrad,Grids));
//////////////////////////////////////
    fprintf(VTK,"SCALARS Omega double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValue_LSE(i,Grids->Variables.NOm,Grids->Variables.OmLSEGrad,Grids));
//////////////////////////////////////
    fprintf(VTK,"SCALARS A double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodalValueAlpha(i,Grids->Variables.NAlpha,Grids));
///////////////////////////////////////
    double *NodeO;
    NodeO=(double *)malloc(sizeof(double)*NN);
    FindBodyNodes(NodeO,Grids);
    fprintf(VTK,"SCALARS O double 1\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e\n",NodeO[i]);
    free(NodeO);
///////////////////////////////////////
    fprintf(VTK,"VECTORS Vel double\n");
    for (i=0;i<NN;i++)
        fprintf(VTK,"%e %e %e\n",NodalValue_LSE(i,Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids),
                                 NodalValue_LSE(i,Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids),
                                 NodalValue_LSE(i,Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids));
///////////////////////////////////////
    fprintf(VTK,"CELL_DATA %d\n",NC);
///////////////////////////////////////
    fprintf(VTK,"SCALARS AC double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        fprintf(VTK,"%e\n",Grids->Variables.Alpha[i]);
///////////////////////////////////////
    fprintf(VTK,"SCALARS P_C double\n");
    fprintf(VTK,"LOOKUP_TABLE default\n");
    for (i=0;i<NC;i++)
        fprintf(VTK,"%e\n",Grids->Variables.NP[i]);
///////////////////////////////////////

    fclose(VTK);
    printf("Write to %s_%d.vtk Completed!\n",Name,Counter);
}


void CreateTHF(char *Name, int Counter)
{
    char path[80],dummyStr[80],dummyStr2[80];
    int i;
    FILE *dummyTHF;
    for(i=0;i<Counter;i++)
    {
        sprintf(dummyStr,"");
        sprintf(dummyStr2,"");
        sprintf(path,"./result/%s_%d.thf",Name,i);
        if( (dummyTHF=fopen(path,"w"))==NULL )
        {
            puts("THF Creating error.................Exiting ( CreateTHF() 000)");
            exit(-1);
        }
        fprintf(dummyTHF,"Time,PFX,PFY,PFZ,VFX,VFY,VFZ,MX,MY,MZ,AX,AY,AZ,AlX,AlY,AlZ,UX,UY,UZ,OmX,OmY,OmZ,XX,XY,XZ,OX,OY,OZ\n");
        fclose(dummyTHF);
        printf("%s_%d.thf Created Succesfully.\n",Name,i);
    }
}
void UpdateTHF(char *Name,struct Grid *Grids, int GroupNo)
{
    int i;
    char path[80];
    FILE *dummyTHF;
    struct MyVector PForce,VForce,Moment,LinAcc,AngAcc,LinVel,AngVel,MassCenter,BodyOr;

    for(i=0;i<GroupNo;i++)
    {
        sprintf(path,"./result/%s_%d.thf",Name,i);
        if( (dummyTHF=fopen(path,"a"))==NULL )
        {
            puts("THF Updating error.................Exiting ( UpdateTHF() 000)");
            exit(-1);
        }
        VForce=CalcViscousForces(Grids,i);
        PForce=CalcPressureForces(Grids,i);
        Moment=Sum(CalcPressureMoments(Grids,i),CalcViscousMoments(Grids,i));
        LinAcc=Grids->DynaVariables.NewLinAccel;
        AngAcc=Grids->DynaVariables.NewAngAccel;
        LinVel=Grids->DynaVariables.NewLinVel;
        AngVel=Grids->DynaVariables.NewAngVel;
        MassCenter=Grids->DynaVariables.NewMassCenter;
        BodyOr=Grids->DynaVariables.NewBodyOr;

        printf("\nPressure Force : ( %e , %e , %e )\n",PForce.X,PForce.Y,PForce.Z);
        printf("Viscous  Force : ( %e , %e , %e )\n",VForce.X,VForce.Y,VForce.Z);
    //    printf("Induced Moment : ( %e , %e , %e )\n",Moment.X,Moment.Y,Moment.Z);
        fprintf(dummyTHF,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n"
                        ,time
                        ,PForce.X,PForce.Y,PForce.Z
                        ,VForce.X,VForce.Y,VForce.Z
                        ,Moment.X,Moment.Y,Moment.Z
                        ,LinAcc.X,LinAcc.Y,LinAcc.Z
                        ,AngAcc.X,AngAcc.Y,AngAcc.Z
                        ,LinVel.X,LinVel.Y,LinVel.Z
                        ,AngVel.X,AngVel.Y,AngVel.Z
                        ,MassCenter.X,MassCenter.Y,MassCenter.Z
                        ,BodyOr.X,BodyOr.Y,BodyOr.Z);
        fclose(dummyTHF);
        printf("%s_%d.thf Updated Succesfully.\n",Name,i);
    }
}

//void CreateTHF(char *Name)
//{
//    char path[80],dummyStr[80],dummyStr2[80];
//    int i;
//    FILE *dummyTHF;
//
//    sprintf(dummyStr,"");
//    sprintf(dummyStr2,"");
//    sprintf(path,"./result/%s.thf",Name);
//    if( (dummyTHF=fopen(path,"w"))==NULL )
//    {
//        puts("THF Creating error.................Exiting ( CreateTHF() 000)");
//        exit(-1);
//    }
//    fprintf(dummyTHF,"Time,PFX,PFY,PFZ,VFX,VFY,VFZ,MX,MY,MZ,AX,AY,AZ,AlX,AlY,AlZ,UX,UY,UZ,OmX,OmY,OmZ,XX,XY,XZ,OX,OY,OZ\n");
//    fclose(dummyTHF);
//    printf("%s.thf Created Succesfully.\n",Name);
//}

//void UpdateTHF(char *Name,struct Grid *Grids)
//{
//    int i,j;
//    char path[80];
//    FILE *dummyTHF;
//
//    sprintf(path,"./result/%s.thf",Name);
//    if( (dummyTHF=fopen(path,"a"))==NULL )
//    {
//        puts("THF Updating error.................Exiting ( UpdateTHF() 000)");
//        exit(-1);
//    }
//
//    struct MyVector PForce,VForce,Moment,LinAcc,AngAcc,LinVel,AngVel,MassCenter,BodyOr;
//    VForce=CalcViscousForces(Grids);
//    PForce=CalcPressureForces(Grids);
//    Moment=Sum(CalcPressureMoments(Grids),CalcViscousMoments(Grids));
//    LinAcc=Grids->DynaVariables.NewLinAccel;
//    AngAcc=Grids->DynaVariables.NewAngAccel;
//    LinVel=Grids->DynaVariables.NewLinVel;
//    AngVel=Grids->DynaVariables.NewAngVel;
//    MassCenter=Grids->DynaVariables.NewMassCenter;
//    BodyOr=Grids->DynaVariables.NewBodyOr;
//
//    printf("Pressure Force : ( %e , %e , %e )\n",PForce.X,PForce.Y,PForce.Z);
//    printf("Viscous  Force : ( %e , %e , %e )\n",VForce.X,VForce.Y,VForce.Z);
////    printf("Induced Moment : ( %e , %e , %e )\n",Moment.X,Moment.Y,Moment.Z);
//    fprintf(dummyTHF,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n"
//                    ,time
//                    ,PForce.X,PForce.Y,PForce.Z
//                    ,VForce.X,VForce.Y,VForce.Z
//                    ,Moment.X,Moment.Y,Moment.Z
//                    ,LinAcc.X,LinAcc.Y,LinAcc.Z
//                    ,AngAcc.X,AngAcc.Y,AngAcc.Z
//                    ,LinVel.X,LinVel.Y,LinVel.Z
//                    ,AngVel.X,AngVel.Y,AngVel.Z
//                    ,MassCenter.X,MassCenter.Y,MassCenter.Z
//                    ,BodyOr.X,BodyOr.Y,BodyOr.Z);
//
//    fclose(dummyTHF);
//    printf("%s.thf Updated Succesfully.\n",Name);
//}
void CreateValue2CSV(char *Name,int i)
{
    char path[80],dummyStr[80],dummyStr2[80];
    FILE *dummyTHF;

        sprintf(dummyStr,"");
        sprintf(dummyStr2,"");
        sprintf(path,"./result/%s_%d.csv",Name,i);
        if( (dummyTHF=fopen(path,"w"))==NULL )
        {
            puts("THF Creating error.................Exiting ( CreateCSV() 000)");
            exit(-1);
        }
        fprintf(dummyTHF,"Time,Pressure\n");
        fclose(dummyTHF);
//        printf("%s_%d.csv Created Succesfully.\n",Name,i);

}
void UpdateCSV(char *Name,struct Grid *Grids,struct MyVector coordinate,int i)
{
    int cellNo;
    char path[80];
    double pressure;
    FILE *dummyCSV;

    cellNo=FindNearestCellTo(coordinate,Grids);
    pressure=Grids->Variables.NP[cellNo];
            +Dot(Minus(coordinate,Grids->Cells[cellNo].CellCenter),Grids->Variables.PressLSEGrad[cellNo]);

    sprintf(path,"./result/%s_%d.csv",Name,i);

    if( (dummyCSV=fopen(path,"a"))==NULL )
    {
        puts("CSV Updating error.................Exiting ( UpdateCSV() 000)");
        exit(-1);
    }

    fprintf(dummyCSV,"%e,%e\n",time,pressure,Grids->Variables.NP[cellNo]);
    fclose(dummyCSV);
//    printf("%s_%d.csv Updated Succesfully.\n",Name,i);
}
