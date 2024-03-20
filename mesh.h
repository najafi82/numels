void MallocFlowVariables(struct Grid *Grids)
{
    Grids->Variables.U  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.V  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.W  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NU = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NV = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NW = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PU = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PV = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PW = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.Alpha   = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NAlpha  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PAlpha  = (double *)malloc(sizeof(double)*Grids->NC);

    ///saman
    Grids->Variables.NEps  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.Eps   = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PEps  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NOm  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.Om   = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.POm  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NK    = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.K     = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PK    = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NuTurb    = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.KLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.EpsLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.OmLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);



//    Grids->Variables.Phi  = (double *)malloc(sizeof(double)*Grids->NC);
//    Grids->Variables.NPhi = (double *)malloc(sizeof(double)*Grids->NC);
//    Grids->Variables.PPhi = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.P  = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.NP = (double *)malloc(sizeof(double)*Grids->NC);
//    Grids->Variables.Dummy = (double *)malloc(sizeof(double)*Grids->NC);
    Grids->Variables.PressLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.ULSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.VLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.WLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
//    Grids->Variables.PhiLSEGrad = (struct MyVector *)malloc(sizeof(struct MyVector)*Grids->NC);
    Grids->Variables.SF = (struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
    Grids->Variables.NSF = (struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
    Grids->Variables.AF = (struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
    Grids->Variables.NAF = (struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
    int PCell,NoOfFace;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        Grids->Variables.SF [PCell].Face=(double*)malloc(sizeof(double)*NoOfFace);
        Grids->Variables.NSF[PCell].Face=(double*)malloc(sizeof(double)*NoOfFace);
        Grids->Variables.AF [PCell].Face=(double*)malloc(sizeof(double)*NoOfFace);
        Grids->Variables.NAF[PCell].Face=(double*)malloc(sizeof(double)*NoOfFace);
    }
}

void InitFlowVariables(struct Grid *Grids)
{
    int PCell,faceNo;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.U[PCell]=0.0;
        Grids->Variables.V[PCell]=0.0;
        Grids->Variables.W[PCell]=0.0;
        Grids->Variables.NU[PCell]=0.0;
        Grids->Variables.NV[PCell]=0.0;
        Grids->Variables.NW[PCell]=0.0;
        Grids->Variables.PU[PCell]=0.0;
        Grids->Variables.PV[PCell]=0.0;
        Grids->Variables.PW[PCell]=0.0;
        Grids->Variables.Alpha[PCell]=0.0;
        Grids->Variables.NAlpha[PCell]=0.0;
        Grids->Variables.PAlpha[PCell]=0.0;
//        Grids->Variables.Phi[PCell]=0.0;
//        Grids->Variables.NPhi[PCell]=0.0;
//        Grids->Variables.PPhi[PCell]=0.0;
        Grids->Variables.P[PCell]=0.0;
        Grids->Variables.NP[PCell]=0.0;
//        Grids->Variables.Dummy[PCell]=0.0;
        Grids->Variables.PressLSEGrad[PCell].X=
        Grids->Variables.PressLSEGrad[PCell].Y=
        Grids->Variables.PressLSEGrad[PCell].Z=0.0;
        Grids->Variables.ULSEGrad[PCell].X=
        Grids->Variables.ULSEGrad[PCell].Y=
        Grids->Variables.ULSEGrad[PCell].Z=0.0;
        Grids->Variables.VLSEGrad[PCell].X=
        Grids->Variables.VLSEGrad[PCell].Y=
        Grids->Variables.VLSEGrad[PCell].Z=0.0;
        Grids->Variables.WLSEGrad[PCell].X=
        Grids->Variables.WLSEGrad[PCell].Y=
        Grids->Variables.WLSEGrad[PCell].Z=0.0;
//        Grids->Variables.PhiLSEGrad[PCell].X=
//        Grids->Variables.PhiLSEGrad[PCell].Y=
//        Grids->Variables.PhiLSEGrad[PCell].Z=0.0;

        for (faceNo=0;faceNo<Grids->Cells[PCell].NoOfFace;faceNo++)
        {
            Grids->Variables.SF [PCell].Face[faceNo]=0.0;
            Grids->Variables.NSF[PCell].Face[faceNo]=0.0;
            Grids->Variables.AF [PCell].Face[faceNo]=0.0;
            Grids->Variables.NAF[PCell].Face[faceNo]=0.0;
        }
    }
}

void InitDynaVariables(struct Grid *Grids)
{
    Grids->DynaVariables.Mass = 0.0;
    Grids->DynaVariables.Inertia[0][0] = 0.0;
    Grids->DynaVariables.Inertia[0][1] = 0.0;
    Grids->DynaVariables.Inertia[0][2] = 0.0;
    Grids->DynaVariables.Inertia[1][0] = 0.0;
    Grids->DynaVariables.Inertia[1][1] = 0.0;
    Grids->DynaVariables.Inertia[1][2] = 0.0;
    Grids->DynaVariables.Inertia[2][0] = 0.0;
    Grids->DynaVariables.Inertia[2][1] = 0.0;
    Grids->DynaVariables.Inertia[2][2] = 0.0;

    Grids->DynaVariables.NewLinAccel.X = 0.0;
    Grids->DynaVariables.NewLinAccel.Y = 0.0;
    Grids->DynaVariables.NewLinAccel.Z = 0.0;

    Grids->DynaVariables.NewAngAccel.X = 0.0;
    Grids->DynaVariables.NewAngAccel.Y = 0.0;
    Grids->DynaVariables.NewAngAccel.Z = 0.0;

    Grids->DynaVariables.NewLinVel.X = 0.0;
    Grids->DynaVariables.NewLinVel.Y = 0.0;
    Grids->DynaVariables.NewLinVel.Z = 0.0;

    Grids->DynaVariables.NewAngVel.X = 0.0;
    Grids->DynaVariables.NewAngVel.Y = 0.0;
    Grids->DynaVariables.NewAngVel.Z = 0.0;

    Grids->DynaVariables.NewMassCenter.X = 0.0;
    Grids->DynaVariables.NewMassCenter.Y = 0.0;
    Grids->DynaVariables.NewMassCenter.Z = 0.0;

    Grids->DynaVariables.NewBodyOr.X = 0.0;
    Grids->DynaVariables.NewBodyOr.Y = 0.0;
    Grids->DynaVariables.NewBodyOr.Z = 0.0;

    Grids->DynaVariables.LinAccel.X = 0.0;
    Grids->DynaVariables.LinAccel.Y = 0.0;
    Grids->DynaVariables.LinAccel.Z = 0.0;

    Grids->DynaVariables.AngAccel.X = 0.0;
    Grids->DynaVariables.AngAccel.Y = 0.0;
    Grids->DynaVariables.AngAccel.Z = 0.0;

    Grids->DynaVariables.LinVel.X = 0.0;
    Grids->DynaVariables.LinVel.Y = 0.0;
    Grids->DynaVariables.LinVel.Z = 0.0;

    Grids->DynaVariables.AngVel.X = 0.0;
    Grids->DynaVariables.AngVel.Y = 0.0;
    Grids->DynaVariables.AngVel.Z = 0.0;

    Grids->DynaVariables.MassCenter.X = 0.0;
    Grids->DynaVariables.MassCenter.Y = 0.0;
    Grids->DynaVariables.MassCenter.Z = 0.0;

    Grids->DynaVariables.BodyOr.X = 0.0;
    Grids->DynaVariables.BodyOr.Y = 0.0;
    Grids->DynaVariables.BodyOr.Z = 0.0;
}

void FillNodes(struct Grid *Grids,char *FileName)//coordinates of nodes has been read
{
    FILE *fpMesh;
    char dummystr[80];
    int i=0,NN;
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("NodeCoord.ebr reading error.................Exiting( ReadNodeCoord 000)");
        exit(-1);
    }
    GetToken(dummystr,fpMesh);
    if ((*Grids).NN!=atoi(dummystr))
    {
        puts("Number of nodes don't match");
        exit(-1);
    }
    for(i=0;i<(*Grids).NN;i++)
    {
        GetToken(dummystr,fpMesh);
        NN=atoi(dummystr)-1;
        GetToken(dummystr,fpMesh);
        (*Grids).Nodes[NN].Pos.X=atof(dummystr);
        GetToken(dummystr,fpMesh);
        (*Grids).Nodes[NN].Pos.Y=atof(dummystr);
        GetToken(dummystr,fpMesh);
        (*Grids).Nodes[NN].Pos.Z=atof(dummystr);
    }
    fclose(fpMesh);
}

void FillConnectivity(struct Grid *aGrid,char *FileName)
{
    FILE *fpMesh;
    char dummystr[80];
    int i=0,j=0,N=0;
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("ElemConnect.ebr reading error.................Exiting (ReadElementConnect 000)");
        exit(-1);
    }

    GetToken(dummystr,fpMesh);
    if (aGrid->NC!=atoi(dummystr))
    {
        puts("Number of elements don't match\n");
        exit(-1);
    }

    for (i=0;i<aGrid->NC;i++)
    {
        GetToken(dummystr,fpMesh);
        aGrid->RCells[i].Type=atoi(dummystr);
        GetToken(dummystr,fpMesh);
        aGrid->RCells[i].NodePerCell=atoi(dummystr);
        aGrid->RCells[i].NodeList=(int *)malloc(sizeof(int)*aGrid->RCells[i].NodePerCell);
        for (j=0;j<aGrid->RCells[i].NodePerCell;j++)
        {
            GetToken(dummystr,fpMesh);
            aGrid->RCells[i].NodeList[j]=atoi(dummystr)-1;
        }
    }
    fclose(fpMesh);
}

void FillRigidBoundary(struct Grid *Grids)
{
    int i,k=0;
    Grids->RigidBoundaryCells=(struct RigidBoundaryCell *)malloc(sizeof(struct RigidBoundaryCell));
    for (i=0;i<Grids->NB;i++)
    {
        if (Grids->BoundaryCells[i].Type==1)
        {
            Grids->RigidBoundaryCells=(struct RigidBoundaryCell*)realloc(Grids->RigidBoundaryCells,sizeof(struct RigidBoundaryCell)*(k+1));
            Grids->RigidBoundaryCells[k].CellNum=Grids->BoundaryCells[i].CellNum;
            Grids->RigidBoundaryCells[k].FaceNum=Grids->BoundaryCells[i].FaceNum;
            k++;
        }
    }
    Grids->NRB=k;

    k=0;
    for (i=0;i<Grids->NRBG;i++)
        k+=Grids->RBNo[i];
    if (Grids->NRB!=k)
    {
        printf("error in Rigid Body Number\n");
        getchar();
    }
}

void FillBoundary(struct Grid *Grids,char *FileName)
{
    FILE *fpMesh;
    char dummystr[80];
    int i=0;
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("Boundary.ebr reading error.................Exiting( ReadBoundaryCondition 000)");
        exit(-1);
    }

    for (i=0;i<Grids->NB;i++)
    {
        GetToken(dummystr,fpMesh);
        Grids->BoundaryCells[i].CellNum=atoi(dummystr)-1;
        GetToken(dummystr,fpMesh);
        Grids->BoundaryCells[i].Type=atoi(dummystr);
        GetToken(dummystr,fpMesh);
        Grids->BoundaryCells[i].FaceNum=atoi(dummystr)-1;
        GetToken(dummystr,fpMesh);
        Grids->BoundaryCells[i].Value=atof(dummystr);
    }
    fclose(fpMesh);
}

void FillNodeSharing(struct Grid *Grids)
{
    int i,j,k,N;
    for (i=0;i<Grids->NN;i++)
        Grids->Nodes[i].NumOfSharedCells=0;

    for (i=0;i<Grids->NC;i++)
        for (j=0;j<Grids->RCells[i].NodePerCell;j++)
            Grids->Nodes[Grids->RCells[i].NodeList[j]].NumOfSharedCells+=1;

    for (i=0;i<Grids->NN;i++)
    {
        Grids->Nodes[i].SharedCells=(int *)malloc(Grids->Nodes[i].NumOfSharedCells*sizeof(int));
        SetIntSVect(Grids->Nodes[i].SharedCells,Grids->Nodes[i].NumOfSharedCells,-1);
    }

    for (i=0;i<Grids->NC;i++)
        for (j=0;j<Grids->RCells[i].NodePerCell;j++)
        {
            N=Grids->RCells[i].NodeList[j];
            for (k=0;k<Grids->Nodes[N].NumOfSharedCells;k++)
                if (Grids->Nodes[N].SharedCells[k]==-1)
                    break;
            Grids->Nodes[N].SharedCells[k]=i;
        }
}

struct MyVector ComputeCellCenter(struct Cell *Cells) //compute based on BLAZIK
{
    struct MyVector CellCenter,up,dummyA,dummyF,dummySVP;
    int i;
    double down,dummyMag;
    down=0.0;
    up.X=up.Y=up.Z=0.0;

    for (i=0;i<Cells->NoOfFace;i++)
    {
        up=Sum(up,SVP((Dot(Cells->FaceCenter[i]
                     ,SVP(1./Mag(Cells->Area[i]),Cells->Area[i])))*Mag(Cells->Area[i])
                ,Cells->FaceCenter[i]));

        down+=Dot(Cells->FaceCenter[i]
                  ,SVP(1./Mag(Cells->Area[i]),Cells->Area[i]))*4.0*Mag(Cells->Area[i]);
    }
    CellCenter=SVP(3./down,up);
    if (Dimension=='2')
        CellCenter.Y=0.5;
    return(CellCenter);
}

void BrickShape(struct Cell *Cells,struct RawCell RCell,struct Node *Nodes)
{
    int faceNo,NodeNum,i;
    int NoFaces=6;
    struct MyVector dummyVect;
    Cells->NoOfFace=NoFaces;
    dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;

    for (i=0;i<RCell.NodePerCell;i++)
        dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[i]].Pos);
    Cells->CellCenter=SVP(1./(double)(RCell.NodePerCell),dummyVect);

    Cells->Neighbor=(int*)malloc(sizeof(int)*NoFaces);
    Cells->Area=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);

    Cells->FaceCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCoincide=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->ProjectedCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->Volume=0.0;
    for (faceNo=0;faceNo<NoFaces;faceNo++)
    {
        NodeNum=BREAK_NODES[faceNo][0];

        struct MyVector FaceCenter1,FaceCenter2;
        double Area1,Area2;
        FaceCenter1=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos,
                                 Sum(Nodes[RCell.NodeList[BREAK_NODES[faceNo][2]]].Pos,
                                     Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos)));
        FaceCenter2=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos,
                                 Sum(Nodes[RCell.NodeList[BREAK_NODES[faceNo][4]]].Pos,
                                     Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos)));
        Area1=Mag(Cross(Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos
                             ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][2]]].Pos)
                       ,Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos
                             ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos)));
        Area2=Mag(Cross(Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos
                             ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][4]]].Pos)
                       ,Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos
                             ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos)));
        Cells->FaceCenter[faceNo]=SVP(1.0/(Area1+Area2)
                                     ,Sum(SVP(Area1,FaceCenter1)
                                         ,SVP(Area2,FaceCenter2)));

        Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][3]]].Pos
                                               ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][1]]].Pos)
                                         ,Minus(Nodes[RCell.NodeList[BREAK_NODES[faceNo][4]]].Pos
                                               ,Nodes[RCell.NodeList[BREAK_NODES[faceNo][2]]].Pos)));
        Cells->Volume+=(1./3.)*Dot(Cells->FaceCenter[faceNo]
                                  ,Cells->Area[faceNo]);
    }
//    Cells->CellCenter=ComputeCellCenter(Cells);
}

void WedgeShape(struct Cell *Cells,struct RawCell RCell,struct Node *Nodes)
{
    int faceNo,NodeNum,i;
    int NoFaces=5;
    struct MyVector dummyVect;
    Cells->NoOfFace=NoFaces;
    dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;

//    for (i=0;i<RCell.NodePerCell;i++)
//        dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[i]].Pos);
//    Cells->CellCenter=SVP(1./(double)(RCell.NodePerCell),dummyVect);

    Cells->Neighbor=(int*)malloc(sizeof(int)*NoFaces);
    Cells->Area=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCoincide=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->ProjectedCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->Volume=0.0;
    for (faceNo=0;faceNo<NoFaces;faceNo++)
    {
        NodeNum=WEDGE_NODES[faceNo][0];
        if ((faceNo!=3)&&(faceNo!=4))
        {
            struct MyVector FaceCenter1,FaceCenter2;
            double Area1,Area2;
            FaceCenter1=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos,
                                     Sum(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][2]]].Pos,
                                         Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos)));
            FaceCenter2=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos,
                                     Sum(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][4]]].Pos,
                                         Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos)));
            Area1=Mag(Cross(Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos
                                 ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][2]]].Pos)
                           ,Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos
                                 ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos)));
            Area2=Mag(Cross(Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos
                                 ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][4]]].Pos)
                           ,Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos
                                 ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos)));
            Cells->FaceCenter[faceNo]=SVP(1.0/(Area1+Area2)
                                         ,Sum(SVP(Area1,FaceCenter1)
                                             ,SVP(Area2,FaceCenter2)));
        }
        else
        {
            dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;
            for (i=0;i<NodeNum;i++)
                dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][i+1]]].Pos);
            Cells->FaceCenter[faceNo]=SVP(1./(double)(NodeNum),dummyVect);
        }

        if ((faceNo!=3)&&(faceNo!=4))
            Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos
                                                   ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos)
                                             ,Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][4]]].Pos
                                                   ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][2]]].Pos)));
        else
            Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][2]]].Pos
                                                   ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos)
                                             ,Minus(Nodes[RCell.NodeList[WEDGE_NODES[faceNo][3]]].Pos
                                                   ,Nodes[RCell.NodeList[WEDGE_NODES[faceNo][1]]].Pos)));
        Cells->Volume+=(1./3.)*Dot(Cells->FaceCenter[faceNo]
                                  ,Cells->Area[faceNo]);
    }

//    double dif;
//    dif=Mag(Cells->Area[3])-Mag(Cells->Area[4]);
//    if(dif>1.0e-20)
//    {
//        printf("area dif is : %e",dif);
//        getchar();
//    }

//    Cells->CellCenter.X=ComputeCellCenter(Cells).X;
//    Cells->CellCenter.Y=ComputeCellCenter(Cells).Y;
//    Cells->CellCenter.Z=ComputeCellCenter(Cells).Z;
    Cells->CellCenter=ComputeCellCenter(Cells);

}

void TetraShape(struct Cell *Cells,struct RawCell RCell,struct Node *Nodes)
{
    int faceNo,NodeNum,i;
    int NoFaces=4;
    struct MyVector dummyVect;
    Cells->NoOfFace=NoFaces;
    dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;

//    for (i=0;i<RCell.NodePerCell;i++)
//        dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[i]].Pos);
//    Cells->CellCenter=SVP(1./(double)(RCell.NodePerCell),dummyVect);

    Cells->Neighbor=(int*)malloc(sizeof(int)*NoFaces);
    Cells->Area=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCoincide=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->ProjectedCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->Volume=0.0;
    for (faceNo=0;faceNo<NoFaces;faceNo++)
    {
        NodeNum=TETRA_NODES[faceNo][0];
        dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;
        for (i=0;i<NodeNum;i++)
            dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[TETRA_NODES[faceNo][i+1]]].Pos);
        Cells->FaceCenter[faceNo]=SVP(1./(double)(NodeNum),dummyVect);

        Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[TETRA_NODES[faceNo][2]]].Pos
                                               ,Nodes[RCell.NodeList[TETRA_NODES[faceNo][1]]].Pos)
                                         ,Minus(Nodes[RCell.NodeList[TETRA_NODES[faceNo][3]]].Pos
                                               ,Nodes[RCell.NodeList[TETRA_NODES[faceNo][1]]].Pos)));
        Cells->Volume+=(1./3.)*Dot(Cells->FaceCenter[faceNo]
                                  ,Cells->Area[faceNo]);
    }
    Cells->CellCenter=ComputeCellCenter(Cells);
}

void PyramidShape(struct Cell *Cells,struct RawCell RCell,struct Node *Nodes)
{
    int faceNo,NodeNum,i;
    int NoFaces=5;
    struct MyVector dummyVect;
    Cells->NoOfFace=NoFaces;
    dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;
//    for (i=0;i<RCell.NodePerCell;i++)
//        dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[i]].Pos);
//    Cells->CellCenter=SVP(1./(double)(RCell.NodePerCell),dummyVect);
    Cells->Neighbor=(int*)malloc(sizeof(int)*NoFaces);
    Cells->Area=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->FaceCoincide=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->ProjectedCenter=(struct MyVector*)malloc(sizeof(struct MyVector)*NoFaces);
    Cells->Volume=0.0;
    for (faceNo=0;faceNo<NoFaces;faceNo++)
    {
        NodeNum=PYRAMID_NODES[faceNo][0];
        dummyVect.X=dummyVect.Y=dummyVect.Z=0.0;
        if (faceNo==0)
        {
            struct MyVector FaceCenter1,FaceCenter2;
            double Area1,Area2;
            FaceCenter1=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos,
                                     Sum(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][2]]].Pos,
                                         Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos)));
            FaceCenter2=SVP(1.0/3.0, Sum(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos,
                                     Sum(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][4]]].Pos,
                                         Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos)));
            Area1=Mag(Cross(Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos
                                 ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][2]]].Pos)
                           ,Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos
                                 ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos)));
            Area2=Mag(Cross(Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos
                                 ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][4]]].Pos)
                           ,Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos
                                 ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos)));
            Cells->FaceCenter[faceNo]=SVP(1.0/(Area1+Area2)
                                         ,Sum(SVP(Area1,FaceCenter1)
                                             ,SVP(Area2,FaceCenter2)));
        }
        else
        {
            for (i=0;i<NodeNum;i++)
                dummyVect=Sum(dummyVect,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][i+1]]].Pos);
            Cells->FaceCenter[faceNo]=SVP(1./(double)(NodeNum),dummyVect);
        }

        if (faceNo==0)
            Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos
                                                   ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos)
                                             ,Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][4]]].Pos
                                                   ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][2]]].Pos)));
        else
            Cells->Area[faceNo]=SVP(0.5,Cross(Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][2]]].Pos
                                                   ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos)
                                             ,Minus(Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][3]]].Pos
                                                   ,Nodes[RCell.NodeList[PYRAMID_NODES[faceNo][1]]].Pos)));
        Cells->Volume+=(1./3.)*Dot(Cells->FaceCenter[faceNo]
                                  ,Cells->Area[faceNo]);
    }
    Cells->CellCenter=ComputeCellCenter(Cells);
}

void  CompCellGeoPara(struct Grid *aGrid)
{
    int PCell;
    for(PCell=0;PCell<aGrid->NC;PCell++)
    {
        switch (aGrid->RCells[PCell].Type)
        {
            case 4:
                BrickShape(&(aGrid->Cells[PCell]),aGrid->RCells[PCell],aGrid->Nodes);
                break;
            case 5:
                WedgeShape(&(aGrid->Cells[PCell]),aGrid->RCells[PCell],aGrid->Nodes);
                break;
            case 6:
                TetraShape(&(aGrid->Cells[PCell]),aGrid->RCells[PCell],aGrid->Nodes);
                break;
            case 7:
                PyramidShape(&(aGrid->Cells[PCell]),aGrid->RCells[PCell],aGrid->Nodes);
                break;
        }
    }

}

void  CompleteFacePara(struct Grid *aGrid)
{
    int PCell,FaceNum,NgbCell;
    for(PCell=0;PCell<aGrid->NC;PCell++)
        for(FaceNum=0;FaceNum<aGrid->Cells[PCell].NoOfFace;FaceNum++)
        {
            NgbCell=aGrid->Cells[PCell].Neighbor[FaceNum];
            if (NgbCell>=0)
                aGrid->Cells[PCell].FaceCoincide[FaceNum] =
                GetCoincidance(aGrid->Cells[PCell].CellCenter,
                               aGrid->Cells[NgbCell].CellCenter,
                               aGrid->Cells[PCell].FaceCenter[FaceNum],
                               Normalized(aGrid->Cells[PCell].Area[FaceNum]));
           else
                aGrid->Cells[PCell].FaceCoincide[FaceNum] =
                aGrid->Cells[PCell].FaceCenter[FaceNum];

//            aGrid->Cells[PCell].FaceCoincide[FaceNum] =
//                                (NgbCell>=0)
//                                ?GetCoincidance(aGrid->Cells[PCell].CellCenter,
//                                                aGrid->Cells[NgbCell].CellCenter,
//                                                aGrid->Cells[PCell].FaceCenter[FaceNum],
//                                                Normalized(aGrid->Cells[PCell].Area[FaceNum]))
//                                :(aGrid->Cells[PCell].FaceCenter[FaceNum]);

            aGrid->Cells[PCell].ProjectedCenter[FaceNum] =
                                GetProjectedCenter(aGrid->Cells[PCell].CellCenter,
                                aGrid->Cells[PCell].FaceCenter[FaceNum],
                                Normalized(aGrid->Cells[PCell].Area[FaceNum]));

        }
}

void SaveNeighbours(struct Grid *Grids,char *FileName)
{

    FILE *fpMesh;
    int i,j;
    if( (fpMesh=fopen(FileName,"w"))==NULL )
    {
        puts("ElemConnect.ebr reading error.................Exiting (ReadElementConnect 000)");
        exit(-1);
    }

    for (i=0;i<Grids->NC;i++)
    {
        for (j=0;j<Grids->Cells[i].NoOfFace;j++)
            fprintf(fpMesh," %8d",Grids->Cells[i].Neighbor[j]);
        fprintf(fpMesh,"\n");
    }
    fclose(fpMesh);
}


int isContain(int a,int *aList,int s)
{
    int i;
    for (i=0;i<s;i++)
        if (a==aList[i])
            return 1;
    return 0;
}

int isMathedFaces(int *FaceNodeList_1, int NumOfNodes_1,
                  int *FaceNodeList_2, int NumOfNodes_2)
{
//  first rapid cheking
    if (NumOfNodes_1!=NumOfNodes_2) return 0;
//  second rapid cheking
    int listSize;
    int totalNodeNum_1=0,totalNodeNum_2=0,k;
    listSize=NumOfNodes_1;
    for (k=0;k<listSize;k++)
    {
        totalNodeNum_1+=FaceNodeList_1[k];
        totalNodeNum_2+=FaceNodeList_2[k];
    }
    if (totalNodeNum_1!=totalNodeNum_2) return 0;
//  final checking
    for (k=0;k<listSize;k++)
    {
        if (!isContain(FaceNodeList_1[k],FaceNodeList_2,listSize))
            return 0;
    }
    return 1;
}

void FindFillNeighbours(struct Grid *Grids)
{
    int i,j,ii,jj,k;
    int counter=0;

    for (i=0;i<Grids->NC;i++)
        for (j=0;j<Grids->Cells[i].NoOfFace;j++)
            Grids->Cells[i].Neighbor[j]=-1;

    int NumOfNodes_1,NumOfNodes_2,
        FaceNodeList_1[4],FaceNodeList_2[4];

    for (i=0;i<Grids->NC;i++)
    {
        for (j=0;j<Grids->Cells[i].NoOfFace;j++)
            if (Grids->Cells[i].Neighbor[j]==-1)
            {
                switch (Grids->RCells[i].Type)
                {
                    case 4:
                        NumOfNodes_1=BREAK_NODES[j][0];
                        for (k=0;k<NumOfNodes_1;k++)
                            FaceNodeList_1[k]=Grids->RCells[i].NodeList[BREAK_NODES[j][k+1]];
                    break;
                    case 6:
                        NumOfNodes_1=TETRA_NODES[j][0];
                        for (k=0;k<NumOfNodes_1;k++)
                            FaceNodeList_1[k]=Grids->RCells[i].NodeList[TETRA_NODES[j][k+1]];
                    break;
                    case 5:
                        NumOfNodes_1=WEDGE_NODES[j][0];
                        for (k=0;k<NumOfNodes_1;k++)
                            FaceNodeList_1[k]=Grids->RCells[i].NodeList[WEDGE_NODES[j][k+1]];
                    break;
                    case 7:
                        NumOfNodes_1=PYRAMID_NODES[j][0];
                        for (k=0;k<NumOfNodes_1;k++)
                            FaceNodeList_1[k]=Grids->RCells[i].NodeList[PYRAMID_NODES[j][k+1]];
                    break;
                }
                for (ii=0;ii<Grids->NC;ii++)
                    if (i!=ii)
                    {
                        int flag=1;
                        for (k=0;k<NumOfNodes_1;k++)
                            if (!isContain(FaceNodeList_1[k],Grids->RCells[ii].NodeList,Grids->RCells[ii].NodePerCell))
                                {
                                    flag=0;
                                    break;
                                }
                        if (flag==1) Grids->Cells[i ].Neighbor[j ]=ii;
                    }
            }
        printf("\r\t");
        printf("%6.3f Percent Completed...",((double)(i+1)/(double)Grids->NC)*100.);
    }
    printf("\n");
}

void FillNeighbours(struct Grid *Grids,char *FileName)
{
    FILE *fpMesh;
    int i,j;
    char dummystr[80];
    if( (fpMesh=fopen(FileName,"r"))==NULL )
    {
        puts("ElemConnect.ebr reading error.................Exiting (ReadElementConnect 000)");
        exit(-1);
    }
    for (i=0;i<Grids->NC;i++)
        for (j=0;j<Grids->Cells[i].NoOfFace;j++)
        {
            GetToken(dummystr,fpMesh);
            Grids->Cells[i].Neighbor[j]=atoi(dummystr);
        }
    fclose(fpMesh);
}

void CheckNodes(struct Grid *Grids, char Dim, double error)
{
    int counterx,countery,counterz,i;
    counterx=countery=counterz=0;
    double dif;
    struct MyVector n0,n1,n2,n3,n4,n5;

    switch (Dimension)
    {
        case '2':
            {
//                for (i=0;i<Grids->NN;i++)
//                {
//                    if (fabs(Grids->Nodes[i].Pos.X)<error)
//                    {
//                        Grids->Nodes[i].Pos.X=0.0;
//                        counterx++;
//                    }
//                    if (fabs(Grids->Nodes[i].Pos.Y)<error)
//                    {
//                        Grids->Nodes[i].Pos.Y=0.0;
//                        countery++;
//                    }
//                    if (fabs(Grids->Nodes[i].Pos.Z)<error)
//                    {
//                        Grids->Nodes[i].Pos.Z=0.0;
//                        counterz++;
//                    }
//                }

                for(i=0;i<Grids->NC;i++)
                {
                    switch(Grids->RCells[i].Type)
                    {
                        case 5: ///wedge
                            {
                                Grids->Nodes[Grids->RCells[i].NodeList[0]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[3]].Pos.X;
                                Grids->Nodes[Grids->RCells[i].NodeList[0]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[3]].Pos.Z;

                                Grids->Nodes[Grids->RCells[i].NodeList[2]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[5]].Pos.X;
                                Grids->Nodes[Grids->RCells[i].NodeList[2]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[5]].Pos.Z;

                                Grids->Nodes[Grids->RCells[i].NodeList[1]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[4]].Pos.X;
                                Grids->Nodes[Grids->RCells[i].NodeList[1]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[4]].Pos.Z;
                            }break;
                        case 4: ///Brick
                            {
//                                Grids->Nodes[Grids->RCells[i].NodeList[0]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[3]].Pos.X;
//                                Grids->Nodes[Grids->RCells[i].NodeList[0]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[3]].Pos.Z;
//
//                                Grids->Nodes[Grids->RCells[i].NodeList[2]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[5]].Pos.X;
//                                Grids->Nodes[Grids->RCells[i].NodeList[2]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[5]].Pos.Z;
//
//                                Grids->Nodes[Grids->RCells[i].NodeList[1]].Pos.X=Grids->Nodes[Grids->RCells[i].NodeList[4]].Pos.X;
//                                Grids->Nodes[Grids->RCells[i].NodeList[1]].Pos.Z=Grids->Nodes[Grids->RCells[i].NodeList[4]].Pos.Z;
                            }break;
                        default:
                            {
                                printf("error\n");
                                getchar();
                            }break;
                    }
                }


            }break;
        case '3':
            {
                printf("3D !!!\n");
                getchar();
            }break;
    }
    printf("Mesh is %cD and %d Nodes fixed successfuly in X\n",Dim,counterx);
    printf("Mesh is %cD and %d Nodes fixed successfuly in Y\n",Dim,countery);
    printf("Mesh is %cD and %d Nodes fixed successfuly in Z\n\n",Dim,counterz);
//    printf("Press Enter to Continue\n");
//    getchar();
}

void CheckMesh(struct Grid *Grids, double TError)
{
    int i;
    double Cellcenter_Error;
    switch (Dimension)
    {
        case '2':
            {
                for(i=0;i<Grids->NC;i++)
                {
                    Cellcenter_Error=fabs(Grids->Cells[i].CellCenter.Y-0.5);
                    if (Cellcenter_Error>TError)
                    {
                        printf("\nAfter Node Modification, there are still some Cellcenter ERROR in 2D mesh\n");
                        printf("it highly reccomanded NOT to continue before remeshing\n");
                        getchar();
                    }
                }
            }break;
        default:
        {

        }break;
    }
}

void FillGrid(struct Grid *Grids, char *FileName)
{
    char NUMNP[10+1]
        ,NELEM[10+1]
        ,NGRPS[10+1]
        ,NBSETS[10+1]
        ,NDFCD[10+1]
        ,NDFVL[10+1]

        ,ElemNum[8+1]
        ,ElemTyp[3+1]
        ,ElemVerNum[3+1]

        ,VerNum[8+1]

        ,dummyLine[80+1];
    int  n,i,NDCSection=0,iii;
    char path[80];
    FILE *NFF,*NodeCoord,*ElemConn,*Bdry;

    Grids->NC=Grids->NN=Grids->NB=0;

    sprintf(path,"./mesh/%s",FileName);
    if( (NFF=fopen(path,"r"))==NULL )
    {
        puts("mymesh.neu reading error.................Exiting ( ebrsMaker 000)");
        exit(-1);
    }
    while (feof(NFF)==0)
    {
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        NDCSection=1;
        fgets(NUMNP,10+1,NFF);
        fgets(NELEM,10+1,NFF);
        fgets(NGRPS,10+1,NFF);
        fgets(NBSETS,10+1,NFF);
        fgets(NDFCD,10+1,NFF);
        fgets(NDFVL,10+1,NFF);

        (*Grids).NN=atoi(NUMNP);
        (*Grids).NC=atoi(NELEM);

        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);

        //  Write Node Coordinate file
        sprintf(path,"./ebrs/%s.NodeCoord.ebr",FileName);
        if( (NodeCoord=fopen(path,"w"))==NULL )
        {
            puts("NodeCoord.ebr writing error.................Exiting ( ebrsMaker 000)");
            exit(-1);
        }
        fprintf(NodeCoord,"%s\n",NUMNP);
        for (i=0;i<(*Grids).NN;i++)
        {
            fgets(dummyLine,80,NFF);
            fprintf(NodeCoord,"%s",dummyLine);
        }
        fclose(NodeCoord);

        fgets(dummyLine,80,NFF);
        fgets(dummyLine,80,NFF);

        //  Write Element Connectivity file
        sprintf(path,"./ebrs/%s.Connectivity.ebr",FileName);
        if( (ElemConn=fopen(path,"w"))==NULL )
        {
            puts("ElemConnect.ebr writing error.................Exiting ( ebrsMaker 000)");
            exit(-1);
        }
        fprintf(ElemConn,"%d\n",atoi(NELEM));
        int NE;
        NE=atoi(NELEM);
        for (n=0;n<NE;n++)
        {
            GetToken(ElemNum,NFF);
            GetToken(ElemTyp,NFF);
            fprintf(ElemConn,"%d ",atoi(ElemTyp));
            GetToken(ElemVerNum,NFF);
            fprintf(ElemConn,"%d ",atoi(ElemVerNum));
            int nn,EVN;
            EVN=atoi(ElemVerNum);
            for (nn=0;nn<EVN;nn++)
            {
                GetToken(dummyLine,NFF);
                fprintf(ElemConn,"%d ",atoi(dummyLine));
            }
            fprintf(ElemConn,"\n");
        }
        fclose(ElemConn);
        do
            GetToken(dummyLine,NFF);
        while (strncmp("BOUNDARY",dummyLine,26)!=0);
        fgets(dummyLine,80,NFF);

        //  Write Element Boundary Conditions file
        sprintf(path,"./ebrs/%s.Boundary.ebr",FileName);
        if( (Bdry=fopen(path,"w"))==NULL )
        {
            puts("Boundary.ebr writing error.................Exiting ( ebrsMaker 000)");
            exit(-1);
        }

        FILE *BdrySet;
        sprintf(path,"./mesh/%s.bdry.set",FileName);
        if( (BdrySet=fopen(path,"r"))==NULL )
        {
            puts("No Boundary Conditions Set File ................................!!!");
            exit(-1);
        }

        char NB[8+1],Sdummy[80],ElemNum2[10+1],BdryType[5+1],FaceNum[5+1];

        GetToken(Sdummy,BdrySet);
        int NoOfBSets=0;

        NoOfBSets=atoi(Sdummy);
        printf("%s      %d\n",path,NoOfBSets);

        int kk;
        kk=0;
        Grids->NRBG=0;
        Grids->RBNo=(double *)malloc(sizeof(double));

        for(iii=0;iii<NoOfBSets;iii++)
        {
            GetToken(Sdummy,NFF);
            GetToken(Sdummy,NFF);
            GetToken(NB,NFF);   //second NO after boundary name in NEU (e.g second NO after XP)
            GetToken(Sdummy,NFF);
            GetToken(Sdummy,NFF);
            Grids->NB+=atoi(NB);
            GetToken(Sdummy,BdrySet);

            if(atoi(Sdummy)==1)
            {
                Grids->NRBG+=1;
                Grids->RBNo=(double *)realloc(Grids->RBNo,sizeof(double)*(Grids->NRBG+1));
                Grids->RBNo[Grids->NRBG-1]=0;
//                printf("RB found !\n");
//                getchar();
            }

            for (n=0;n<atoi(NB);n++)
            {
                fgets(ElemNum2,10+1,NFF);
                fprintf(Bdry,"%s",ElemNum2);
                fgets(BdryType,5+1,NFF);
                //    BdryType
                fprintf(Bdry,"%s%d","    ",atoi(Sdummy));

                if(atoi(Sdummy)==1)
                    Grids->RBNo[Grids->NRBG-1]+=1;

                //    BdryFaceNo
                fgets(FaceNum,5+1,NFF);
                fprintf(Bdry,"%s",FaceNum);
                //    BdryValue
                fprintf(Bdry,"%20.10f\n",0.0);
                fgets(dummyLine,80,NFF);
            }
            GetToken(Sdummy,NFF);
            GetToken(Sdummy,NFF);
            GetToken(Sdummy,NFF);
            GetToken(Sdummy,NFF);
        }
        fclose(Bdry);
        fclose(BdrySet);
        break;
    }
    fclose(NFF);

    char dummyPath[80];
    sprintf(dummyPath,"./ebrs/%s.NodeCoord.ebr",FileName);
    Grids->Nodes=(struct Node*) malloc(sizeof(struct Node)*Grids->NN);
    FillNodes(Grids,dummyPath);
    Grids->Cells=(struct Cell*) malloc(sizeof(struct Cell)*Grids->NC);
    Grids->RCells=(struct RawCell*) malloc(sizeof(struct RawCell)*Grids->NC);
    sprintf(dummyPath,"./ebrs/%s.Connectivity.ebr",FileName);
    FillConnectivity(Grids,dummyPath);

//    CheckNodes(Grids,Dimension,Mesh_Truncation);

    CompCellGeoPara(Grids);
/// Comment three following lines if the NGBs were written ///
    FindFillNeighbours(Grids);
    sprintf(dummyPath,"./ebrs/%s.Neighbour.ebr",FileName);
    SaveNeighbours(Grids,dummyPath);

    sprintf(dummyPath,"./ebrs/%s.Neighbour.ebr",FileName);
    FillNeighbours(Grids,dummyPath);
///////////////////////////////////////////

    CompleteFacePara(Grids);
    Grids->BoundaryCells=(struct BoundaryCell*) malloc(sizeof(struct BoundaryCell)*Grids->NB);
    sprintf(dummyPath,"./ebrs/%s.Boundary.ebr",FileName);
    FillBoundary(Grids,dummyPath);
    FillRigidBoundary(Grids);
    FillNodeSharing(Grids);
    MallocFlowVariables(Grids);
    InitFlowVariables(Grids);
    InitDynaVariables(Grids);

//    CheckMesh(Grids,Mesh_Truncation);
}

void CellProp(struct Grid *Grids, int PCell)
{
    int i,j;

    printf("the Properties of cell No: %d is:\n",PCell);
    printf("the No of faces: %d\n",Grids->Cells[PCell].NoOfFace);
    switch (Grids->RCells[PCell].Type)
    {
        case 4:     ///Brick
        {
//            printf("the Cell is Brick and has 6 Faces.\n");
//            printf("the cell center is: %f ; %f ; %f\n",
//                   Grids->Cells[PCell].CellCenter.X,
//                   Grids->Cells[PCell].CellCenter.Y,
//                   Grids->Cells[PCell].CellCenter.Z);
//            printf("\n");
            for (i=0;i<6;i++)
            {
                printf("the face No %d has the following nodes:\n",i);
                for (j=0;j<BREAK_NODES[i][0];j++)
                {
//                    printf("node No: %d \n",
//                           Grids->RCells[PCell].NodeList[BREAK_NODES[i][j+1]]+1);
                    printf("node position: %f  ; %f  ;  %f \n",
                           Grids->Nodes[Grids->RCells[PCell].NodeList[BREAK_NODES[i][j+1]]].Pos.X,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[BREAK_NODES[i][j+1]]].Pos.Y,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[BREAK_NODES[i][j+1]]].Pos.Z);
                }
                printf("Area face no %d is: %f  ; %f  ;  %f \n",i,
                       Grids->Cells[PCell].Area[i].X,
                       Grids->Cells[PCell].Area[i].Y,
                       Grids->Cells[PCell].Area[i].Z);
//                printf("Center of face no %d is: %f  ; %f  ;  %f \n",i,
//                       Grids->Cells[PCell].FaceCenter[i].X,
//                       Grids->Cells[PCell].FaceCenter[i].Y,
//                       Grids->Cells[PCell].FaceCenter[i].Z);
                printf("\n");
            }

        }break;
        case 5:     ///Wedge
        {
//            printf("the Cell is Wedge and has 5 Faces.\n");
//            printf("the cell center is: %f ; %f ; %f\n",
//                   Grids->Cells[PCell].CellCenter.X,
//                   Grids->Cells[PCell].CellCenter.Y,
//                   Grids->Cells[PCell].CellCenter.Z);
//            printf("\n");
            for (i=0;i<5;i++)
            {
                printf("the face No %d has the following nodes:\n",i);
                for (j=0;j<WEDGE_NODES[i][0];j++)
                {
//                    printf("node No: %d \n",
//                           Grids->RCells[PCell].NodeList[WEDGE_NODES[i][j+1]]+1);
                    printf("node position: %f  ; %f  ;  %f \n",
                           Grids->Nodes[Grids->RCells[PCell].NodeList[WEDGE_NODES[i][j+1]]].Pos.X,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[WEDGE_NODES[i][j+1]]].Pos.Y,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[WEDGE_NODES[i][j+1]]].Pos.Z);
                }
                printf("Area face no %d is: %f  ; %f  ;  %f \n",i,
                       Grids->Cells[PCell].Area[i].X,
                       Grids->Cells[PCell].Area[i].Y,
                       Grids->Cells[PCell].Area[i].Z);
                printf("\n");
            }
        }break;
        case 6:     ///Tetra
        {
            printf("the Cell is Tetra and has 4 Faces.\n");
            printf("the cell center is: %f ; %f ; %f\n",
                   Grids->Cells[PCell].CellCenter.X,
                   Grids->Cells[PCell].CellCenter.Y,
                   Grids->Cells[PCell].CellCenter.Z);
            printf("\n");
            for (i=0;i<4;i++)
            {
                printf("the face No %d has the following nodes:\n",i);
                for (j=0;j<TETRA_NODES[i][0];j++)
                {
//                    printf("node No: %d \n",
//                           Grids->RCells[PCell].NodeList[TETRA_NODES[i][j+1]]+1);
                    printf("node position: %f  ; %f  ;  %f \n",
                           Grids->Nodes[Grids->RCells[PCell].NodeList[TETRA_NODES[i][j+1]]].Pos.X,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[TETRA_NODES[i][j+1]]].Pos.Y,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[TETRA_NODES[i][j+1]]].Pos.Z);
                }
                printf("Area face no %d is: %f  ; %f  ;  %f \n",i,
                       Grids->Cells[PCell].Area[i].X,
                       Grids->Cells[PCell].Area[i].Y,
                       Grids->Cells[PCell].Area[i].Z);
                printf("\n");
            }
        }break;
        case 7:     ///Pyramid
        {
            printf("the Cell is Pyramid and has 5 Faces.\n");
            printf("the cell center is: %f ; %f ; %f\n",
                   Grids->Cells[PCell].CellCenter.X,
                   Grids->Cells[PCell].CellCenter.Y,
                   Grids->Cells[PCell].CellCenter.Z);
            printf("\n");
            for (i=0;i<5;i++)
            {
                printf("the face No %d has the following nodes:\n",i);
                for (j=0;j<PYRAMID_NODES[i][0];j++)
                {
//                    printf("node No: %d \n",
//                           Grids->RCells[PCell].NodeList[PYRAMID_NODES[i][j+1]]+1);
                    printf("node position: %f  ; %f  ;  %f \n",
                           Grids->Nodes[Grids->RCells[PCell].NodeList[PYRAMID_NODES[i][j+1]]].Pos.X,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[PYRAMID_NODES[i][j+1]]].Pos.Y,
                           Grids->Nodes[Grids->RCells[PCell].NodeList[PYRAMID_NODES[i][j+1]]].Pos.Z);
                }
                printf("Area face no %d is: %f  ; %f  ;  %f \n",i,
                       Grids->Cells[PCell].Area[i].X,
                       Grids->Cells[PCell].Area[i].Y,
                       Grids->Cells[PCell].Area[i].Z);
                printf("\n");
            }break;
        }
    }
}
