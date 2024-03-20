//void UpdatePhiLSEGradient(struct Grid *Grids)
//{
//    int PCell;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//        Grids->Variables.PhiLSEGrad[PCell]=LSEGradient(Grids->Variables.NPhi,PCell,Grids);
//}

void SetLSBC(struct MatrixCoefficient *LMC,
             double *LST,
             struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector Cell2Face,FaceCenter,V,W,Area;
    double Flux,Value,PhiF;
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        Value   = Grids->BoundaryCells[k].Value;
        FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
        Cell2Face=Minus(FaceCenter,Grids->Cells[PCell].CellCenter);
        PhiF=Grids->Variables.Phi[PCell]+Dot(Grids->Variables.PhiLSEGrad[PCell],Cell2Face);
        Area=Grids->Cells[PCell].Area[FaceNum];
        W=RigidBodyVel(Grids,FaceCenter);
        Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(W,Area);
        LST[PCell]+=-Flux*PhiF;
    }
}


void MakeLSMatrixCoeff(struct MatrixCoefficient *LMC,double *LST,struct Grid *Grids,
                       char TemporalTermDiscretizationMethodForLS,
                       char TimeDiscretizationMethodForLS)
{
    int PCell,s;
    double *NH,*H,*PH
          ,Omega0=1.0
          ,Omega1=1.0
          ,Omega2=0.0
          ,CGamma0=1.0;
    double Vol,ac[8];
    struct MyVector *LSEGrad;

    switch (TemporalTermDiscretizationMethodForLS)
    {
        case 'A' : {Omega0 =     1.0;Omega1 =     1.0 ; Omega2 =     0.0;}break;
        case 'B' : {Omega0 = 3.0/2.0;Omega1 = 4.0/2.0 ; Omega2 = 1.0/2.0;}break;
    }

    switch (TimeDiscretizationMethodForLS)
    {
        case 'I' : CGamma0 = 0.0; break;
        case 'C' : CGamma0 = 0.5; break;
        case 'E' : CGamma0 = 1.0; break;
    }

    NH=Grids->Variables.NPhi;
    PH=Grids->Variables.PPhi;
    H=Grids->Variables.Phi;
    LSEGrad=Grids->Variables.PhiLSEGrad;

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        int NoOfFace;
        Vol = Grids->Cells[PCell].Volume;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        SetSVect(ac,8,0.0);
        if (TimeDiscretizationMethodForLS=='E')
            FindConvCoeffLagged(ac,Grids,PCell,H,ConvectionMethod);
        else
            FindConvCoeff_AxNode(ac,Grids,PCell,NH,LSEGrad,'C'/*ConvectionMethod*/);

        LMC->Elem[PCell][0]=(Vol/dt*Omega0+ac[0]*(1.0-CGamma0));
        LST[PCell]=(H[PCell]*Omega1-PH[PCell]*Omega2)*Vol/dt;

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                LMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0);

        LST[PCell]+=-(ac[0]*CGamma0)*H[PCell]+(ac[NoOfFace+1]);
        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                LST[PCell]+=-((ac[s+1]*CGamma0)*H[Grids->Cells[PCell].Neighbor[s]]);
    }
}


void MallocLSMatrixes(struct Grid *Grids
                      ,QMatrix *LaspackLMC
                      ,struct MatrixCoefficient *LMC
                      ,double **MyLST
                      ,Vector *LST
                      ,Vector *L)
{
    int NC=0,NCol=0,PCell;
    NC=Grids->NC;

    (*LMC).Elem=(double**) malloc(sizeof(double*)*NC);
    for(PCell=0;PCell<NC;PCell++)
    {

        NCol=Grids->Cells[PCell].NoOfFace+1;
        (*LMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
    }
    (*MyLST)  =(double*) malloc(sizeof(double)*NC);

    Q_Constr(LaspackLMC,"GLMC",NC,False,Rowws,Normal,True);
    V_Constr(LST,"LST",NC,Normal,True);
    V_Constr(L,"L",NC,Normal,True);
}

void FreeLSMatrixes( struct Grid *Grids
                    ,QMatrix *LaspackLMC
                    ,struct MatrixCoefficient *LMC
                    ,double **MyLST
                    ,Vector *LST
                    ,Vector *L)
{
    int NC=0,PCell;
    NC=Grids->NC;
    for(PCell=0;PCell<NC;PCell++)
        free((*LMC).Elem[PCell]);
    free((*MyLST));

    Q_Destr(LaspackLMC);
    V_Destr(LST);
    V_Destr(L);
}

void LSSolver(struct Grid *Grids
             ,char TemporalTermDiscretizationMethodForLS
             ,char TimeDiscretizationMethodForLS)
{
    struct MatrixCoefficient LMC;
    double *MyLST;
    int i;
    QMatrix LaspackLMC;
    Vector LST;
    Vector L;

//    Memory allocation
    MallocVoFMatrixes(Grids, &LaspackLMC, &LMC, &MyLST, &LST, &L);
//    printf("VoF started ...\n");
    int dInt=0;
    do
    {
        dInt=0;
        UpdatePhiLSEGradient(Grids);
        MakeLSMatrixCoeff(&LMC,MyLST,Grids,TemporalTermDiscretizationMethodForLS,TimeDiscretizationMethodForLS);
        SetLSBC(&LMC,MyLST,Grids);
        CopyToLaspackVector(&L,Grids->Variables.NPhi,Grids);
        CopyToLaspackVector(&LST,MyLST,Grids);
        CopyToLaspackMatrix(&LaspackLMC,&LMC,Grids);
        SetRTCAccuracy(1e-6);
        L=*BiCGSTABIter(&LaspackLMC,&L,&LST,100,ILUPrecond,1.0);
        dInt+=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NPhi,&L,Grids);
        printf("########################  LS Solver Iterations %d\n",dInt);
    }
    while(dInt>0);
    FreeLSMatrixes(Grids, &LaspackLMC, &LMC, &MyLST, &LST, &L);
}

