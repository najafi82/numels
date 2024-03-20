void MallocKOMatrixes(struct Grid *Grids
                   ,QMatrix *LaspackMC
                   ,struct MatrixCoefficient *OmMC
                   ,struct MatrixCoefficient *KMC
                   ,double **OmST
                   ,double **KST
                   ,Vector *ST
                   ,Vector *Fee)
{
    int NC=0,NCol=0,PCell;
    NC=Grids->NC;

    (*OmMC).Elem =(double**) malloc(sizeof(double*)*NC);
    (*KMC).Elem =(double**) malloc(sizeof(double*)*NC);

    for(PCell=0;PCell<NC;PCell++)
    {
        NCol=Grids->Cells[PCell].NoOfFace+1;
        (*OmMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
        (*KMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
    }

    *OmST=(double*) malloc(sizeof(double)*NC);
    *KST=(double*) malloc(sizeof(double)*NC);

    Q_Constr(LaspackMC,"LaspackMC",NC,False,Rowws,Normal,True);
    V_Constr(ST,"ST",NC,Normal,True);
    V_Constr(Fee,"Fee",NC,Normal,True);
}

void FreeKOMatrixes(struct Grid *Grids
                   ,QMatrix *LaspackMC
                   ,struct MatrixCoefficient *OmSMC
                   ,struct MatrixCoefficient *KSMC
                   ,double **OmST
                   ,double **KST
                   ,Vector *ST
                   ,Vector *Fee)
{
    int NC=0,PCell;
    NC=Grids->NC;
    for(PCell=0;PCell<NC;PCell++)
    {
        free((*OmSMC).Elem[PCell]);
        free((*KSMC).Elem[PCell]);
    }

    free((*OmSMC).Elem);
    free((*KSMC).Elem);

    free(*OmST);
    free(*KST);

    Q_Destr(LaspackMC);
    V_Destr(ST);
    V_Destr(Fee);
}

void UpdateKOmLSEGrad(struct Grid *Grids)
{
    int PCell;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.OmLSEGrad[PCell]=LSEGradient(Grids->Variables.NOm,PCell,Grids);
        Grids->Variables.KLSEGrad[PCell]=LSEGradient(Grids->Variables.NK,PCell,Grids);
    }
}

void GetKoNormalGradParts(double *aP,double *aN,double *ST,
                        double *Fee,struct MyVector *LSEGrad,
                        int PCell,int FaceNum,struct Grid *Grids)
{
    struct MyVector PPrime,NPrime,CPVector,CNVector;
    double delta;
    int NgbCell,i,MatchedFaceNum;
    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*ST=0.0;

    PPrime = Grids->Cells[PCell].ProjectedCenter[FaceNum];
    CPVector = Minus(PPrime,Grids->Cells[PCell].CellCenter);
    if (NgbCell>=0)
    {
        for(i=0;i<Grids->Cells[NgbCell].NoOfFace;i++)
            if (Mag(Sum(Grids->Cells[  PCell].Area[FaceNum]
                       ,Grids->Cells[NgbCell].Area[      i]))<geps)
            {
                MatchedFaceNum=i;
                break;
            }
        NPrime = Grids->Cells[NgbCell].ProjectedCenter[MatchedFaceNum];
//        NPrime = Grids->Cells[NgbCell].ProjectedCenter[FaceNum];
        CNVector = Minus(NPrime,Grids->Cells[NgbCell].CellCenter);
        delta = Mag(Minus(PPrime,NPrime));
        *aP = -1.0/delta;
        *aN = +1.0/delta;
        *ST = (1.0/delta)*(Dot(LSEGrad[NgbCell],CNVector)-Dot(LSEGrad[PCell],CPVector));
    }
}

void FindKoDiffCoeff_AxNode(double *a1,struct Grid *Grids,int PCell)
{
    int i;
    int NoOfFace;
    double ViscF,DensF,NuF,AMag,aP,aN,ST,Nut;
    double w1,w2;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
        if (Grids->Cells[PCell].Neighbor[i]>=0)
        {
            GetKoNormalGradParts(&aP,&aN,&ST,Grids->Variables.NK,Grids->Variables.KLSEGrad,PCell,i,Grids);
            AMag=Mag(Grids->Cells[PCell].Area[i]);
            ViscF=0.5*(ViscCal(Grids->Variables.PAlpha[PCell])+ViscCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));
            DensF=0.5*(DensCal(Grids->Variables.PAlpha[PCell])+DensCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));

            ///test
            w1=Mag(Minus(Grids->Cells[PCell].CellCenter
                 ,Grids->Cells[PCell].FaceCenter[i]));
            w2=Mag(Minus(Grids->Cells[Grids->Cells[PCell].Neighbor[i]].CellCenter
                 ,Grids->Cells[PCell].FaceCenter[i]));
            Nut=(w2*(Grids->Variables.NuTurb[PCell])
                     +w1*(Grids->Variables.NuTurb[Grids->Cells[PCell].Neighbor[i]]))/(w1+w2);

//            Nut=0.5*((Grids->Variables.NuTurb[PCell])+(Grids->Variables.NuTurb[Grids->Cells[PCell].Neighbor[i]]));

            NuF=(ViscF/DensF)+(Nut*sigstar);

            a1[0]  += NuF*aP*AMag;
            a1[i+1] = NuF*aN*AMag;
            a1[NoOfFace+1]  += NuF*ST*AMag;
        }
}


void FindOmDiffCoeff_AxNode(double *a1,struct Grid *Grids,int PCell)
{
    int i;
    int NoOfFace;
    double ViscF,DensF,NuF,AMag,aP,aN,ST,Nut;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
        if (Grids->Cells[PCell].Neighbor[i]>=0)
        {
            GetKoNormalGradParts(&aP,&aN,&ST,Grids->Variables.NOm,Grids->Variables.OmLSEGrad,PCell,i,Grids);
            AMag=Mag(Grids->Cells[PCell].Area[i]);
            ViscF=0.5*(ViscCal(Grids->Variables.PAlpha[PCell])+ViscCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));
            DensF=0.5*(DensCal(Grids->Variables.PAlpha[PCell])+DensCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));
            Nut=0.5*((Grids->Variables.NuTurb[PCell])+(Grids->Variables.NuTurb[Grids->Cells[PCell].Neighbor[i]]));
            NuF=(ViscF/DensF)+(Nut*sig);
            a1[0]  += NuF*aP*AMag;
            a1[i+1] = NuF*aN*AMag;
            a1[NoOfFace+1]  += NuF*ST*AMag;
        }
}


void FindKOConvCoeff_AxNode(double *a2,struct Grid *Grids,int PCell,double *H,struct MyVector *LSEGrad,char Method)
{
    double Flux,G=0.7,aP,aN,ST,kapa=0.1,betaf,AU,AN;
    struct MyVector Area,FaceCenter,W,r;
    int s,NoOfFace,ACell,DCell;
    struct DonnerAcceptor DA;

    NoOfFace=Grids->Cells[PCell].NoOfFace;
    a2[NoOfFace+1]=0.0;
    DonnAccCell(&DA,Grids,PCell);

    for (s=0;s<NoOfFace;s++)
    {
        if (Grids->Cells[PCell].Neighbor[s]>=0)
        {
            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
            Area=Grids->Cells[PCell].Area[s];
            W=RigidBodyVel(Grids,FaceCenter);
            Flux=Grids->Variables.NSF[PCell].Face[s]-Dot(W,Area);
      //Cell Matrix Coefficients Calculation with different interpolation schemes
            switch (Method)
            {
                case 'C': // Linear Interpolation (CDS)
                {
                    GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,s,Grids);
                    a2[0]=a2[0]+aP*Flux;
                    a2[s+1]=aN*Flux;
                    a2[NoOfFace+1]+= -ST*Flux;
                }break;
                case 'T': // Interpolation (UDS+CDS)
                {
                    a2[0]=a2[0]+G*0.5*Flux+(1.0-G)*max(Flux,0.0);
                    a2[s+1]=G*0.5*Flux+(1.0-G)*(-max(-Flux,0.0));
                    a2[NoOfFace+1]=0.0;
                }break;

                case 'U': // Upwind Interpolation (UDS)
                {
                    a2[0]=a2[0]+max(Flux,0.0);
                    a2[s+1]=-max(-Flux,0.0);
                    a2[NoOfFace+1]=0.0;
                }break;

                case 'D': // Deffered CDS and UDS Interpolation
                {
                    GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,s,Grids);
                    a2[0]=a2[0]+max(Flux,0.0);
                    a2[s+1]=-max(-Flux,0.0);
                    a2[NoOfFace+1]=+(max(Flux,0.0)*H[PCell]-max(-Flux,0.0)*H[Grids->Cells[PCell].Neighbor[s]])
                                   -(aP*H[PCell]+aN*H[Grids->Cells[PCell].Neighbor[s]]+ST)*Flux;
                }break;
                case 'G': // Gamma Differencing Scheme (Jassac Method)
                {
                    DCell=DA.Face[s].d;
                    ACell=DA.Face[s].a;
                    r=Sum(Grids->Cells[ACell].CellCenter,MinusVec(Grids->Cells[DCell].CellCenter));
                    AU=H[ACell]-2.0*Dot(LSEGrad[DCell],r);
                    AN=(H[DCell]-AU)/(H[ACell]-AU);
                    if (isequal(H[ACell],AU))
                    {
                        AN=100.0;
                        if (isequal(H[DCell],AU))
                            AN=1.0;
                    }

                    if ((AN>=kapa)&&(AN<1.0)) betaf=0.5;//CDS
                    else if ((AN>0.0)&&(AN<kapa)) betaf=1.0-(AN/(2.0*kapa));//Hybrid
                    else betaf=1.0;//UDS

                    if (PCell==DCell)
                    {
                        a2[0]+=betaf*Flux;
                        a2[s+1]=(1.0-betaf)*Flux;
                    }
                    else
                    {
                        a2[0]+=(1.0-betaf)*Flux;
                        a2[s+1]=betaf*Flux;
                    }
                }break;
            }
        }
    }
}

void MakeKoMatrix(struct MatrixCoefficient *SMC,double *KST,struct Grid *Grids,char ConvectionMethod
                 ,char TimeDiscretizationMethodForConvection
                 ,char TimeDiscretizationMethodForDiffusion)
{
    int PCell,s;
    double *NH,*H
          ,Omega0=1.0
          ,DGamma0=1.0
          ,CGamma0=1.0;
    double ad[8]
          ,ac[8];
    struct MyVector *LSEGrad;

    double dummy;

    switch (TimeDiscretizationMethodForConvection)
    {
        case 'I' : CGamma0 = 0.0; break;
        case 'C' : CGamma0 = 0.5; break;
        case 'E' : CGamma0 = 1.0; break;
    }

    switch (TimeDiscretizationMethodForDiffusion)
    {
        case 'I' : DGamma0 = 0.0; break;
        case 'C' : DGamma0 = 0.5; break;
        case 'E' : DGamma0 = 1.0; break;
    }

    NH=Grids->Variables.NK;
    H=Grids->Variables.K;
    LSEGrad=Grids->Variables.KLSEGrad;

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        SetSVect(ad,8,0.0);
        SetSVect(ac,8,0.0);

        FindKOConvCoeff_AxNode(ac,Grids,PCell,NH,LSEGrad,ConvectionMethod);

        FindKoDiffCoeff_AxNode(ad,Grids,PCell);

        dummy=Grids->Cells[PCell].Volume;

        SMC->Elem[PCell][0]=(Grids->Cells[PCell].Volume/dt*Omega0-ad[0]*(1.0-DGamma0)+ac[0]*(1.0-CGamma0));

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                SMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0)-ad[s+1]*(1.0-DGamma0);

        KST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*H[PCell])+(ac[NoOfFace+1]+ad[NoOfFace+1]);

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                KST[PCell]+= -((ac[s+1]*CGamma0+ad[s+1]*DGamma0)*H[Grids->Cells[PCell].Neighbor[s]]);

    }
}

void MakeOmMatrix(struct MatrixCoefficient *OmSMC,double *OmST,struct Grid *Grids,char ConvectionMethod
                 ,char TimeDiscretizationMethodForConvection
                 ,char TimeDiscretizationMethodForDiffusion)
{
    int PCell,s;
    double *NH,*H
          ,Omega0=1.0
          ,DGamma0=1.0
          ,CGamma0=1.0;
    double ad[8]
          ,ac[8];
    double ux,uy,uz,vx,vy,vz,wx,wy,wz,pt,vol1;
    struct MyVector *LSEGrad;

    switch (TimeDiscretizationMethodForConvection)
    {
        case 'I' : CGamma0 = 0.0; break;
        case 'C' : CGamma0 = 0.5; break;
        case 'E' : CGamma0 = 1.0; break;
    }

    switch (TimeDiscretizationMethodForDiffusion)
    {
        case 'I' : DGamma0 = 0.0; break;
        case 'C' : DGamma0 = 0.5; break;
        case 'E' : DGamma0 = 1.0; break;
    }

    NH=Grids->Variables.NOm;
    H=Grids->Variables.Om;
    LSEGrad=Grids->Variables.OmLSEGrad;

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        vol1=Grids->Cells[PCell].Volume;
        SetSVect(ad,8,0.0);
        SetSVect(ac,8,0.0);

        FindKOConvCoeff_AxNode(ac,Grids,PCell,NH,LSEGrad,ConvectionMethod);
        FindOmDiffCoeff_AxNode(ad,Grids,PCell);

        ux=Grids->Variables.ULSEGrad[PCell].X;
        uy=Grids->Variables.ULSEGrad[PCell].Y;
        uz=Grids->Variables.ULSEGrad[PCell].Z;

        vx=Grids->Variables.VLSEGrad[PCell].X;
        vy=Grids->Variables.VLSEGrad[PCell].Y;
        vz=Grids->Variables.VLSEGrad[PCell].Z;

        wx=Grids->Variables.WLSEGrad[PCell].X;
        wy=Grids->Variables.WLSEGrad[PCell].Y;
        wz=Grids->Variables.WLSEGrad[PCell].Z;

        pt=Grids->Variables.NuTurb[PCell]*(2.0*pow(ux,2.0)+2.0*pow(vy,2.0)+2.0*pow(wz,2.0)+pow(uy+vx,2.0)+pow(uz+wx,2.0)+pow(vz+wy,2.0));

        OmSMC->Elem[PCell][0]=(vol1/dt*Omega0-ad[0]*(1.0-DGamma0)+ac[0]*(1.0-CGamma0));

        switch (Omequation)
        {
            case 'A' : { OmSMC->Elem[PCell][0]-=(alf*pt*vol1/Grids->Variables.K[PCell]); } ;break;

            case 'B' : { OmSMC->Elem[PCell][0]+=-(alf*pt*vol1/Grids->Variables.K[PCell])
                                                +(beta*vol1*Grids->Variables.Om[PCell]); } ;break;
            case 'C' : { OmSMC->Elem[PCell][0]+=(beta*vol1*Grids->Variables.Om[PCell]); } ;break;
        }

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                OmSMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0)-ad[s+1]*(1.0-DGamma0);

        OmST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*H[PCell])+(ac[NoOfFace+1]+ad[NoOfFace+1]);

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                OmST[PCell]+= -((ac[s+1]*CGamma0+ad[s+1]*DGamma0)*H[Grids->Cells[PCell].Neighbor[s]]);
    }
}

void SetKOSourceTerms(double *KST,double *OmST,struct Grid *Grids)
{
    int PCell;
    double Vol,Dens;
    struct MyVector GradP,Vel,Added,GradK;
    double ux,uy,uz,vx,vy,vz,wx,wy,wz,pt;

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        ux=Grids->Variables.ULSEGrad[PCell].X;
        uy=Grids->Variables.ULSEGrad[PCell].Y;
        uz=Grids->Variables.ULSEGrad[PCell].Z;
        vx=Grids->Variables.VLSEGrad[PCell].X;
        vy=Grids->Variables.VLSEGrad[PCell].Y;
        vz=Grids->Variables.VLSEGrad[PCell].Z;
        wx=Grids->Variables.WLSEGrad[PCell].X;
        wy=Grids->Variables.WLSEGrad[PCell].Y;
        wz=Grids->Variables.WLSEGrad[PCell].Z;

        pt=Grids->Variables.NuTurb[PCell]*(2.0*pow(ux,2.0)+2.0*pow(vy,2.0)+2.0*pow(wz,2.0)+pow(uy+vx,2.0)+pow(uz+wx,2.0)+pow(vz+wy,2.0));

        Vol=Grids->Cells[PCell].Volume;
        KST[PCell]+=(Grids->Variables.K[PCell])*Vol/dt+pt*Vol-betastar*(Grids->Variables.Om[PCell])*(Grids->Variables.K[PCell])*Vol;

        switch (Omequation)
        {
            case 'A' : {OmST[PCell]+=(Grids->Variables.Om[PCell])*Vol/dt-beta*Vol*pow(Grids->Variables.Om[PCell],2.0);} break;
            case 'B' : {OmST[PCell]+=(Grids->Variables.Om[PCell])*Vol/dt;} break;
            case 'C' : {OmST[PCell]+=(Grids->Variables.Om[PCell])*Vol/dt+alf*pt*Vol*(Grids->Variables.Om[PCell])/(Grids->Variables.K[PCell]);} break;
        }
    }
}

void SetKOBC(struct MatrixCoefficient *KSMC
                  ,struct MatrixCoefficient *OmSMC
                  ,double *KST
                  ,double *OmST
                  ,struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector n,Area,FaceCenter,V,W,rC,Pf,PC;
    double Amag,Flux,Value,dn,nu,nuk,nuo,vmag;
    double Omegainlet,l,Kinlet;
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        Value   = Grids->BoundaryCells[k].Value;
        switch (Type)
        {
            case 1:// Rigid Body Wall Flow Function;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
                nuk=nu+(Grids->Variables.NuTurb[PCell])*sigstar;
                nuo=nu+(Grids->Variables.NuTurb[PCell])*sig;
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);

                //  Set Diffusive Flux Boundary
                KSMC->Elem[PCell][0]+=nuk/dn*Amag;
                OmSMC->Elem[PCell][0]+=nuo/dn*Amag;

                Omegainlet=60.*nu/(0.75*pow(dn,2.0));
                KST[PCell]+=nuk/dn*Amag*(0.0-Dot(Grids->Variables.KLSEGrad[PCell],rC));
                OmST[PCell]+=nuo/dn*Amag*(Omegainlet-Dot(Grids->Variables.OmLSEGrad[PCell],rC));

                //  Set Convective Flux Boundary
                OmST[PCell]+=-Flux*Omegainlet;
            }break;
            case 3:// Out-Flow;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                Pf=Minus(FaceCenter,Grids->Cells[PCell].CellCenter);
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
                PC=Minus(FaceCenter,Grids->Cells[PCell].CellCenter);
                W=RigidBodyVel(Grids,FaceCenter);
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
                nuk=nu+(Grids->Variables.NuTurb[PCell])*sigstar;
                nuo=nu+(Grids->Variables.NuTurb[PCell])*sig;
                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);;

                //  Set Diffusive Flux Boundary
                KST[PCell]+=nuk/dn*Amag*Dot(Grids->Variables.KLSEGrad[PCell],PC);
                OmST[PCell]+=nuo/dn*Amag*Dot(Grids->Variables.OmLSEGrad[PCell],PC);

                //  Set Convective Flux Boundary
                KSMC->Elem[PCell][0]+=Flux;
                OmSMC->Elem[PCell][0]+=Flux;

                KST[PCell]+=-Flux*Dot(Grids->Variables.KLSEGrad[PCell],PC);
                OmST[PCell]+=-Flux*Dot(Grids->Variables.OmLSEGrad[PCell],PC);
            }break;
            case 4:// Symmetry;
            {
            }break;
            case 5:// Dirichlet Flow Function;
            {
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].CellCenter));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
                nuk=nu+(Grids->Variables.NuTurb[PCell])*sigstar;
                nuo=nu+(Grids->Variables.NuTurb[PCell])*sig;
                W=RigidBodyVel(Grids,FaceCenter);
                V=FlowFunction(FaceCenter,time);
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);
                vmag=Mag(V);

                //  Set Diffusive Flux Boundary
                KSMC->Elem[PCell][0]+=nuk/dn*Amag;
                OmSMC->Elem[PCell][0]+=nuo/dn*Amag;

                Kinlet=1.5*pow(Ti*vmag,2.0);
                l=1.0;
                Omegainlet=vmag/l;  //l: characteristic length or computational length

                KST[PCell]+=nuk/dn*Amag*(Kinlet-Dot(Grids->Variables.KLSEGrad[PCell],rC));
                OmST[PCell]+=nuo/dn*Amag*(Omegainlet-Dot(Grids->Variables.OmLSEGrad[PCell],rC));

                //  Set Convective Flux Boundary
                KST[PCell]+=-Flux*Kinlet;
//                EpsST[PCell]+=-Flux*0.09*pow(1.5*pow(Ti*vmag,2.0),2.0)/(Grids->Variables.NuTurb[PCell]);
                OmST[PCell]+=-Flux*Omegainlet;

            }break;
            case 6:// Static Wall;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
                nuk=nu+(Grids->Variables.NuTurb[PCell])*sigstar;
                nuo=nu+(Grids->Variables.NuTurb[PCell])*sig;
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);

                //  Set Diffusive Flux Boundary
                KSMC->Elem[PCell][0]+=nuk/dn*Amag;
                OmSMC->Elem[PCell][0]+=nuo/dn*Amag;

                Omegainlet=60.*nu/(0.75*pow(dn,2.0));
                KST[PCell]+=nuk/dn*Amag*(0.0-Dot(Grids->Variables.KLSEGrad[PCell],rC));
                OmST[PCell]+=nuo/dn*Amag*(Omegainlet-Dot(Grids->Variables.OmLSEGrad[PCell],rC));

                //  Set Convective Flux Boundary
                OmST[PCell]+=-Flux*Omegainlet;
            }break;
            case 7:// Free Slip Wall;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
                nuk=nu+(Grids->Variables.NuTurb[PCell])*sigstar;
                nuo=nu+(Grids->Variables.NuTurb[PCell])*sig;
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);

                //  Set Diffusive Flux Boundary
                KSMC->Elem[PCell][0]+=nuk/dn*Amag;
                OmSMC->Elem[PCell][0]+=nuo/dn*Amag;

                Omegainlet=60.*nu/(0.75*pow(dn,2.0));
                KST[PCell]+=nuk/dn*Amag*(0.0-Dot(Grids->Variables.KLSEGrad[PCell],rC));
                OmST[PCell]+=nuo/dn*Amag*(Omegainlet-Dot(Grids->Variables.OmLSEGrad[PCell],rC));

                //  Set Convective Flux Boundary
                OmST[PCell]+=-Flux*Omegainlet;
            }break;
        }
    }
}

//void ComputeWallFunction(struct Grid *Grids)
//{
//    int PCell,k,Type,FaceNum;
//    struct MyVector n,Area,FaceCenter,V,W,rC,Pf,PC;
//    double Amag,Flux,Value,dn,nu,nuk,nue,vmag;
//    double yp,des;
//    double ux,uy,uz,vx,vy,vz,wx,wy,wz,pt,lm;
//
//    for(k=0;k<Grids->NB;k++)
//    {
//        Type = Grids->BoundaryCells[k].Type;
//        if (Type==1)
//        {
//            PCell = Grids->BoundaryCells[k].CellNum;
//            FaceNum = Grids->BoundaryCells[k].FaceNum;
//            FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//            dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
//            nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
//            yp=dn*pow(0.09,0.25)*pow(Grids->Variables.NK[PCell],0.5)/nu;
//            switch (WallFunction)
//            {
//                case 'A' :
//                    {
//                        if(yp<=11.5)
//                        {
//                            Grids->Variables.NuTurb[PCell]=0.0;
//                        }
//                        else
//                        {
//                            des=DensCal(Grids->Variables.PAlpha[PCell]);
//                            Grids->Variables.NuTurb[PCell]=des*0.419*dn*pow(0.09,0.25)*pow(Grids->Variables.NK[PCell],0.5)/log(9.79*yp);
//                        }
//                    }break;
//                case 'B' :
//                    {
//                        ux=Grids->Variables.ULSEGrad[PCell].X;
//                        uy=Grids->Variables.ULSEGrad[PCell].Y;
//                        uz=Grids->Variables.ULSEGrad[PCell].Z;
//                        vx=Grids->Variables.VLSEGrad[PCell].X;
//                        vy=Grids->Variables.VLSEGrad[PCell].Y;
//                        vz=Grids->Variables.VLSEGrad[PCell].Z;
//                        wx=Grids->Variables.WLSEGrad[PCell].X;
//                        wy=Grids->Variables.WLSEGrad[PCell].Y;
//                        wz=Grids->Variables.WLSEGrad[PCell].Z;
//                        pt=(2.0*pow(ux,2.0)+2.0*pow(vy,2.0)+2.0*pow(wz,2.0)+pow(uy+vx,2.0)+pow(uz+wx,2.0)+pow(vz+wy,2.0));
//                        lm=0.419*dn*(1.0-exp(-yp/26.0));
//                        Grids->Variables.NuTurb[PCell]=pow(lm,2.0)*pow(pt,0.5);
//                    }break;
//            }
//        }
//    }
//}

void SolveKOEquation(struct Grid *Grids
                    ,char ConvectionMethod)
{
    struct MatrixCoefficient OmSMC,KSMC;
    double *KST,*OmST;
    int i,j,dummyCNT,dInt,dInte,dIntk,PCell,c,k,Type,FaceNum;
    double dn,yp,nu,des;
    struct MyVector FaceCenter;
    double ux,uy,uz,vx,vy,vz,wx,wy,wz,pt,lm;

    QMatrix LaspackMC;
    Vector ST;
    Vector Fee;

    MallocKOMatrixes(Grids,&LaspackMC,&OmSMC,&KSMC,&OmST,&KST,&ST,&Fee);

    i=0;
    do
    {
        UpdateKOmLSEGrad(Grids);
        dInt=dInte=dIntk=0;
        i++;

        MakeKMatrix(&KSMC,KST,Grids,ConvectionMethod,
                    TimeDiscretizationMethodForConvection,
                    TimeDiscretizationMethodForDiffusion);

        MakeOmMatrix(&OmSMC,OmST,Grids,ConvectionMethod,
                    TimeDiscretizationMethodForConvection,
                    TimeDiscretizationMethodForDiffusion);

        SetKOSourceTerms(KST,OmST,Grids);
        SetKOBC(&KSMC,&OmSMC,KST,OmST,Grids);

//        for(c=0;c<Grids->NC;c++)
//        {
//            printf("\n");
//            printf("KSMC[ %d , %d ] =     %e\n",c,c,KSMC.Elem[c][0]);
//
//            for(i=0;i<Grids->Cells[c].NoOfFace;i++)
//                if(Grids->Cells[c].Neighbor[i]>=0)
//                    printf("KSMC[ %d , %d ] =     %e\n",c,Grids->Cells[c].Neighbor[i],KSMC.Elem[c][i+1]);
//
//            printf("KST[ %d ] =     %e\n",c,KST[c]);
//        }
//
//        printf("#############################################");

//        for(c=0;c<Grids->NC;c++)
//        {
//            printf("\n");
//            printf("ESMC[ %d , %d ] =     %e\n",c,c,EpsSMC.Elem[c][0]);
//            for(i=0;i<Grids->Cells[c].NoOfFace;i++)
//                if(Grids->Cells[c].Neighbor[i]>=0)
//                    printf("ESMC[ %d , %d ] =     %e\n",c,Grids->Cells[c].Neighbor[i],EpsSMC.Elem[c][i]);
//            printf("EST[ %d ] =     %e\n",c,EpsST[c]);
//        }
//        getchar();

        SetRTCAccuracy(1e-8);
        ResetLaspackMatrix(&LaspackMC,Grids);
        CopyToLaspackMatrix(&LaspackMC,&KSMC,Grids);
        CopyToLaspackVector(&ST ,KST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NK,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
        dIntk=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NK,&Fee,Grids);

        ResetLaspackMatrix(&LaspackMC,Grids);
        CopyToLaspackMatrix(&LaspackMC,&OmSMC,Grids);
        CopyToLaspackVector(&ST ,OmST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NOm,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
        dInte=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NOm,&Fee,Grids);

        printf("##################  K Solver Iterations %d\n",dIntk);
        printf("##################  E Solver Iterations %d\n",dInte);
    }while((dInt>0)&&(i<4));

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.NuTurb[PCell]=Grids->Variables.NK[PCell]/Grids->Variables.NOm[PCell];
    }

    //ComputeWallFunction(Grids);

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.KLSEGrad[PCell]=Gradient_Disc_AxNode(Grids->Variables.NK,Grids->Variables.KLSEGrad,PCell,Grids);
    }

    FreeKOMatrixes(Grids,&LaspackMC,&OmSMC,&KSMC,&OmST,&KST,&ST,&Fee);
}

