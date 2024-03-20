double Upwind(double *Alpha,struct Grid *Grids,struct DonnerAcceptor DA,int FaceNum)
{
    struct MyVector D2A;
    double result;
    D2A=Sum(         Grids->Cells[DA.Face[FaceNum].a].CellCenter
           ,MinusVec(Grids->Cells[DA.Face[FaceNum].d].CellCenter));
    result = Alpha[DA.Face[FaceNum].a]
            -2.0*Dot(Gradient(Alpha,DA.Face[FaceNum].d,Grids),D2A);
    result=max(min(result,1.0),0.0);
    return result;
}

//calculation of normalized face donner alpha
double NormDonnerAlpha(double AU,struct DonnerAcceptor DA,double *Alpha,int FaceNum)
{
    double NDAlpha;
    NDAlpha=(Alpha[DA.Face[FaceNum].d]-AU)/(Alpha[DA.Face[FaceNum].a]-AU);
    if (isequal(Alpha[DA.Face[FaceNum].a],AU))
    {
        NDAlpha=100.;
        if (isequal(Alpha[DA.Face[FaceNum].d],AU))
            NDAlpha=1.;
    }
    return NDAlpha;
}

double CourantNum(struct Grid *Grids,int PCell)
{
    double CoNum=0.0,Flux;
    int FaceNum;
    struct MyVector FaceCenter,W,Area;
    int NoOfFace;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (FaceNum=0;FaceNum<NoOfFace;FaceNum++)
    {
        Area=Grids->Cells[PCell].Area[FaceNum];
        FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
        W=RigidBodyVel(Grids,FaceCenter);
        Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(W,Area);
        CoNum+=max(0.0,(-Flux*dt/Grids->Cells[PCell].Volume));
    }
    return CoNum;
}

//Face alpha_CBC (convection boundedness criteria) calculation
double CBCfAlphaf(double NDAlpha,double CoNum)
{
    double AlphaCBCf;
    if ((NDAlpha>=0.0)&&(NDAlpha<=1.0))
        AlphaCBCf=min(1.0,(NDAlpha/CoNum));
    else
        AlphaCBCf=NDAlpha;
    return AlphaCBCf;
}

//Face alpha_UQ (ULTIMATE-QUICKEST  ) calculation
double UQAlphaf(double AlphaCBCf,double NDAlpha,double CoNum)
{
    double AlphaUQf;
    if ((NDAlpha>=0.)&&(NDAlpha<=1.))
        AlphaUQf=min((((8.*(CoNum)*(NDAlpha))+((1.- CoNum)*(6.*NDAlpha+3.)))/8.),AlphaCBCf);
    else
        AlphaUQf=NDAlpha;
    return AlphaUQf;
}


double FThetha(struct MyVector Grad,struct MyVector d)
{
    double tethaf;
    double GradVal,dVal;
    GradVal=sqrt(Dot(Grad,Grad));
    dVal=sqrt(Dot(d,d));
    tethaf=acos(Dot(Grad,d)/(GradVal*dVal));
    if (isequal(GradVal*dVal,0.0)||(fabs(Dot(Grad,d)/(GradVal*dVal))>=1.0))
        tethaf=0.0;
    return tethaf;
}

//Calculation of the weight factor
double FGamma(double tethaf,double k)
{
 return min((1.0*(cos(2.0*tethaf)+1.0)/2.0),1.0);
}

double NormfaceAlphaCICSAM(struct Grid *Grids,double AlphaCBCf,double AlphaUQf,struct DonnerAcceptor DA,int FaceNum)
{
    double TethaF,GammaF,k;
    double NfAlpha;
    struct MyVector D2A;
    D2A = Sum(         Grids->Cells[DA.Face[FaceNum].a].CellCenter
             ,MinusVec(Grids->Cells[DA.Face[FaceNum].d].CellCenter));
    TethaF = FThetha(Gradient(Grids->Variables.Alpha,DA.Face[FaceNum].d,Grids),D2A);
    GammaF = FGamma(TethaF,k);
    NfAlpha=(GammaF*AlphaCBCf)+((1.-GammaF)*AlphaUQf);
    return NfAlpha;
}

double AlphaJ(double NDAlpha)
{
    double Alphaj;
    if ((NDAlpha<0.0) || (NDAlpha>1.0)) Alphaj=NDAlpha;
    else if ((0.0<=NDAlpha) && (NDAlpha<=0.5)) Alphaj=2.0*NDAlpha;
    else  Alphaj=1.0;
    return Alphaj;
}

double NormfaceAlphaHRIC(struct Grid *Grids,struct DonnerAcceptor DA,double NDAlpha,double CoNum,int FaceNum)
{
    double HRICfAlpha;
    double Alphaj,SAlphaj,TethaF;
    struct MyVector D2A;
    D2A = Sum(         Grids->Cells[DA.Face[FaceNum].a].CellCenter
             ,MinusVec(Grids->Cells[DA.Face[FaceNum].d].CellCenter));
    Alphaj=AlphaJ(NDAlpha);
    TethaF = FThetha(Gradient(Grids->Variables.Alpha,DA.Face[FaceNum].d,Grids),D2A);
    SAlphaj=(Alphaj*sqrt(fabs(cos(TethaF))))+(NDAlpha*(1.-sqrt(fabs(cos(TethaF)))));
    if (CoNum<0.3) HRICfAlpha=SAlphaj;
    else if (CoNum>0.7) HRICfAlpha=NDAlpha;
    else HRICfAlpha=NDAlpha+((SAlphaj-NDAlpha)*((0.7-CoNum)/0.4));
    return HRICfAlpha;
}

double NormfaceAlphaSTACS(struct Grid *Grids,struct DonnerAcceptor DA,double NDAlpha,double CoNum,int FaceNum)
{
    double STACSfAlpha;
    double SBEEfAlpha,STOICfAlpha,TethaF,cos4;
    struct MyVector D2A;
    D2A = Sum(         Grids->Cells[DA.Face[FaceNum].a].CellCenter
             ,MinusVec(Grids->Cells[DA.Face[FaceNum].d].CellCenter));
    SBEEfAlpha = STOICfAlpha = NDAlpha;
    if      ((NDAlpha>0.0)&&(NDAlpha<1.0)) SBEEfAlpha = 1.0;
    if      ((NDAlpha>0.0)&&(NDAlpha<0.5)) STOICfAlpha = 0.5+0.5*NDAlpha;
    else if ((NDAlpha>0.5)&&(NDAlpha<5.0/6.0)) STOICfAlpha = 3.0/8.0+0.75*NDAlpha;
    else if ((NDAlpha>5.0/6.0)&&(NDAlpha<1.0)) STOICfAlpha = 1.0;
    TethaF = FThetha(Gradient(Grids->Variables.Alpha,DA.Face[FaceNum].d,Grids),D2A);
    cos4 = cos(TethaF)*cos(TethaF)*cos(TethaF)*cos(TethaF);
    STACSfAlpha = SBEEfAlpha*cos4 + STOICfAlpha*(1.0-cos4);
    return STACSfAlpha;
}

double GetFBBDNormalFaceValue(double normDonnerPhi,double courantNum)
{
    double Result;

    if ((normDonnerPhi>0.0)&&(normDonnerPhi<=min(courantNum,1.0/3.0)))
        Result = (max((1.0/courantNum),3.0))*normDonnerPhi;
    else if ((normDonnerPhi>min(courantNum,1.0/3.0))&&(normDonnerPhi<1.0))
        Result = 1.0;
    else
        Result = normDonnerPhi;

    return Result;
}

double GetFBHRNormalFaceValue(double normDonnerPhi,double courantNum)
{
    double Result;

    if ((normDonnerPhi>0.0)&&(normDonnerPhi<=min((courantNum/(4.0*(1.0-courantNum))),1.0/8.0)))
        Result = (max((1.0/courantNum),3.0))*normDonnerPhi;
    else if ((normDonnerPhi>min((courantNum/(4.0*(1.0-courantNum))),1.0/8.0))&&(normDonnerPhi<=3.0/4.0))
        Result = normDonnerPhi + 1.0/4.0;
    else if ((normDonnerPhi>3.0/4.0)&&(normDonnerPhi<1.0))
        Result = 1.0;
    else
        Result = normDonnerPhi;
    return Result;
}

double NormfaceAlphaFBICS(struct Grid *Grids,struct DonnerAcceptor DA,double NDAlpha,double CoNum,int FaceNum)
{
    double FBICSfAlpha;
    double SBEEfAlpha,STOICfAlpha,TethaF,cos4;
    double PhiFBBD,PhiFBHR,GammaF;
    struct MyVector D2A;

    D2A = Sum(         Grids->Cells[DA.Face[FaceNum].a].CellCenter
             ,MinusVec(Grids->Cells[DA.Face[FaceNum].d].CellCenter));

    TethaF = FThetha(Gradient(Grids->Variables.Alpha,DA.Face[FaceNum].d,Grids),D2A);

    PhiFBBD=GetFBBDNormalFaceValue(NDAlpha,CoNum);
    PhiFBHR=GetFBHRNormalFaceValue(NDAlpha,CoNum);
    GammaF=pow(cos(TethaF),4.0);
    FBICSfAlpha=(GammaF*PhiFBBD)+((1.-GammaF)*PhiFBHR);
    return FBICSfAlpha;
}

//Calculation of the Betha value
double FaceBetha(double Nfalpha,double NDAlpha)
{
    double result;
    result=(Nfalpha-NDAlpha)/(1.0-NDAlpha);
 //occures when donner and acceptor alpha are equal.
    if (isequal(1.0,NDAlpha))
        result=1.0;
 //if there was no flux it must equal into zero for
    return result;
}

void MakeFaceAlpha(struct Grid *Grids,char VoFMethod)
{
    struct DonnerAcceptor DA;
    struct MyVector W,Area,FaceCenter;
    int PCell,s;
    double a[8],Betha[6],Vol,Flux;
    double AlphaU,NDAlpha,AlphaCBCf,AlphaUQf,NfAlpha;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        DonnAccCell(&DA,Grids,PCell);
        SetSVect(Betha,6,0.0);
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(s=0;s<NoOfFace;s++)
        {
            Area=Grids->Cells[PCell].Area[s];
            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
            W=RigidBodyVel(Grids,FaceCenter);
            Flux=Grids->Variables.NSF[PCell].Face[s]-Dot(W,Area);
            if (Grids->Cells[PCell].Neighbor[s]>=0)
            {
                //Upwind cell alpha for of current face
                AlphaU = min(max(Upwind(Grids->Variables.NAlpha,Grids,DA,s),0.0),1.0);
                //Donner normalized alpha value of current face
                NDAlpha = NormDonnerAlpha(AlphaU,DA,Grids->Variables.NAlpha,s);
                double CoNum;
                CoNum = CourantNum(Grids,PCell);
                //Face convection boundedness criteria alpha
                AlphaCBCf = CBCfAlphaf(NDAlpha,CoNum);
                //Face ULTIMATE QUICK alpha
                AlphaUQf = UQAlphaf(AlphaCBCf,NDAlpha,CoNum);
                //calculation of normalized face alpha used in face betha evaluation
                switch (VoFMethod)
                {
                    case 'C':
                        //CICSAM face value
                        NfAlpha = NormfaceAlphaCICSAM(Grids,AlphaCBCf,AlphaUQf,DA,s);
                    break;

                    case 'H':
                        //HRIC face value
                        NfAlpha = NormfaceAlphaHRIC(Grids,DA,NDAlpha,CoNum,s);
                    break;

                    case 'S':
                        //STACS face value
                        NfAlpha = NormfaceAlphaSTACS(Grids,DA,NDAlpha,CoNum,s);
                    break;

                    case 'F':
                        //FBICS face value
                        NfAlpha = NormfaceAlphaFBICS(Grids,DA,NDAlpha,CoNum,s);
                    break;
                }

                //Interpolation parameter for face alpha value
                Betha[s]=FaceBetha(NfAlpha,NDAlpha);
                //PCell Equation coefficients according to current cell face
            }


            //face value calculation for the next itteration source term calculation
            Grids->Variables.NAF[PCell].Face[s] = ((Grids->Cells[PCell].Neighbor[s]>=0)
                                                  ?( Grids->Variables.NAlpha[DA.Face[s].d]*(1.0-Betha[s])
                                                    +Grids->Variables.NAlpha[DA.Face[s].a]*Betha[s])
                                                  :(Grids->Variables.NAlpha[PCell]));
        }
    }
}

void MakeVoFMatrixCoeff(struct MatrixCoefficient *AMC,double *AST,struct Grid *Grids
                       ,char VoFMethod
                       ,char TemporalTermDiscretizationMethodForVoF
                       ,char TimeDiscretizationMethodForVoF)
{
    struct DonnerAcceptor DA;
    struct MyVector W,LaggedW,Area,FaceCenter;
    int PCell,s;
    double a[8],Betha[6],Vol,Flux,LaggedFlux;
    double AlphaU,NDAlpha,AlphaCBCf,AlphaUQf,NfAlpha;
    double Omega0,Omega1,Omega2,Gamma;
    switch (TemporalTermDiscretizationMethodForVoF)
    {
        case 'A' : {Omega0 =     1.0 ; Omega1 =     1.0 ; Omega2 =     0.0;}break;
        case 'B' : {Omega0 = 3.0/2.0 ; Omega1 = 4.0/2.0 ; Omega2 = 1.0/2.0;}break;
    }

    switch (TimeDiscretizationMethodForVoF)
    {
        case 'I': Gamma=0.0;break;
        case 'C': Gamma=0.5;break;
        case 'E': Gamma=1.0;break;
    }

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        DonnAccCell(&DA,Grids,PCell);
        SetSVect(a,8,0.0);
        Vol=Grids->Cells[PCell].Volume;
        a[0]=Vol/dt*Omega0;
        a[NoOfFace+1]=(Grids->Variables.Alpha[PCell]*Omega1 - Grids->Variables.PAlpha[PCell]*Omega2)*Vol/dt;
        SetSVect(Betha,6,0.0);

        for(s=0;s<NoOfFace;s++)
        {
            Area=Grids->Cells[PCell].Area[s];
            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
            W=RigidBodyVel(Grids,FaceCenter);
            Flux=Grids->Variables.NSF[PCell].Face[s]-Dot(W,Area);
            LaggedW=LaggedRigidBodyVel(Grids,FaceCenter);
            LaggedFlux=Grids->Variables.SF[PCell].Face[s]-Dot(LaggedW,Area);
            if (Grids->Cells[PCell].Neighbor[s]>=0)
            {
                //Upwind cell alpha for of current face
                AlphaU = min(max(Upwind(Grids->Variables.Alpha,Grids,DA,s),0.0),1.0);
                //Donner normalized alpha value of current face
                NDAlpha = NormDonnerAlpha(AlphaU,DA,Grids->Variables.Alpha,s);
                double CoNum;
                CoNum=CourantNum(Grids,PCell);
                //Face convection boundedness criteria alpha
                AlphaCBCf= CBCfAlphaf(NDAlpha,CoNum);
                //Face ULTIMATE QUICK alpha
                AlphaUQf = UQAlphaf(AlphaCBCf,NDAlpha,CoNum);
                //calculation of normalized face alpha used in face betha evaluation
                switch (VoFMethod)
                {
                    case 'C':
                        //CICSAM face value
                        NfAlpha = NormfaceAlphaCICSAM(Grids,AlphaCBCf,AlphaUQf,DA,s);
                    break;

                    case 'H':
                        //HRIC face value
                        NfAlpha = NormfaceAlphaHRIC(Grids,DA,NDAlpha,CoNum,s);
                    break;

                    case 'S':
                        //STACS face value
                        NfAlpha = NormfaceAlphaSTACS(Grids,DA,NDAlpha,CoNum,s);
                    break;

                    case 'F':
                        //FBICS face value
                        NfAlpha = NormfaceAlphaFBICS(Grids,DA,NDAlpha,CoNum,s);
                    break;
                }

                //Interpolation parameter for face alpha value
                Betha[s]=FaceBetha(NfAlpha,NDAlpha);
                //PCell Equation coefficients according to current cell face

                if (PCell==DA.Face[s].d)
                {
                    a[0]+=(1.0-Gamma)*(1.0-Betha[s])*Flux;
                    a[s+1]=(1.0-Gamma)*Betha[s]*Flux;
                }
                else
                {
                    a[0]+=(1.0-Gamma)*Betha[s]*Flux;
                    a[s+1]=(1.0-Gamma)*(1.0-Betha[s])*Flux;
                }
            }
            else
                a[0]+=Flux;


            a[NoOfFace+1]+=-Gamma*((Grids->Cells[PCell].Neighbor[s]>=0)
                                 ?((Grids->Variables.Alpha[DA.Face[s].d]*(1.0-Betha[s])+Grids->Variables.Alpha[DA.Face[s].a]*Betha[s])*Flux)
                                 :(0.0));

        }

        AMC->Elem[PCell][0]=a[0];
        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                AMC->Elem[PCell][s+1]=a[s+1];
        AST[PCell]=a[NoOfFace+1];

    }
}

void MallocVoFMatrixes(struct Grid *Grids
                      ,QMatrix *LaspackAMC
                      ,struct MatrixCoefficient *AMC
                      ,double **MyAST
                      ,Vector *AST
                      ,Vector *A)
{
    int NC=0,NCol=0,PCell;
    NC=Grids->NC;

    (*AMC).Elem=(double**) malloc(sizeof(double*)*NC);
    for(PCell=0;PCell<NC;PCell++)
    {

        NCol=Grids->Cells[PCell].NoOfFace+1;
        (*AMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
    }
    (*MyAST) = (double*) malloc(sizeof(double)*NC);

    Q_Constr(LaspackAMC,"GAMC",NC,False,Rowws,Normal,True);
    V_Constr(AST,"AST",NC,Normal,True);
    V_Constr(A,"A",NC,Normal,True);
}

void FreeVoFMatrixes(struct Grid *Grids
                    ,QMatrix *LaspackAMC
                    ,struct MatrixCoefficient *AMC
                    ,double **MyAST
                    ,Vector *AST
                    ,Vector *A)
{
    int NC=0,PCell;
    NC=Grids->NC;
    for(PCell=0;PCell<NC;PCell++)
        free((*AMC).Elem[PCell]);
    free((*AMC).Elem);

    free(*MyAST);

    Q_Destr(LaspackAMC);
    V_Destr(AST);
    V_Destr(A);
}

void FilterAlpha(struct Grid *Grids)
{
    int i;
    for(i=0;i<Grids->NC;i++)
        if (Grids->Variables.NAlpha[i]<0.0) Grids->Variables.NAlpha[i]=0.0;
        else if (Grids->Variables.NAlpha[i]>1.0) Grids->Variables.NAlpha[i]=1.0;
}

int IsClose2InterFace(struct Grid *Grids,int PCell)
{
    double beta=0.00075;
    struct MyVector NewGrad;
    NewGrad=Gradient(Grids->Variables.NAlpha,PCell,Grids);
    return ((sqrt(Dot(NewGrad,NewGrad))*pow(Grids->Cells[PCell].Volume,1./3.))>beta);
}

void FilterNonClose2InterFaceAlpha(struct Grid *Grids)
{
    int PCell;
    for(PCell=0;PCell<Grids->NC;PCell++)
        if (!IsClose2InterFace(Grids,PCell))
            Grids->Variables.NAlpha[PCell]=((Grids->Variables.NAlpha[PCell]<0.5)?(0.0):(1.0));
}

void SetVOFBC(struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector n,Area,FaceCenter,V,W,rC,Pf,PC;
    double Amag,Flux,Value,dn,nu;
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
            }break;
            case 3:// Out-Flow;
            {
            }break;
            case 4:// Symmetry;
            {
            }break;
            case 5:// Dirichlet Flow Function;
            {
                Grids->Variables.NAlpha[PCell]=1.0;
            }break;
            case 6:// Static Wall;
            {
            }break;
        }
    }
}
//Free Surface Equation Solving results in Volume Fraction(Alpha) distribution
//void FreeSurface(double *Alpha,struct SurfFlux *Alphaf,struct SurfFlux *SF,double *ao,char VoFMethod)

void VoFSolver(struct Grid *Grids
              ,char VoFMethod
              ,char TemporalTermDiscretizationMethodForVoF
              ,char TimeDiscretizationMethodForVoF)
{
    struct MatrixCoefficient AMC;
    double *MyAST;
    int i;
    QMatrix LaspackAMC;
    Vector AST;
    Vector A;

//    printf("VoF started ...\n");
    //Memory allocation
    MallocVoFMatrixes(Grids, &LaspackAMC, &AMC,&MyAST, &AST, &A);

    MakeVoFMatrixCoeff(&AMC,MyAST,Grids,VoFMethod
                                 ,TemporalTermDiscretizationMethodForVoF
                                 ,TimeDiscretizationMethodForVoF);
    CopyToLaspackMatrix(&LaspackAMC,&AMC,Grids);
    CopyToLaspackVector(&A,Grids->Variables.NAlpha,Grids);
    CopyToLaspackVector(&AST,MyAST,Grids);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SetRTCAccuracy(1e-8);
    A=*BiCGSTABIter(&LaspackAMC,&A,&AST,100,ILUPrecond,1.0);
    int dInt;
    dInt=GetLastNoIter();
    printf("##################  VoF Solver Iterations %d\n",dInt);
    CopyFromLaspackVector(Grids->Variables.NAlpha,&A,Grids);
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//    SetVOFBC(Grids);

    FilterAlpha(Grids);
//    FilterNonClose2InterFaceAlpha(Grids);
    MakeFaceAlpha(Grids,VoFMethod);
    FreeVoFMatrixes(Grids, &LaspackAMC, &AMC, &MyAST, &AST, &A);
//    printf("VoF ended ...\n");
}
