double FindAbsMax(double *A,int NC)
{
    int i;
    double maxVal;
    maxVal=fabs(A[0]);
    for (i=1;i<NC;i++)
        if (fabs(A[i])>maxVal)
            maxVal=fabs(A[i]);
    return maxVal;
}


struct MyVector Sum(struct MyVector A ,struct MyVector B)
{
    A.X+=B.X;
    A.Y+=B.Y;
    A.Z+=B.Z;
    return(A);
}

struct MyVector Minus(struct MyVector A,struct MyVector B)
{
    A.X-=B.X;
    A.Y-=B.Y;
    A.Z-=B.Z;
    return(A);
}

struct MyVector MinusVec(struct MyVector A)
{
  A.X=-A.X;
  A.Y=-A.Y;
  A.Z=-A.Z;
  return(A);
}

struct MyVector Cross(struct MyVector A ,struct MyVector B)
{
    struct MyVector C;
    C.X=A.Y*B.Z-A.Z*B.Y;
    C.Y=A.Z*B.X-A.X*B.Z;
    C.Z=A.X*B.Y-A.Y*B.X;
    return(C);
}

struct MyVector SVP(double a,struct MyVector SUM)
{
    SUM.X*=a;
    SUM.Y*=a;
    SUM.Z*=a;
    return(SUM);
}

inline double Mag(struct MyVector Vec)
{
    return sqrt(Dot(Vec,Vec));
}

inline struct MyVector Normalized(struct MyVector Vec)
{
    return SVP(1.0/Mag(Vec),Vec);
}

double SelectDens(struct MyVector V,struct MyVector Added,double PDens)
{
    double result,MagV;
    result = PDens;
//    if (time>dt)
    MagV=Mag(V);
    if (MagV>0.1)
        result=((Mag(Added)/MagV)>2.0)?DensW:PDens;
    return result;
}

void EqualSurfFlux(struct SurfFlux *SF,struct SurfFlux *NSF,struct Grid *Grids)
{
    int i,PCell,NC,NoOfFace;
    NC=Grids->NC;
    for(PCell=0;PCell<NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(i=0;i<NoOfFace;i++)
            SF[PCell].Face[i]=NSF[PCell].Face[i];
    }
}

void GetToken(char *string, FILE *fpMesh)
{
    char dummychar;
    sprintf(string,"");
    dummychar=getc(fpMesh);
    if (feof(fpMesh)==0)
    {
        while ((dummychar==' ')||(dummychar=='\n'))
          dummychar=getc(fpMesh);
        do
        {
          sprintf(string,"%s%c",string,dummychar);
          dummychar=getc(fpMesh);
        }while ((dummychar!=' ')&&(dummychar!='\n'));
    }
}

/*
void TakeAWord(char *result,FILE *filePointer)
{
    char ch;
    int IsOver=0;
    ch = getc( filePointer );
    sprintf(result,"");
    while (!(ch==EOF))
    {
        if ((ch==32)&&(!IsOver))  IsOver=0;
        else if (!(ch==32))
        {
            sprintf(result,"%s%c",result,ch);
            IsOver=1;
        }
        else break;
        ch = getc( filePointer );
    }
}
*/

void EqualSVect(double *Vect,double *NVect,int r)
{
  int i;
  for(i=0;i<r;i++)
    Vect[i]=NVect[i];
}

void SumSVect(double *Vect1,double *Vect2,double w2,double *Vect3,double w3,int r)
{
  int i;
  for(i=0;i<r;i++)
    Vect1[i]=w2*Vect2[i]+w3*Vect3[i];
}

void SetSVect(double *Vect,int r,double a)
{
  int i;
  for(i=0;i<r;i++)
    Vect[i]=a;
}

void SetIntSVect(int *Vect,int r,int val)
{
  int i;
  for(i=0;i<r;i++)
    Vect[i]=val;
}

struct MyVector RigidBodyVel(struct Grid *Grids,struct MyVector Position)
{
    struct MyVector V,R;
    R=Sum(MinusVec(Grids->DynaVariables.NewMassCenter),Position);
    V=Sum(Grids->DynaVariables.NewLinVel,Cross(Grids->DynaVariables.NewAngVel,R));
    return V;
}

struct MyVector LaggedRigidBodyVel(struct Grid *Grids,struct MyVector Position)
{
    struct MyVector V,R;
    R=Sum(MinusVec(Grids->DynaVariables.MassCenter),Position);
    V=Sum(Grids->DynaVariables.LinVel,Cross(Grids->DynaVariables.AngVel,R));
    return V;
}

int FindNearestCellTo(struct MyVector Target,struct Grid *Grids)
{
    int result=0,PCell;
    struct MyVector R;
    double RMag,IRMag;
    R = Sum(MinusVec(Target),Grids->Cells[result].CellCenter);
    IRMag= sqrt(Dot(R,R));
    for (PCell=1;PCell<Grids->NC;PCell++)
    {
        R = Sum(MinusVec(Target),Grids->Cells[PCell].CellCenter);
        RMag= sqrt(Dot(R,R));
        if (RMag<=IRMag)
        {
          result=PCell;
          IRMag=RMag;
        }
    }
    return result;
}

double Weight(struct Grid *Grids,int PCell,int FaceNum)
{
    double result=1.0;
    struct MyVector C2C,C2F;
    if (Grids->Cells[PCell].Neighbor[FaceNum]>=0)
    {
      C2C=Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[Grids->Cells[PCell].Neighbor[FaceNum]].CellCenter));
      C2F=Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[PCell].FaceCoincide[FaceNum]));
      result=1.0-sqrt(Dot(C2F,C2F)/Dot(C2C,C2C));
    }
    return result;
}

double GetVolumeFraction(double Phi,double Vol)
{
    double H;
    double epsilon;
    epsilon=0.5*pow(Vol,1./3.);
    if (Phi<-epsilon)
        H=0.0;
    else if (Phi>epsilon)
        H=1.0;
    else
        H=0.5+Phi/(2.0*epsilon)+1.0/(2.0*pi)*sin(pi*Phi/epsilon);
    return (1.0-H);
}


struct MyVector GetCoincidance(struct MyVector Node1,
                               struct MyVector Node2,
                               struct MyVector FaceCenter,
                               struct MyVector Normal)
{
// Line Equation :
//   (X-X1)   (Y-Y1)   (Z-Z1)
//   ------ = ------ = ------ = Landa
//     LX       LY       LZ

// Surface Equation :
// Nx * (X-X0) + Ny * (Y-Y0) + Nz * (Z-Z0) = 0

// Nx * (LX * Landa + X1 - X0) +
// Ny * (LY * Landa + Y1 - Y0) +
// Nz * (LZ * Landa + Z1 - Z0) = 0

//            Nx*(X1-X0) + Ny*(Y1-Y0) + Nz*(Z1-Z0)
// Landa = - --------------------------------------
//                  Nx*LX + Ny*LY + Nz*LZ
    struct MyVector Result;
    double LX,LY,LZ,Landa;
    LX=Node1.X-Node2.X;
    LY=Node1.Y-Node2.Y;
    LZ=Node1.Z-Node2.Z;

    Landa = - (Normal.X*(Node1.X-FaceCenter.X)+
               Normal.Y*(Node1.Y-FaceCenter.Y)+
               Normal.Z*(Node1.Z-FaceCenter.Z))/
              (Normal.X*LX+Normal.Y*LY+Normal.Z*LZ);
    Result.X = Landa * LX + Node1.X;
    Result.Y = Landa * LY + Node1.Y;
    Result.Z = Landa * LZ + Node1.Z;
    return Result;
}

struct MyVector GetProjectedCenter(struct MyVector Node,
                                   struct MyVector LineNode,
                                   struct MyVector LineVector)
{
    double s;
    struct MyVector Result,Diff;
    Diff = Minus(Node,LineNode);
    s=Dot(Diff,LineVector);
    Result = Sum (LineNode,SVP(s,LineVector));
    return Result;
}

void GetFaceValueParts(double *aP,double *aN,double *ST
                        ,struct MyVector *LSEGrad,int PCell,int FaceNum,struct Grid *Grids)
/* Fee_f = aP * Fee_P + aP * Fee_N + ST */
{
    struct MyVector PPrime,NPrime,P2F,N2F,FaceCoincide,FaceCenter,CVector;
    double delta,deltaP,G;
    double FeeP,FeeNgb,Feef;
    int NgbCell;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*ST=0.0;
    FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
    FaceCoincide = Grids->Cells[PCell].FaceCoincide[FaceNum];
    CVector = Minus(FaceCenter,FaceCoincide);
    PPrime = Sum(CVector,Grids->Cells[PCell].CellCenter);
    P2F = Minus(FaceCenter,PPrime);
    deltaP = Mag(P2F);
    if (Grids->Cells[PCell].Neighbor[FaceNum]>=0)
    {
        NPrime = Sum(CVector,Grids->Cells[NgbCell].CellCenter);
        N2F = Minus(FaceCenter,NPrime);
        delta = Mag(N2F);
        G = (delta)/(delta+deltaP);
        *aP = G;
        *aN = 1.0-G;
        *ST =      G  * Dot(LSEGrad[  PCell],CVector)
             +(1.0-G) * Dot(LSEGrad[NgbCell],CVector);
    }
    else
    {
        *aP = 1.0;
        *ST = Dot(LSEGrad[PCell],P2F);
    }
}

GetAllUpwindFaceValueParts(double *aP,double *aN,
                           double *UST,double *VST,double *WST,
                           double *U,struct MyVector *LSEGradU,
                           double *V,struct MyVector *LSEGradV,
                           double *W,struct MyVector *LSEGradW,
                           int PCell,int FaceNum,struct Grid *Grids)
{
    struct MyVector PPrime,NPrime,P2F,N2F,FaceCoincide,FaceCenter,CVector;
    double Flux,G;
    double FeeP,FeeNgb,Feef;
    int NgbCell;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*UST=*VST=*WST=0.0;


    FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
    FaceCoincide = Grids->Cells[PCell].FaceCoincide[FaceNum];
    CVector = Minus(FaceCenter,FaceCoincide);
    PPrime = Sum(CVector,Grids->Cells[PCell].CellCenter);
    P2F = Minus(FaceCenter,PPrime);
    Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-
         Dot(RigidBodyVel(Grids,FaceCenter),Grids->Cells[PCell].Area[FaceNum]);

    if (Grids->Cells[PCell].Neighbor[FaceNum]>=0)
    {
        G=(Flux>0.0)?1.0:0.0;
        *aP = G;
        *aN = 1.0-G;
        *UST=      G  * Dot(LSEGradU[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradU[NgbCell],CVector);
        *VST=      G  * Dot(LSEGradV[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradV[NgbCell],CVector);
        *WST=      G  * Dot(LSEGradW[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradW[NgbCell],CVector);
    }
    else
    {
        *aP = 1.0;
        *UST = Dot(LSEGradU[PCell],P2F);
        *VST = Dot(LSEGradV[PCell],P2F);
        *WST = Dot(LSEGradW[PCell],P2F);
    }
}

void GetAllFaceValueParts(double *aP,double *aN,
                          double *UST,double *VST,double *WST,
                          double *U,struct MyVector *LSEGradU,
                          double *V,struct MyVector *LSEGradV,
                          double *W,struct MyVector *LSEGradW,
                          int PCell,int FaceNum,struct Grid *Grids)
/* Fee_f = aP * Fee_P + aP * Fee_N + ST */
{
    struct MyVector PPrime,NPrime,P2F,N2F,FaceCoincide,FaceCenter,CVector;
    double delta,deltaP,G;
    double FeeP,FeeNgb,Feef;
    int NgbCell;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*UST=*VST=*WST=0.0;


    FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
    FaceCoincide = Grids->Cells[PCell].FaceCoincide[FaceNum];
    CVector = Minus(FaceCenter,FaceCoincide);
    PPrime = Sum(CVector,Grids->Cells[PCell].CellCenter);
    P2F = Minus(FaceCenter,PPrime);
    deltaP = Mag(P2F);
    if (Grids->Cells[PCell].Neighbor[FaceNum]>=0)
    {
        NPrime = Sum(CVector,Grids->Cells[NgbCell].CellCenter);
        N2F = Minus(FaceCenter,NPrime);
        delta = Mag(N2F);
        G = (delta)/(delta+deltaP);
        *aP = G;
        *aN = 1.0-G;
        *UST=      G  * Dot(LSEGradU[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradU[NgbCell],CVector);
        *VST=      G  * Dot(LSEGradV[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradV[NgbCell],CVector);
        *WST=      G  * Dot(LSEGradW[  PCell],CVector)
             +(1.0-G) * Dot(LSEGradW[NgbCell],CVector);
    }
    else
    {
        *aP = 1.0;
        *UST = Dot(LSEGradU[PCell],P2F);
        *VST = Dot(LSEGradV[PCell],P2F);
        *WST = Dot(LSEGradW[PCell],P2F);
    }
}

double GetFaceValueParts_Disc(double *aP,double *aN,double *ST
                             ,double *Fee,struct MyVector *LSEGrad,int PCell,int FaceNum,struct Grid *Grids)
/* Fee_f = aP * Fee_P + aP * Fee_N + ST */
{
    struct MyVector PPrime,NPrime,FaceCenter,CPVector,CNVector,P2F;
    double deltaP,deltaN
          ,gammaP,gammaN
          ,DensP,DensN;
    int NgbCell,i,MatchedFaceNum;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*ST=0.0;

    FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];

    PPrime = Grids->Cells[PCell].ProjectedCenter[FaceNum];
    CPVector = Minus(PPrime,Grids->Cells[PCell].CellCenter);
    deltaP = Mag(Minus(PPrime,FaceCenter));
    DensP = DensCal(Grids->Variables.Alpha[PCell]);
    gammaP = deltaP * DensP;

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
        CNVector = Minus(NPrime,Grids->Cells[NgbCell].CellCenter);
        deltaN = Mag(Minus(NPrime,FaceCenter));
        DensN = DensCal(Grids->Variables.Alpha[NgbCell]);
        gammaN = deltaN * DensN;


        *aP = gammaN/(gammaN+gammaP);
        *aN = gammaP/(gammaN+gammaP);
        *ST = (+Dot(LSEGrad[  PCell],CPVector)*gammaN
               +Dot(LSEGrad[NgbCell],CNVector)*gammaP)/(gammaN+gammaP);

    }
    else
    {
        P2F = Minus(FaceCenter,Grids->Cells[PCell].CellCenter);
        *aP = 1.0;
        *ST = Dot(LSEGrad[PCell],P2F);
    }
}

double FaceValue_AxNode(double *Fee,struct MyVector *LSEGrad,struct Grid *Grids,int PCell,int FaceNum)
{
    double aP,aN,ST,res;
    int NgbCell;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,FaceNum,Grids);
    if (NgbCell>=0)
        res = aP*Fee[PCell] + aN*Fee[NgbCell] + ST;
    else
        res = aP*Fee[PCell] + ST;
    return res;
}

double FaceValue_Boundary(double *Fee,struct MyVector *LSEGrad,struct Grid *Grids,int PCell,int FaceNum)
{
    double res;
    struct MyVector CC,FC,PF;
    CC=Grids->Cells[PCell].CellCenter;
    FC=Grids->Cells[PCell].FaceCenter[FaceNum];
    PF=Minus(FC,CC);
    res = Fee[PCell] + Dot(LSEGrad[PCell],PF);
    return res;
}

double FaceValue(double *Fee,struct Grid *Grids,int PCell,int FaceNum)
{
    struct MyVector PCoord,NgbCoord,FaceCenter;
    double delta,deltaP,G;
    double FeeP,FeeNgb,Feef;
    FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
    FeeP=Fee[PCell];
    PCoord=Grids->Cells[PCell].CellCenter;

    if (Grids->Cells[PCell].Neighbor[FaceNum]>=0)
    {
        NgbCoord =Grids->Cells[Grids->Cells[PCell].Neighbor[FaceNum]].CellCenter;
        FeeNgb = Fee[Grids->Cells[PCell].Neighbor[FaceNum]];
        deltaP = sqrt(Dot(Sum(PCoord,MinusVec(FaceCenter)),Sum(PCoord,MinusVec(FaceCenter))));
        delta = sqrt(Dot(Sum(NgbCoord,MinusVec(FaceCenter)),Sum(NgbCoord,MinusVec(FaceCenter))));
        G = (delta)/(delta+deltaP);
        Feef = (G*FeeP+(1.0-G)*FeeNgb);
    }
    else
        Feef=FeeP;
    return Feef;
}

struct MyVector SurfaceGradient(struct SurfFlux *AF,int PCell,struct Grid *Grids)
{
    int i,NoOfFace;
    struct MyVector Result,Area;
    double Vol,Hf;
    Result.X=Result.Y=Result.Z=0.0;
    Vol=Grids->Cells[PCell].Volume;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
    {
        Area=Grids->Cells[PCell].Area[i];
        Hf=AF[PCell].Face[i];
        Result.X+=Hf*Area.X/Vol;
        Result.Y+=Hf*Area.Y/Vol;
        Result.Z+=Hf*Area.Z/Vol;
    }
    if (Dimension=='2')
        Result.Y=0.0;
    return Result;
}

struct MyVector Gradient(double *H,int PCell,struct Grid *Grids)
{
    int i,NoOfFace;
    struct MyVector Result,FaceCenter,Area;
    double Vol,Hf;
    Result.X=Result.Y=Result.Z=0.0;
    Vol=Grids->Cells[PCell].Volume;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
    {
        FaceCenter=Grids->Cells[PCell].FaceCenter[i];
        Area=Grids->Cells[PCell].Area[i];
        Hf=FaceValue(H,Grids,PCell,i);
        Result.X+=Hf*Area.X/Vol;
        Result.Y+=Hf*Area.Y/Vol;
        Result.Z+=Hf*Area.Z/Vol;
    }
    if (Dimension=='2')
        Result.Y=0.0;
    return Result;
}

//struct MyVector SpecialGradient(double *H,int PCell,struct Grid *Grids)
//{
//    int i,NgbCell,NoOfFace;
//    struct MyVector Result,FaceCenter,Area,C2F,Grad;
//    double Vol,Hf,W,NewW,DensP,DensNgb;
//    Result.X=Result.Y=Result.Z=0.0;
//    Vol=Grids->Cells[PCell].Volume;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//
//    for (i=0;i<NoOfFace;i++)
//    {
//        NgbCell=((Grids->Cells[PCell].Neighbor[i]>=0)?(Grids->Cells[PCell].Neighbor[i]):PCell);
//        FaceCenter=Grids->Cells[PCell].FaceCenter[i];
//        Area=Grids->Cells[PCell].Area[i];
//
//        W=Weight(Grids,PCell,i);
//        DensP=DensCal(Grids->Variables.Alpha[PCell]);
//        DensNgb=DensCal(Grids->Variables.Alpha[NgbCell]);
//
//        NewW=(DensNgb*(1.0-W))/(DensNgb*(1.0-W)+DensP*W);
//        Hf=NewW*H[PCell]+(1.0-NewW)*H[NgbCell];
//
//        double AddedPressureX=0.0,AddedPressureY=0.0,AddedPressureZ=0.0;
//
//        if (Grids->Cells[PCell].Neighbor[i]<0)
//        {
////            Grad=Gradient(H,PCell,Grids);
////            C2F=Sum(MinusVec(Grids->Cells[PCell].CellCenter),FaceCenter);
////            Hf=H[PCell]+Dot(Grad,C2F);
//            Hf=H[PCell];
//            AddedPressureX=DensP*Gravity.X*(FaceCenter.X-Grids->Cells[PCell].CellCenter.X);
//            AddedPressureY=DensP*Gravity.Y*(FaceCenter.Y-Grids->Cells[PCell].CellCenter.Y);
//            AddedPressureZ=DensP*Gravity.Z*(FaceCenter.Z-Grids->Cells[PCell].CellCenter.Z);
//        }
//        Result.X+=(Hf+AddedPressureX)*Area.X/Vol;
//        Result.Y+=(Hf+AddedPressureY)*Area.Y/Vol;
//        Result.Z+=(Hf+AddedPressureZ)*Area.Z/Vol;
//    }
//
//    return Result;
//}

struct MyVector Gradient_Disc_AxNode(double *H,struct MyVector *LSEGrad,int PCell,struct Grid *Grids)
{
    int i,NgbCell,NoOfFace;
    struct MyVector Result,FaceCenter,Area,C2F,Grad;
    double Vol,Hf,W,NewW,DensP,DensNgb,aP,aN,ST;
    Result.X=Result.Y=Result.Z=0.0;
    Vol=Grids->Cells[PCell].Volume;
    NoOfFace=Grids->Cells[PCell].NoOfFace;

    for (i=0;i<NoOfFace;i++)
    {
        NgbCell=((Grids->Cells[PCell].Neighbor[i]>=0)?(Grids->Cells[PCell].Neighbor[i]):PCell);
        FaceCenter=Grids->Cells[PCell].FaceCenter[i];
        Area=Grids->Cells[PCell].Area[i];

        GetFaceValueParts_Disc(&aP,&aN,&ST,H,LSEGrad,PCell,i,Grids);
        Hf=aP*H[PCell]+aN*H[NgbCell]+ST;

        Result.X+=(Hf)*Area.X/Vol;
        Result.Y+=(Hf)*Area.Y/Vol;
        Result.Z+=(Hf)*Area.Z/Vol;
    }
    if (Dimension=='2')
        Result.Y=0.0;
    return Result;
}

struct MyVector Gradient_AxNode(double *H,struct MyVector *LSEGrad,int PCell,struct Grid *Grids)
{
    int i,NgbCell,NoOfFace;
    struct MyVector Result,FaceCenter,Area,C2F,Grad;
    double Vol,Hf,W,NewW,DensP,DensNgb,aP,aN,ST;
    Result.X=Result.Y=Result.Z=0.0;
    Vol=Grids->Cells[PCell].Volume;
    NoOfFace=Grids->Cells[PCell].NoOfFace;

    for (i=0;i<NoOfFace;i++)
    {
        NgbCell=((Grids->Cells[PCell].Neighbor[i]>=0)?(Grids->Cells[PCell].Neighbor[i]):PCell);
        FaceCenter=Grids->Cells[PCell].FaceCenter[i];
        Area=Grids->Cells[PCell].Area[i];

        GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,i,Grids);

        Hf=aP*H[PCell]+aN*H[NgbCell]+ST;

        Result.X+=(Hf)*Area.X/Vol;
        Result.Y+=(Hf)*Area.Y/Vol;
        Result.Z+=(Hf)*Area.Z/Vol;
    }
    if (Dimension=='2')
        Result.Y=0.0;
    return Result;
}


void GetAllNormalGradParts(double *aP,double *aN,
                           double *UST,double *VST,double *WST,
                           double *U,struct MyVector *LSEGradU,
                           double *V,struct MyVector *LSEGradV,
                           double *W,struct MyVector *LSEGradW,
                           int PCell,int FaceNum,struct Grid *Grids)

{
    struct MyVector PPrime,NPrime,FaceCenter,CPVector,CNVector,P2F;
    double delta;
    double FeeP,FeeNgb,Feef;
    int NgbCell,i,MatchedFaceNum;
//  dfee/dn = aP * Fee (PCell) + aN * Fee (NgbCell) + ST
    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*UST=*VST=*WST=0.0;

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
        CNVector = Minus(NPrime,Grids->Cells[NgbCell].CellCenter);
        delta = Mag(Minus(PPrime,NPrime));
        *aP = -1.0/delta;
        *aN = +1.0/delta;
        *UST = (Dot(LSEGradU[NgbCell],CNVector)
               -Dot(LSEGradU[  PCell],CPVector))/delta;
        *VST = (Dot(LSEGradV[NgbCell],CNVector)
               -Dot(LSEGradV[  PCell],CPVector))/delta;
        *WST = (Dot(LSEGradW[NgbCell],CNVector)
               -Dot(LSEGradW[  PCell],CPVector))/delta;

    }
//    else
//    {
//        dfee/dn=(Fee(PPrime)-Fee(FC))/deltaF
//        dfee/dn=(Fee(P)+g(P)*(rPPrime-rP)-Fee(P)-g(P)*(rF-rP))/deltaF
//        dfee/dn=(g(P)*(rPPrime-rF))/deltaF
//        P2F=Minus(PPrime,Grids->Cells[PCell].FaceCenter[FaceNum]);
//        *ST = 0.0;
//    }
}


void GetNormalGradParts(double *aP,double *aN,double *ST,
                        double *Fee,struct MyVector *LSEGrad,
                        int PCell,int FaceNum,struct Grid *Grids)
{
    struct MyVector PPrime,NPrime,FaceCenter,CPVector,CNVector,P2F;
    double delta;
    double FeeP,FeeNgb,Feef;
    int NgbCell,i,MatchedFaceNum;
//  dfee/dn = aP * Fee (PCell) + aN * Fee (NgbCell) + ST
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
        CNVector = Minus(NPrime,Grids->Cells[NgbCell].CellCenter);
        delta = Mag(Minus(PPrime,NPrime));
        *aP = -1.0/delta;
        *aN = +1.0/delta;
        *ST = (Dot(LSEGrad[NgbCell],CNVector)
              -Dot(LSEGrad[  PCell],CPVector))/delta;

    }
//    else
//    {
//        dfee/dn=(Fee(PPrime)-Fee(FC))/deltaF
//        dfee/dn=(Fee(P)+g(P)*(rPPrime-rP)-Fee(P)-g(P)*(rF-rP))/deltaF
//        dfee/dn=(g(P)*(rPPrime-rF))/deltaF
//        P2F=Minus(PPrime,Grids->Cells[PCell].FaceCenter[FaceNum]);
//        *ST = 0.0;
//    }
}

void GetNormalGradParts_Disc_AxNode(double *aP,double *aN,double *ST,
                                    double *Fee,struct MyVector *LSEGrad,
                                    int PCell,int FaceNum,struct Grid *Grids)
/* Fee_f = aP * Fee_P + aP * Fee_N + ST */
{
    struct MyVector PPrime,FaceCenter,CPVector,P2F;
    double deltaP,gammaP,DensP;
    int NgbCell,i,MatchedFaceNum;

    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    *aP=*aN=*ST=0.0;

    FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
    PPrime = Grids->Cells[PCell].ProjectedCenter[FaceNum];
    CPVector = Minus(PPrime,Grids->Cells[PCell].CellCenter);
    deltaP = Mag(Minus(PPrime,FaceCenter));
    DensP = DensCal(Grids->Variables.Alpha[PCell]);
    gammaP = deltaP * DensP;

    if (NgbCell>=0)
    {
        GetFaceValueParts_Disc(aP,aN,ST,Fee,LSEGrad,PCell,FaceNum,Grids);

        *aP = (*aP - 1.0)/gammaP;
        *aN = (*aN      )/gammaP;
        *ST = (*ST - Dot(LSEGrad[PCell],CPVector))/gammaP;
    }
//    else
//    {
//
//    }
}

void GetNormalGradParts_Disc_OverRelaxed(double *aP,double *aN,double *ST
                                        ,double *Fee,struct MyVector *LSEGrad,int PCell,int FaceNum,struct Grid *Grids)
/* Fee_f = aP * Fee_P + aP * Fee_N + ST */
{
    int NgbCell;
    double DensF,minDens,maxDens,NuF,W,Amag;
    struct MyVector S,K,Delta,D,GradP,GradNgb,GradF;
//    GradP = Gradient_Disc_AxNode(Fee,LSEGrad,PCell,Grids);
    GradP=LSEGrad[PCell];
    NgbCell=Grids->Cells[PCell].Neighbor[FaceNum];
    (*aP)=(*aN)=(*ST)=0.0;
    if (NgbCell>=0)
    {
        minDens=min(DensCal(Grids->Variables.Alpha[PCell]),DensCal(Grids->Variables.Alpha[NgbCell]));
        maxDens=max(DensCal(Grids->Variables.Alpha[PCell]),DensCal(Grids->Variables.Alpha[NgbCell]));
//        GradNgb = Gradient_Disc_AxNode(Fee,LSEGrad,NgbCell,Grids);
        GradNgb=LSEGrad[NgbCell];
        W=Weight(Grids,PCell,FaceNum);
        GradF=Sum(SVP((W),GradP),SVP((1.0-W),GradNgb));
//        DensF=0.5*(DensCal(Grids->Variables.Alpha[PCell])+DensCal(Grids->Variables.Alpha[NgbCell]));
        if (fabs((maxDens-minDens)-fabs(DensW-DensA))<0.01*maxDens)
            DensF=0.5*(minDens+maxDens);
        else
            DensF=minDens;
        NuF=1.0/DensF;
        D=Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[NgbCell].CellCenter));
        S=Grids->Cells[PCell].Area[FaceNum];
//Minimum Correction Approach
//            Delta=SVP(Dot(D,S)/Dot(D,D),D);
//Orthogonal Correction Approach
//            Delta=SVP(sqrt(Dot(S,S)/Dot(D,D)),D);
//Over-relaxed Approach
        Delta=SVP(Dot(S,S)/Dot(D,S),D);
        K=Sum(S,MinusVec(Delta));
        Amag=Mag(S);
        *aP += -(NuF*sqrt(Dot(Delta,Delta)/Dot(D,D)))/Amag;
        *aN  =  (NuF*sqrt(Dot(Delta,Delta)/Dot(D,D)))/Amag;
        *ST +=  (NuF*Dot(K,GradF))/Amag;
    }
}

double FSElevationSet(double x,double y,double t)
{
  double z;
  double CWE=0.0;  //Calm water elavation is set (Position of free surface in global coordination)
//Calm Water Free surface is set as below
  z=CWE;
  return(z);
}

struct MyVector FlowFunction(struct MyVector Pos,double t)
{
  struct MyVector NewVel;
  double x,y,z;
  x=Pos.X;y=Pos.Y,z=Pos.Z;
  NewVel.X=1.0;
  NewVel.Y=0.0;// 2.0*sin(pi*z)*cos(pi*z)*sqr(sin(pi*y));
  NewVel.Z=1.0e-10;//-2.0*sin(pi*y)*cos(pi*y)*sqr(sin(pi*z));

  return(NewVel);
}

void CopyToLaspackMatrix(QMatrix *LaspackMC,
                         struct MatrixCoefficient *SMC,
                         struct Grid *Grids)
{
    int i,s,PCell,CellNo,BdrSum;
    double Val;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        BdrSum=0;
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for (s=0;s<NoOfFace;s++)
            if (Grids->Cells[PCell].Neighbor[s]<0)
                BdrSum++;
        Q_SetLen(LaspackMC,PCell+1,NoOfFace+1-BdrSum);
        i=0;
        Q_SetEntry(LaspackMC,PCell+1,i++,PCell+1,SMC->Elem[PCell][0]);
        for(s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
            {
                CellNo = Grids->Cells[PCell].Neighbor[s];
                Q_SetEntry(LaspackMC,PCell+1,i++,CellNo+1,SMC->Elem[PCell][s+1]);
            }
    }
}

void CopyToLaspackVector(Vector *GlobU,double *U,struct Grid *Grids)
{
    int PCell;
    for (PCell=0;PCell<Grids->NC;PCell++)
        V_SetCmp(GlobU,PCell+1,U[PCell]);
}

void CopyFromLaspackVector(double *U,Vector *GlobU,struct Grid *Grids)
{
    int PCell;
    for (PCell=0;PCell<Grids->NC;PCell++)
        U[PCell]=V_GetCmp(GlobU,PCell+1);
}

void CopyFromLaspackVectorByOffset(double *U,Vector *GlobU,double Offset,struct Grid *Grids)
{
    int PCell;
    for (PCell=0;PCell<Grids->NC;PCell++)
        U[PCell]=V_GetCmp(GlobU,PCell+1)-Offset;
}

int isInNode(struct PosValueList *PV)
{
    int result=1,i;
    for(i=1;i<PV->size;i++)
        if (Mag(Minus(PV->Elem[0].Pos,PV->Elem[i].Pos))>geps)
        {
            result = 0;
            break;
        }
    return result;
}

int isInLine(struct PosValueList *PV,struct MyVector *LL)
{

// Line Equation :
//   (X-X1)   (Y-Y1)   (Z-Z1)
//   ------ = ------ = ------ = Landa
//     LX       LY       LZ

    struct MyVector P,P0,L;
    int result=1,findNode=0,i;

    P0 = PV->Elem[0].Pos;
    for(i=1;i<PV->size;i++)
        if (Mag(Minus(P0,PV->Elem[i].Pos))>geps)
        {
            findNode=1;
            L = Normalized(Minus(P0,PV->Elem[i].Pos));
            break;
        }
    LL->X=L.X;
    LL->Y=L.Y;
    LL->Z=L.Z;

    for(i=1;i<PV->size;i++)
    {
        P=PV->Elem[i].Pos;
        if((fabs(L.X*(P.Y-P0.Y)-L.Y*(P.X-P0.X))+
            fabs(L.X*(P.Z-P0.Z)-L.Z*(P.X-P0.X))+
            fabs(L.Y*(P.Z-P0.Z)-L.Z*(P.Y-P0.Y)))>1e-3)
        {
            result = 0;
            break;
        }
    }
    return (result&&findNode);
}

int isInSurface(struct PosValueList *PV,struct MyVector *NN)
{
// Surface Equation :
// Nx * (X-X0) + Ny * (Y-Y0) + Nz * (Z-Z0) = 0

    struct MyVector P,P0,P1,P2,N;
    int result=1,i;
    int firstFound=0,secondFound=0;
    P0 = PV->Elem[0].Pos;

    for(i=1;i<PV->size;i++)
    {
        if (!firstFound)
        {
            if (Mag(Minus(P0,PV->Elem[i].Pos))>geps)
            {
                P1 = PV->Elem[i].Pos;
                firstFound=1;
            };
        }
        else if (Mag(Cross(Minus(P0,P1),Minus(P0,PV->Elem[i].Pos)))>geps)
        {
            P2 = PV->Elem[i].Pos;
            secondFound=1;
            break;
        };
    }

    if (firstFound&&secondFound)
    {
        N =Normalized(Cross(Minus(P0,P1),Minus(P0,P2)));

        NN->X=N.X;
        NN->Y=N.Y;
        NN->Z=N.Z;
        for(i=1;i<PV->size;i++)
        {
            P=PV->Elem[i].Pos;
            if (fabs(N.X*(P.X-P0.X)+
                     N.Y*(P.Y-P0.Y)+
                     N.Z*(P.Z-P0.Z))>1e-3)
            {
                result = 0;
                break;
            }
        }
    }
    return (result&&firstFound&&secondFound);
}

void ShiftAdd(struct PosValueList *PV,struct MyVector L)
{
    int i;
    struct PosValueList result;
    result.size=PV->size*2;
    result.Elem=(struct PosValue*)malloc(sizeof(struct PosValue)*result.size);
    for (i=0;i<PV->size;i++)
    {
        result.Elem[i].Pos=PV->Elem[i].Pos;
        result.Elem[i+PV->size].Pos=Sum(PV->Elem[i].Pos,L);
        result.Elem[i].Value=result.Elem[i+PV->size].Value=PV->Elem[i].Value;
        result.Elem[i].Weight=result.Elem[i+PV->size].Weight=PV->Elem[i].Weight;
    }
    PV->size=result.size;
    free(PV->Elem);
    PV->Elem=(struct PosValue*)malloc(sizeof(struct PosValue)*PV->size);

    for (i=0;i<PV->size;i++)
    {
        PV->Elem[i].Pos=result.Elem[i].Pos;
        PV->Elem[i].Value=result.Elem[i].Value;
        PV->Elem[i].Weight=result.Elem[i].Weight;
    }

    free(result.Elem);

}

/*
struct MyVector getLSEGrad3D(struct PosValueList *PV)
{
    double x2bar=0.0,xybar=0.0,xzbar=0.0,xbar=0.0,
           y2bar=0.0,yzbar=0.0,ybar=0.0,
           z2bar=0.0,zbar=0.0,
           xfeebar=0.0,yfeebar=0.0,zfeebar=0.0,feebar=0.0,
           xbar2=0.0,ybar2=0.0,zbar2=0.0,
           xzbar2=0.0, xybar2=0.0, yzbar2=0.0,
           denominator=0.0,num;
    struct MyVector result;
    int n;
    num = (double)PV->size;
    for (n=0;n<PV->size;n++)
    {
        xbar+=PV->Elem[n].Pos.X;
        ybar+=PV->Elem[n].Pos.Y;
        zbar+=PV->Elem[n].Pos.Z;

        x2bar+=sqr(PV->Elem[n].Pos.X);
        y2bar+=sqr(PV->Elem[n].Pos.Y);
        z2bar+=sqr(PV->Elem[n].Pos.Z);

        xybar+=PV->Elem[n].Pos.X*PV->Elem[n].Pos.Y;
        xzbar+=PV->Elem[n].Pos.X*PV->Elem[n].Pos.Z;
        yzbar+=PV->Elem[n].Pos.Y*PV->Elem[n].Pos.Z;

        feebar+=PV->Elem[n].Value;
        xfeebar+=PV->Elem[n].Pos.X*PV->Elem[n].Value;
        yfeebar+=PV->Elem[n].Pos.Y*PV->Elem[n].Value;
        zfeebar+=PV->Elem[n].Pos.Z*PV->Elem[n].Value;
    }

    xbar/=num;
    ybar/=num;
    zbar/=num;

    x2bar/=num;
    y2bar/=num;
    z2bar/=num;

    xybar/=num;
    xzbar/=num;
    yzbar/=num;

    feebar/=num;
    xfeebar/=num;
    yfeebar/=num;
    zfeebar/=num;

    xbar2=sqr(xbar);
    ybar2=sqr(ybar);
    zbar2=sqr(zbar);

    xzbar2=sqr(xzbar);
    xybar2=sqr(xybar);
    yzbar2=sqr(yzbar);

    denominator =(  x2bar*y2bar*z2bar - x2bar*y2bar*zbar2 - x2bar*yzbar2 + 2.0*x2bar*yzbar*ybar*zbar - x2bar*ybar2*z2bar -
                     xybar2*z2bar + xybar2*zbar2 + 2.0*xybar*yzbar*xzbar - 2.0*xybar*yzbar*xbar*zbar -
                     2.0*xybar*ybar*xzbar*zbar + 2.0*xybar*ybar*xbar*z2bar - y2bar*xzbar2 + 2.0*xzbar*y2bar*xbar*zbar +
                     xzbar2*ybar2 - 2.0*xzbar*ybar*xbar*yzbar - y2bar*xbar2*z2bar + xbar2*yzbar2);
    result.X = - ( - xfeebar*y2bar*z2bar + xfeebar*y2bar*zbar2 + xfeebar*yzbar2 - 2.0*xfeebar*yzbar*ybar*zbar +
                     xfeebar*ybar2*z2bar + yfeebar*xybar*z2bar - yfeebar*xybar*zbar2 - yfeebar*yzbar*xzbar +
                     yfeebar*yzbar*xbar*zbar + yfeebar*ybar*xzbar*zbar - yfeebar*ybar*xbar*z2bar - zfeebar*xybar*yzbar +
                     zfeebar*xybar*ybar*zbar + zfeebar*y2bar*xzbar - zfeebar*y2bar*xbar*zbar - zfeebar*xzbar*ybar2 +
                     zfeebar*ybar*xbar*yzbar + feebar*xybar*yzbar*zbar - feebar*xybar*ybar*z2bar - feebar*y2bar*xzbar*zbar +
                     feebar*y2bar*xbar*z2bar + feebar*yzbar*xzbar*ybar - feebar*xbar*yzbar2) / denominator;
    result.Y = - (   xfeebar*xybar*z2bar - xfeebar*xybar*zbar2 - xfeebar*yzbar*xzbar + xfeebar*yzbar*xbar*zbar +
                     xfeebar*ybar*xzbar*zbar - xfeebar*ybar*xbar*z2bar - yfeebar*x2bar*z2bar + yfeebar*x2bar*zbar2 +
                     yfeebar*xzbar2 - 2.0*yfeebar*xzbar*xbar*zbar + yfeebar*xbar2*z2bar + zfeebar*x2bar*yzbar -
                     zfeebar*x2bar*ybar*zbar - zfeebar*xzbar*xybar + zfeebar*xzbar*xbar*ybar + zfeebar*xbar*xybar*zbar -
                     zfeebar*xbar2*yzbar - feebar*x2bar*yzbar*zbar + feebar*x2bar*ybar*z2bar + feebar*xzbar*xybar*zbar -
                     feebar*xzbar2*ybar - feebar*xbar*xybar*z2bar + feebar*xbar*yzbar*xzbar) / denominator;
    result.Z = - ( - xfeebar*xybar*yzbar + xfeebar*xybar*ybar*zbar + xfeebar*y2bar*xzbar - xfeebar*y2bar*xbar*zbar -
                     xfeebar*xzbar*ybar2 + xfeebar*ybar*xbar*yzbar + yfeebar*x2bar*yzbar - yfeebar*x2bar*ybar*zbar -
                     yfeebar*xzbar*xybar + yfeebar*xzbar*xbar*ybar + yfeebar*xbar*xybar*zbar - yfeebar*xbar2*yzbar -
                     zfeebar*x2bar*y2bar + zfeebar*x2bar*ybar2 + zfeebar*xybar2 - 2.0*zfeebar*xybar*xbar*ybar +
                     zfeebar*xbar2*y2bar + feebar*x2bar*y2bar*zbar - feebar*x2bar*ybar*yzbar - feebar*xybar2*zbar +
                     feebar*xybar*xzbar*ybar + feebar*xbar*xybar*yzbar - feebar*xbar*y2bar*xzbar) / denominator;
    return result;
}
*/

struct MyVector getLSEGrad3D(struct PosValueList *PV)
{
    double x2bar=0.0,xybar=0.0,xzbar=0.0,xbar=0.0,
           y2bar=0.0,yzbar=0.0,ybar=0.0,
           z2bar=0.0,zbar=0.0,
           xfeebar=0.0,yfeebar=0.0,zfeebar=0.0,feebar=0.0,
           xbar2=0.0,ybar2=0.0,zbar2=0.0,
           xzbar2=0.0, xybar2=0.0, yzbar2=0.0,
           denominator=0.0,num=0.0;
    struct MyVector result;
    int n;
    for (n=0;n<PV->size;n++)
    {
        num+=PV->Elem[n].Weight;

        xbar+=PV->Elem[n].Weight*PV->Elem[n].Pos.X;
        ybar+=PV->Elem[n].Weight*PV->Elem[n].Pos.Y;
        zbar+=PV->Elem[n].Weight*PV->Elem[n].Pos.Z;

        x2bar+=PV->Elem[n].Weight*sqr(PV->Elem[n].Pos.X);
        y2bar+=PV->Elem[n].Weight*sqr(PV->Elem[n].Pos.Y);
        z2bar+=PV->Elem[n].Weight*sqr(PV->Elem[n].Pos.Z);

        xybar+=PV->Elem[n].Weight*PV->Elem[n].Pos.X*PV->Elem[n].Pos.Y;
        xzbar+=PV->Elem[n].Weight*PV->Elem[n].Pos.X*PV->Elem[n].Pos.Z;
        yzbar+=PV->Elem[n].Weight*PV->Elem[n].Pos.Y*PV->Elem[n].Pos.Z;

        feebar+=PV->Elem[n].Weight*PV->Elem[n].Value;
        xfeebar+=PV->Elem[n].Weight*PV->Elem[n].Pos.X*PV->Elem[n].Value;
        yfeebar+=PV->Elem[n].Weight*PV->Elem[n].Pos.Y*PV->Elem[n].Value;
        zfeebar+=PV->Elem[n].Weight*PV->Elem[n].Pos.Z*PV->Elem[n].Value;
    }

    xbar/=num;
    ybar/=num;
    zbar/=num;

    x2bar/=num;
    y2bar/=num;
    z2bar/=num;

    xybar/=num;
    xzbar/=num;
    yzbar/=num;

    feebar/=num;
    xfeebar/=num;
    yfeebar/=num;
    zfeebar/=num;

    xbar2=sqr(xbar);
    ybar2=sqr(ybar);
    zbar2=sqr(zbar);

    xzbar2=sqr(xzbar);
    xybar2=sqr(xybar);
    yzbar2=sqr(yzbar);

    denominator =(  x2bar*y2bar*z2bar - x2bar*y2bar*zbar2 - x2bar*yzbar2 + 2.0*x2bar*yzbar*ybar*zbar - x2bar*ybar2*z2bar -
                     xybar2*z2bar + xybar2*zbar2 + 2.0*xybar*yzbar*xzbar - 2.0*xybar*yzbar*xbar*zbar -
                     2.0*xybar*ybar*xzbar*zbar + 2.0*xybar*ybar*xbar*z2bar - y2bar*xzbar2 + 2.0*xzbar*y2bar*xbar*zbar +
                     xzbar2*ybar2 - 2.0*xzbar*ybar*xbar*yzbar - y2bar*xbar2*z2bar + xbar2*yzbar2);

    result.X = - ( - xfeebar*y2bar*z2bar + xfeebar*y2bar*zbar2 + xfeebar*yzbar2 - 2.0*xfeebar*yzbar*ybar*zbar +
                     xfeebar*ybar2*z2bar + yfeebar*xybar*z2bar - yfeebar*xybar*zbar2 - yfeebar*yzbar*xzbar +
                     yfeebar*yzbar*xbar*zbar + yfeebar*ybar*xzbar*zbar - yfeebar*ybar*xbar*z2bar - zfeebar*xybar*yzbar +
                     zfeebar*xybar*ybar*zbar + zfeebar*y2bar*xzbar - zfeebar*y2bar*xbar*zbar - zfeebar*xzbar*ybar2 +
                     zfeebar*ybar*xbar*yzbar + feebar*xybar*yzbar*zbar - feebar*xybar*ybar*z2bar - feebar*y2bar*xzbar*zbar +
                     feebar*y2bar*xbar*z2bar + feebar*yzbar*xzbar*ybar - feebar*xbar*yzbar2) / denominator;
    result.Y = - (   xfeebar*xybar*z2bar - xfeebar*xybar*zbar2 - xfeebar*yzbar*xzbar + xfeebar*yzbar*xbar*zbar +
                     xfeebar*ybar*xzbar*zbar - xfeebar*ybar*xbar*z2bar - yfeebar*x2bar*z2bar + yfeebar*x2bar*zbar2 +
                     yfeebar*xzbar2 - 2.0*yfeebar*xzbar*xbar*zbar + yfeebar*xbar2*z2bar + zfeebar*x2bar*yzbar -
                     zfeebar*x2bar*ybar*zbar - zfeebar*xzbar*xybar + zfeebar*xzbar*xbar*ybar + zfeebar*xbar*xybar*zbar -
                     zfeebar*xbar2*yzbar - feebar*x2bar*yzbar*zbar + feebar*x2bar*ybar*z2bar + feebar*xzbar*xybar*zbar -
                     feebar*xzbar2*ybar - feebar*xbar*xybar*z2bar + feebar*xbar*yzbar*xzbar) / denominator;
    result.Z = - ( - xfeebar*xybar*yzbar + xfeebar*xybar*ybar*zbar + xfeebar*y2bar*xzbar - xfeebar*y2bar*xbar*zbar -
                     xfeebar*xzbar*ybar2 + xfeebar*ybar*xbar*yzbar + yfeebar*x2bar*yzbar - yfeebar*x2bar*ybar*zbar -
                     yfeebar*xzbar*xybar + yfeebar*xzbar*xbar*ybar + yfeebar*xbar*xybar*zbar - yfeebar*xbar2*yzbar -
                     zfeebar*x2bar*y2bar + zfeebar*x2bar*ybar2 + zfeebar*xybar2 - 2.0*zfeebar*xybar*xbar*ybar +
                     zfeebar*xbar2*y2bar + feebar*x2bar*y2bar*zbar - feebar*x2bar*ybar*yzbar - feebar*xybar2*zbar +
                     feebar*xybar*xzbar*ybar + feebar*xbar*xybar*yzbar - feebar*xbar*y2bar*xzbar) / denominator;
//!Debug
//    if (denominator<1e-10)
//    {
//        printf("there is a Bug!\n");
//        exit(-1);
//    }
//!Debug
    if (Dimension=='2')
        result.Y=0.0;
    return result;
}


struct MyVector getLSEGrad1D(struct PosValueList *PV,struct MyVector L)
{
    struct MyVector result,N2L;
    N2L.X=N2L.Y=N2L.Z=0.0;
    // N2L is not a normal to L vector
    if (fabs(L.X)<min(fabs(L.Y),fabs(L.Z))) N2L.X=1.0;
    else if (fabs(L.Y)<fabs(L.Z)) N2L.Y=1.0;
    else N2L.Z=1.0;
    ShiftAdd(PV,N2L);
    ShiftAdd(PV,Cross(N2L,L));
    result = getLSEGrad3D (PV);
    return result;
}

struct MyVector getLSEGrad2D(struct PosValueList *PV,struct MyVector N)
{
    struct MyVector result;
    ShiftAdd(PV,N);
    result = getLSEGrad3D (PV);
    return result;
}

struct MyVector MakeLSEGradient(struct PosValueList *PV)
{
    struct MyVector result,L;

    if (isInNode(PV))
        result.X=result.Y=result.Z=0.0;
    else if (isInLine(PV,&L))
        result = getLSEGrad1D(PV,L);
    else if (isInSurface(PV,&L))
        result = getLSEGrad2D(PV,L);
    else
        result = getLSEGrad3D(PV);
    if (Dimension=='2')
        result.Y=0.0;
    return result;
}


struct MyVector LSEGradient(double *fee,int PCell,struct Grid *Grids)
{
    struct PosValueList PV;
    int list[500],counter=0;
    int FaceNum,j,NgbNo,
        NodeNumber,n;
    struct MyVector result;

    PV.size=0;

#ifdef LSE_NodeShare
    for (n=0;n<Grids->RCells[PCell].NodePerCell;n++)
    {
        NodeNumber=Grids->RCells[PCell].NodeList[n];
        for(j=0;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
        {
            NgbNo=Grids->Nodes[NodeNumber].SharedCells[j];
#endif
#ifdef LSE_Neighbour
    for (j=0;j<Grids->Cells[PCell].NoOfFace;j++)
    {
        if (Grids->Cells[PCell].Neighbor[j]>=0)
        {
            NgbNo=Grids->Cells[PCell].Neighbor[j];
#endif

            list[counter]=NgbNo;
            PV.size++;
            counter++;
        }
    }
    list[counter]=-10;

    PV.size++;
    PV.Elem=(struct PosValue*)malloc(sizeof(struct PosValue)*PV.size);
    n=0;counter=0;
    while (list[n]!=-10)
    {
        if (list[n]!=-1)
        {
            PV.Elem[counter].Pos=Grids->Cells[list[n]].CellCenter;
            PV.Elem[counter].Value=fee[list[n]];
            PV.Elem[counter].Weight=1.0;
            counter++;
        }
        n++;
    }
    PV.Elem[PV.size-1].Pos=Grids->Cells[PCell].CellCenter;
    PV.Elem[PV.size-1].Value=fee[PCell];
    PV.Elem[PV.size-1].Weight=1.0;

    result = MakeLSEGradient(&PV);
    free(PV.Elem);
    if (Dimension=='2')
        result.Y=0.0;
    return result;
}


struct MyVector LSEGradient_Disc(double *fee,int PCell,struct Grid *Grids)
{
    struct PosValueList PV;
    int list[500],counter=0;
    int FaceNum,j,NgbNo,
        NodeNumber,n;
    struct MyVector result;
    double MaxDens,MinDens,gamma=0.01;
    int LightDensBase=0,HeavyDensBase=0,MidDensBase=0,isWeighted=0;
    MaxDens=max(DensW,DensA);
    MinDens=min(DensW,DensA);
    PV.size=0;


    if (fabs(DensCal(Grids->Variables.Alpha[PCell]))<5.0)
        LightDensBase=1;
    else if (fabs(DensCal(Grids->Variables.Alpha[PCell]))>900.0)
        HeavyDensBase=1;
    else
        MidDensBase=1;

#ifdef LSE_NodeShare
    for (n=0;n<Grids->RCells[PCell].NodePerCell;n++)
    {
        NodeNumber=Grids->RCells[PCell].NodeList[n];
        for(j=0;j<Grids->Nodes[NodeNumber].NumOfSharedCells;j++)
        {
            NgbNo=Grids->Nodes[NodeNumber].SharedCells[j];
#endif
#ifdef LSE_Neighbour
    for (j=0;j<Grids->Cells[PCell].NoOfFace;j++)
    {
        if (Grids->Cells[PCell].Neighbor[j]>=0)
        {
            NgbNo=Grids->Cells[PCell].Neighbor[j];
#endif
            if (LightDensBase)
            {
                if (fabs(DensCal(Grids->Variables.Alpha[NgbNo])-DensCal(Grids->Variables.Alpha[PCell]))<5.1)
                {
                    list[counter]=NgbNo;
                    PV.size++;
                }
                else
                    list[counter]=-1;
            }
            else if (HeavyDensBase)
            {
                if (fabs(DensCal(Grids->Variables.Alpha[NgbNo])-DensCal(Grids->Variables.Alpha[PCell]))<100.1)
                {
                    list[counter]=NgbNo;
                    PV.size++;
                }
                else
                    list[counter]=-1;
            }
            else
            {
//                if (fabs(DensCal(Grids->Variables.Alpha[NgbNo]))>1.1)
//                {
//                    isWeighted=1;
//                    list[counter]=NgbNo;
//                    PV.size++;
//                }
//                else
                    list[counter]=-1;
            }

            counter++;
        }
    }
    list[counter]=-10;

    PV.size++;
    PV.Elem=(struct PosValue*)malloc(sizeof(struct PosValue)*PV.size);
    n=0;counter=0;
    if (isWeighted)
    {
        while (list[n]!=-10)
        {
            if (list[n]!=-1)
            {
                PV.Elem[counter].Pos=Grids->Cells[list[n]].CellCenter;
                PV.Elem[counter].Value=fee[list[n]];
                PV.Elem[counter].Weight=1.0/DensCal(Grids->Variables.Alpha[list[n]]);
                counter++;
            }
            n++;
        }
        PV.Elem[PV.size-1].Pos=Grids->Cells[PCell].CellCenter;
        PV.Elem[PV.size-1].Value=fee[PCell];
        PV.Elem[PV.size-1].Weight=1.0/DensCal(Grids->Variables.Alpha[PCell]);
    }
    else
    {
        while (list[n]!=-10)
        {
            if (list[n]!=-1)
            {
                PV.Elem[counter].Pos=Grids->Cells[list[n]].CellCenter;
                PV.Elem[counter].Value=fee[list[n]];
                PV.Elem[counter].Weight=1.0;
                counter++;
            }
            n++;
        }
        PV.Elem[PV.size-1].Pos=Grids->Cells[PCell].CellCenter;
        PV.Elem[PV.size-1].Value=fee[PCell];
        PV.Elem[PV.size-1].Weight=1.0;
    }


    result = MakeLSEGradient(&PV);
    free(PV.Elem);
    if (Dimension=='2')
        result.Y=0.0;
    return result;
}
