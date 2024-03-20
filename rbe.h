void InitialConditionForRigidBodyVariables(struct Grid *Grids)
{
    Grids->DynaVariables.Mass=1.0;
    Grids->DynaVariables.Inertia[0][0]=1.0; Grids->DynaVariables.Inertia[0][1]=0.0; Grids->DynaVariables.Inertia[0][2]=0.0;
    Grids->DynaVariables.Inertia[1][0]=0.0; Grids->DynaVariables.Inertia[1][1]=1.0; Grids->DynaVariables.Inertia[1][2]=0.0;
    Grids->DynaVariables.Inertia[2][0]=0.0; Grids->DynaVariables.Inertia[2][1]=0.0; Grids->DynaVariables.Inertia[2][2]=1.0;

    Grids->DynaVariables.NewLinAccel.X  =Grids->DynaVariables.LinAccel.X  =0.0;
    Grids->DynaVariables.NewLinAccel.Y  =Grids->DynaVariables.LinAccel.Y  =0.0;
    Grids->DynaVariables.NewLinAccel.Z  =Grids->DynaVariables.LinAccel.Z  =0.0;

    Grids->DynaVariables.NewAngAccel.X  =Grids->DynaVariables.AngAccel.X  =0.0;
    Grids->DynaVariables.NewAngAccel.Y  =Grids->DynaVariables.AngAccel.Y  =0.0;
    Grids->DynaVariables.NewAngAccel.Z  =Grids->DynaVariables.AngAccel.Z  =0.0;

    Grids->DynaVariables.NewLinVel.X    =Grids->DynaVariables.LinVel.X    =0.0;
    Grids->DynaVariables.NewLinVel.Y    =Grids->DynaVariables.LinVel.Y    =0.0;
    Grids->DynaVariables.NewLinVel.Z    =Grids->DynaVariables.LinVel.Z    =0.0;

    Grids->DynaVariables.NewAngVel.X    =Grids->DynaVariables.AngVel.X    =0.0;
    Grids->DynaVariables.NewAngVel.Y    =Grids->DynaVariables.AngVel.Y    =0.0;
    Grids->DynaVariables.NewAngVel.Z    =Grids->DynaVariables.AngVel.Z    =0.0;

    Grids->DynaVariables.NewMassCenter.X=Grids->DynaVariables.MassCenter.X=0.0;
    Grids->DynaVariables.NewMassCenter.Y=Grids->DynaVariables.MassCenter.Y=0.0;
    Grids->DynaVariables.NewMassCenter.Z=Grids->DynaVariables.MassCenter.Z=0.0;

    Grids->DynaVariables.NewBodyOr.X    =Grids->DynaVariables.BodyOr.X    =0.0;
    Grids->DynaVariables.NewBodyOr.Y    =Grids->DynaVariables.BodyOr.Y    =0.0;
    Grids->DynaVariables.NewBodyOr.Z    =Grids->DynaVariables.BodyOr.Z    =0.0;
}

struct MyVector CalcPressureForces(struct Grid *Grids, int GroupNo)
{
    int FaceNum,PCell,k,i;
    int c;
    double Pf;
    struct MyVector FaceCenter,Area,C2F;
    struct MyVector F;
    F.X=F.Y=F.Z=0.0;
    c=0;
    if (GroupNo==0)
        c=0;
    else
    {
        for (i=0;i<GroupNo;i++)
            c+=Grids->RBNo[i];
    }
//    printf("c = %d\n",c);
//    getchar();
    for(k=c;k<c+Grids->RBNo[GroupNo];k++)
    {
            PCell=Grids->RigidBoundaryCells[k].CellNum;
            FaceNum=Grids->RigidBoundaryCells[k].FaceNum;
//            printf("cell is: %d   & faceNum is: %d\n",PCell,FaceNum);
            Area=Grids->Cells[PCell].Area[FaceNum];
            FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
            C2F=Sum(FaceCenter,MinusVec(Grids->Cells[PCell].CellCenter));
            Pf=Grids->Variables.P[PCell]+Dot(Grids->Variables.PressLSEGrad[PCell],C2F);
            F=Sum(F,SVP(Pf,Area));
    }
    return F;
}

struct MyVector CalcPressureMoments(struct Grid *Grids, int GroupNo)
{
    int FaceNum,PCell,k,c,i;
    double Pf;
    struct MyVector FaceCenter,Area,C2F;
    struct MyVector M;
    M.X=M.Y=M.Z=0.0;
    c=0;
    if (GroupNo==0)
        c=0;
    else
    {
        for (i=0;i<GroupNo;i++)
            c+=Grids->RBNo[i];
    }
    for(k=c;k<c+Grids->RBNo[GroupNo];k++)
    {
        PCell=Grids->RigidBoundaryCells[k].CellNum;
        FaceNum=Grids->RigidBoundaryCells[k].FaceNum;
        Area=Grids->Cells[PCell].Area[FaceNum];
        FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
        C2F=Sum(FaceCenter,MinusVec(Grids->Cells[PCell].CellCenter));
        Pf=Grids->Variables.NP[PCell]+Dot(Grids->Variables.PressLSEGrad[PCell],C2F);
        M=Sum(M,Cross(Sum(MinusVec(Grids->DynaVariables.MassCenter),FaceCenter),SVP(Pf,Area)));
    }
    return M;
}

struct MyVector CalcViscousForces(struct Grid *Grids, int GroupNo)
{
    int FaceNum,PCell,k,c,i;
    double dn,Amag,mu;
    struct MyVector FaceCenter,Area,n,DV,C2F;
    struct MyVector F;
    F.X=F.Y=F.Z=0.0;
    c=0;
    if (GroupNo==0)
        c=0;
    else
    {
        for (i=0;i<GroupNo;i++)
            c+=Grids->RBNo[i];
    }
    for(k=c;k<c+Grids->RBNo[GroupNo];k++)
    {
        PCell=Grids->RigidBoundaryCells[k].CellNum;
        FaceNum=Grids->RigidBoundaryCells[k].FaceNum;
        Area=Grids->Cells[PCell].Area[FaceNum];
        FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
        Amag=sqrt(Dot(Area,Area));
        n=SVP((1./Amag),Area);
        dn=fabs((Dot(Sum(Grids->Cells[PCell].CellCenter,MinusVec(FaceCenter)),n)));
        mu=ViscCal(Grids->Variables.Alpha[PCell]);
        DV.X=Grids->Variables.NU[PCell];
        DV.Y=Grids->Variables.NV[PCell];
        DV.Z=Grids->Variables.NW[PCell];
        struct MyVector t,VP;
        VP=Sum(DV,MinusVec(SVP(Dot(DV,n),n)));
        t=SVP(1./(sqrt(Dot(VP,VP))+geps),VP);
        F=Sum(F,SVP((mu*Amag/dn*Dot(DV,t)),t));
    }
    return F;
}

struct MyVector CalcViscousMoments(struct Grid *Grids, int GroupNo)
{
    int FaceNum,PCell,k,c,i;
    double dn,Amag,mu;
    struct MyVector FaceCenter,Area,n,DV,C2F;
    struct MyVector M;
    M.X=M.Y=M.Z=0.0;
    c=0;
    if (GroupNo==0)
        c=0;
    else
    {
        for (i=0;i<GroupNo;i++)
            c+=Grids->RBNo[i];
    }
    for(k=c;k<c+Grids->RBNo[GroupNo];k++)
    {
        PCell=Grids->RigidBoundaryCells[k].CellNum;
        FaceNum=Grids->RigidBoundaryCells[k].FaceNum;
        Area=Grids->Cells[PCell].Area[FaceNum];
        FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];

        Amag=sqrt(Dot(Area,Area));
        n=SVP((1./Amag),Area);

        dn=fabs((Dot(Sum(Grids->Cells[PCell].CellCenter,MinusVec(FaceCenter)),n)));
        mu=ViscCal(Grids->Variables.Alpha[PCell]);
        DV.X=Grids->Variables.U[PCell];
        DV.Y=Grids->Variables.V[PCell];
        DV.Z=Grids->Variables.W[PCell];
        struct MyVector t,VP;
        VP=Sum(DV,MinusVec(SVP(Dot(DV,n),n)));
        t=SVP(1./sqrt(Dot(VP,VP)+geps),VP);
        M=Sum(M,Cross(Sum(MinusVec(Grids->DynaVariables.MassCenter),FaceCenter),SVP((mu*Amag/dn*Dot(DV,t)),t)));

    }
    return M;
}

void InvMatrix(double InvA[3][3],double A[3][3])
{
    double Det;
    Det=(-A[0][0]*A[1][1]*A[2][2]+A[0][0]*A[2][1]*A[1][2]+A[1][0]*A[0][1]*A[2][2]-A[1][0]*A[2][1]*A[0][2]-A[2][0]*A[0][1]*A[1][2]+A[2][0]*A[1][1]*A[0][2]);
    InvA[0][0]=-(A[1][1]*A[2][2]-A[2][1]*A[1][2])/Det;
    InvA[0][1]= (A[0][1]*A[2][2]-A[2][1]*A[0][2])/Det;
    InvA[0][2]=-(A[0][1]*A[1][2]-A[1][1]*A[0][2])/Det;
    InvA[1][0]= (A[1][0]*A[2][2]-A[2][0]*A[1][2])/Det;
    InvA[1][1]=-(A[0][0]*A[2][2]-A[2][0]*A[0][2])/Det;
    InvA[1][2]= (A[0][0]*A[1][2]-A[1][0]*A[0][2])/Det;
    InvA[2][0]=-(A[1][0]*A[2][1]-A[2][0]*A[1][1])/Det;
    InvA[2][1]= (A[0][0]*A[2][1]-A[2][0]*A[0][1])/Det;
    InvA[2][2]=-(A[0][0]*A[1][1]-A[1][0]*A[0][1])/Det;
}

struct MyVector MatrixDotVector(double Matrix[3][3],struct MyVector V)
{
    struct MyVector U;
    U.X=V.X*Matrix[0][0]+V.Y*Matrix[0][1]+V.Z*Matrix[0][2];
    U.Y=V.X*Matrix[1][0]+V.Y*Matrix[1][1]+V.Z*Matrix[1][2];
    U.Z=V.X*Matrix[2][0]+V.Y*Matrix[2][1]+V.Z*Matrix[2][2];
    return U;
}

void CEA(double R[3][3],struct MyVector Or)
{
    double f,t;
    f=Or.X;t=Or.Y;
    R[0][0]=1.0;
    R[0][1]=sin(f)*tan(t);
    R[0][2]=cos(f)*tan(t);
    R[1][0]=0.0;
    R[1][1]=cos(f);
    R[1][2]=-sin(f);
    R[2][0]=0.0;
    R[2][1]=sin(f)/cos(t);
    R[2][2]=cos(f)/cos(t);
}

void DCM(double R[3][3],struct MyVector Or)
{
    double f,t,s;
    f=Or.X;t=Or.Y;s=Or.Z;
    R[0][0]=cos(t)*cos(s);
    R[0][1]=cos(t)*sin(s);
    R[0][2]=-sin(t);
    R[1][0]=sin(f)*sin(t)*cos(s)-cos(f)*sin(s);
    R[1][1]=sin(f)*sin(t)*sin(s)+cos(f)*cos(s);
    R[1][2]=sin(f)*cos(t);
    R[2][0]=cos(f)*sin(t)*cos(s)+sin(f)*sin(s);
    R[2][1]=cos(f)*sin(t)*sin(s)-sin(f)*cos(s);
    R[2][2]=cos(f)*cos(t);
}

void UpdateCellCenters(struct Grid *Grids)
{
    int PCell;
    struct MyVector R,dummyVect,NewPoint,Omega,Vc,MassC;
    double err;

    Vc=Grids->DynaVariables.NewLinVel;
    Omega=Grids->DynaVariables.NewAngVel;
    MassC=Grids->DynaVariables.NewMassCenter;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NewPoint=Grids->Cells[PCell].CellCenter;
//        do
//        {
            dummyVect=NewPoint;
            R=Sum(MinusVec(MassC),NewPoint);
            NewPoint=Sum(Grids->Cells[PCell].CellCenter,SVP(dt,Sum(Vc,Cross(Omega,R))));
            err=sqrt(Dot(Sum(MinusVec(dummyVect),NewPoint),Sum(MinusVec(dummyVect),NewPoint)));
//        }while (err>=geps);
        Grids->Cells[PCell].CellCenter=NewPoint;
    }
}

void UpdateFaceCenters(struct Grid *Grids)
{
    int PCell,FaceNum,NoOfFace;
    struct MyVector R,dummyVect,NewPoint,Omega,Vc,MassC;
    double err;

    Vc=Grids->DynaVariables.NewLinVel;
    Omega=Grids->DynaVariables.NewAngVel;
    MassC=Grids->DynaVariables.NewMassCenter;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(FaceNum=0;FaceNum<NoOfFace;FaceNum++)
        {
            NewPoint=Grids->Cells[PCell].FaceCenter[FaceNum];
//            do
//            {
                dummyVect=NewPoint;
                R=Sum(MinusVec(MassC),NewPoint);
                NewPoint=Sum(Grids->Cells[PCell].FaceCenter[FaceNum],SVP(dt,Sum(Vc,Cross(Omega,R))));
                err=sqrt(Dot(Sum(MinusVec(dummyVect),NewPoint),Sum(MinusVec(dummyVect),NewPoint)));
//            }while (err>=geps);
            Grids->Cells[PCell].FaceCenter[FaceNum]=NewPoint;
        }
    }
}

void UpdateFaceCoincides(struct Grid *Grids)
{
    int PCell,FaceNum,NoOfFace;
    struct MyVector R,dummyVect,NewPoint,Omega,Vc,MassC;
    double err;

    Vc=Grids->DynaVariables.NewLinVel;
    Omega=Grids->DynaVariables.NewAngVel;
    MassC=Grids->DynaVariables.NewMassCenter;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(FaceNum=0;FaceNum<NoOfFace;FaceNum++)
        {
            NewPoint=Grids->Cells[PCell].FaceCoincide[FaceNum];
//            do
//            {
                dummyVect=NewPoint;
                R=Sum(MinusVec(MassC),NewPoint);
                NewPoint=Sum(Grids->Cells[PCell].FaceCoincide[FaceNum],SVP(dt,Sum(Vc,Cross(Omega,R))));
                err=sqrt(Dot(Sum(MinusVec(dummyVect),NewPoint),Sum(MinusVec(dummyVect),NewPoint)));
//            }while (err>=geps);
            Grids->Cells[PCell].FaceCoincide[FaceNum]=NewPoint;
        }
    }
}

void UpdateProjectedCenters(struct Grid *Grids)
{
    int PCell,FaceNum,NoOfFace;
    struct MyVector R,dummyVect,NewPoint,Omega,Vc,MassC;
    double err;

    Vc=Grids->DynaVariables.NewLinVel;
    Omega=Grids->DynaVariables.NewAngVel;
    MassC=Grids->DynaVariables.NewMassCenter;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(FaceNum=0;FaceNum<NoOfFace;FaceNum++)
        {
            NewPoint=Grids->Cells[PCell].ProjectedCenter[FaceNum];
//            do
//            {
                dummyVect=NewPoint;
                R=Sum(MinusVec(MassC),NewPoint);
                NewPoint=Sum(Grids->Cells[PCell].ProjectedCenter[FaceNum],SVP(dt,Sum(Vc,Cross(Omega,R))));
                err=sqrt(Dot(Sum(MinusVec(dummyVect),NewPoint),Sum(MinusVec(dummyVect),NewPoint)));
//            }while (err>=geps);
            Grids->Cells[PCell].ProjectedCenter[FaceNum]=NewPoint;
        }
    }
}

void UpdateNodePositions(struct Grid *Grids)
{
    int NNum;
    struct MyVector R,dummyVect,NewPoint,Omega,Vc,MassC;
    double err;

    Vc=Grids->DynaVariables.NewLinVel;
    Omega=Grids->DynaVariables.NewAngVel;
    MassC=Grids->DynaVariables.NewMassCenter;
    for (NNum=0;NNum<Grids->NN;NNum++)
    {
        NewPoint=Grids->Nodes[NNum].Pos;
//        do
//        {
            dummyVect=NewPoint;
            R=Sum(MinusVec(MassC),NewPoint);
            NewPoint=Sum(Grids->Nodes[NNum].Pos,SVP(dt,Sum(Vc,Cross(Omega,R))));
            err=sqrt(Dot(Sum(MinusVec(dummyVect),NewPoint),Sum(MinusVec(dummyVect),NewPoint)));
//        }while (err>=geps);
        Grids->Nodes[NNum].Pos=NewPoint;
    }
}

void UpdateFaceAreas(struct Grid *Grids)
{
    int PCell,FaceNum,NoOfFace;
    struct MyVector AVect;
    double R[3][3],InvR[3][3];

    DCM(R,Grids->DynaVariables.NewBodyOr);
    InvMatrix(InvR,R);
    DCM(R,Grids->DynaVariables.BodyOr);
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;

        for (FaceNum=0;FaceNum<NoOfFace;FaceNum++)
            Grids->Cells[PCell].Area[FaceNum]=MatrixDotVector(InvR,MatrixDotVector(R,Grids->Cells[PCell].Area[FaceNum]));
    }
}


void AdvanceGrids(struct Grid *Grids)
{
    double R[3][3];

    Grids->DynaVariables.NewMassCenter = Sum(Grids->DynaVariables.MassCenter,SVP(dt,Grids->DynaVariables.NewLinVel));
    CEA(R,Grids->DynaVariables.BodyOr);
    Grids->DynaVariables.NewBodyOr=Sum(SVP(dt,(MatrixDotVector(R,Grids->DynaVariables.NewAngVel))),Grids->DynaVariables.BodyOr);
//    printf("\delta is : ( %e , %e , %e )\n",SVP(dt,Grids->DynaVariables.NewLinVel).X
//                                           ,SVP(dt,Grids->DynaVariables.NewLinVel).Y
//                                           ,SVP(dt,Grids->DynaVariables.NewLinVel).Z);
    UpdateCellCenters(Grids);
    UpdateFaceCenters(Grids);
    UpdateFaceCoincides(Grids);
    UpdateProjectedCenters(Grids);
    UpdateNodePositions(Grids);
    UpdateFaceAreas(Grids);

    Grids->DynaVariables.LinAccel   = Grids->DynaVariables.NewLinAccel;
    Grids->DynaVariables.AngAccel   = Grids->DynaVariables.NewAngAccel;
    Grids->DynaVariables.LinVel     = Grids->DynaVariables.NewLinVel;
    Grids->DynaVariables.AngVel     = Grids->DynaVariables.NewAngVel;
    Grids->DynaVariables.MassCenter = Grids->DynaVariables.NewMassCenter;
    Grids->DynaVariables.BodyOr     = Grids->DynaVariables.NewBodyOr;

    printf("\nMassCenter is : ( %e , %e , %e )\n",Grids->DynaVariables.MassCenter.X,Grids->DynaVariables.MassCenter.Y,Grids->DynaVariables.MassCenter.Z);
    printf("BodyOrient is : ( %e , %e , %e )\n",Grids->DynaVariables.BodyOr.X,Grids->DynaVariables.BodyOr.Y,Grids->DynaVariables.BodyOr.Z);

}

void SixDOFSolver(struct MyVector CF,struct MyVector CM,struct MyVector Force,struct MyVector Moment,struct Grid *Grids)
{
    struct MyVector tteta,teta0,dummyVect;
    double InvI[3][3],R[3][3],InvR[3][3],err;
    double aau=.5 ,lau=.5, aau2=.5 ,lau2=.5;

  //Projection of control force and moment from local CS to global CS
    DCM(R,Grids->DynaVariables.NewBodyOr);//We need to use InvR if control forces wich applide under local coordination, works.
    InvMatrix(InvR,R);
    Force=Sum(MatrixDotVector(InvR,CF),Sum(Force,SVP(Grids->DynaVariables.Mass,Gravity)));
    Moment=Sum(CM,MatrixDotVector(R,Moment));

    Grids->DynaVariables.AngAccel=MatrixDotVector(R,Grids->DynaVariables.AngAccel);
    Grids->DynaVariables.AngVel=MatrixDotVector(R,Grids->DynaVariables.AngVel);
    Grids->DynaVariables.NewAngAccel=Grids->DynaVariables.AngAccel;
    Grids->DynaVariables.NewAngVel=Grids->DynaVariables.AngVel;

    Force.X=0.0;
    Force.Y=0.0;
//    Force.Z=0.0;
    Moment.X=0.0;
    Moment.Y=0.0;
    Moment.Z=0.0;

    InvMatrix(InvI,Grids->DynaVariables.Inertia);

//    aau=lau=0.6;
//    aau2=lau2=0.4;

    do
    {
        dummyVect=Grids->DynaVariables.NewAngAccel;
        Grids->DynaVariables.NewAngAccel=Sum(SVP(aau*aau2,MatrixDotVector(InvI,(Sum(Moment,MinusVec(Cross(Grids->DynaVariables.NewAngVel,(MatrixDotVector(Grids->DynaVariables.Inertia,Grids->DynaVariables.NewAngVel))))))))
                     ,Sum(SVP(aau2*(1.-aau),Grids->DynaVariables.NewAngAccel),SVP((1.-aau2),Grids->DynaVariables.AngAccel)));
        Grids->DynaVariables.NewAngVel=Sum(Grids->DynaVariables.AngVel,SVP(dt,Grids->DynaVariables.NewAngAccel));
        err=sqrt(Dot(Sum(MinusVec(dummyVect),Grids->DynaVariables.NewAngAccel),Sum(MinusVec(dummyVect),Grids->DynaVariables.NewAngAccel)));
    }while (err>=geps);

    Grids->DynaVariables.NewBodyOr=Grids->DynaVariables.BodyOr;

    dummyVect=Grids->DynaVariables.NewBodyOr;
    CEA(R,Grids->DynaVariables.NewBodyOr);
    Grids->DynaVariables.NewBodyOr=Sum(SVP(dt,(MatrixDotVector(R,Grids->DynaVariables.NewAngVel))),Grids->DynaVariables.BodyOr);

    err=sqrt(Dot(Sum(MinusVec(Grids->DynaVariables.NewBodyOr),dummyVect),Sum(MinusVec(Grids->DynaVariables.NewBodyOr),dummyVect)));

//    DCM(R,Grids->DynaVariables.BodyOr);
//    InvMatrix(InvR,R);

//    Grids->DynaVariables.AngAccel=MatrixDotVector(InvR,Grids->DynaVariables.AngAccel);
//    Grids->DynaVariables.AngVel=MatrixDotVector(InvR,Grids->DynaVariables.AngVel);

    DCM(R,Grids->DynaVariables.NewBodyOr);
    InvMatrix(InvR,R);

    Grids->DynaVariables.NewAngAccel=MatrixDotVector(InvR,Grids->DynaVariables.NewAngAccel);
    Grids->DynaVariables.NewAngVel=MatrixDotVector(InvR,Grids->DynaVariables.NewAngVel);

    Grids->DynaVariables.NewLinAccel=Sum(SVP(lau*lau2/Grids->DynaVariables.Mass,Force),Sum(SVP(lau2*(1.-lau),Grids->DynaVariables.NewLinAccel),SVP(1.-lau2,Grids->DynaVariables.LinAccel)));
    Grids->DynaVariables.NewLinVel=Sum(Grids->DynaVariables.LinVel,SVP(dt,Grids->DynaVariables.NewLinAccel));
}

void SolveRBEquation(double *RBEError,struct Grid *Grids)
{
    struct MyVector CF,CM //Control Force and Moment
                   ,dummyVect,FluidForce,FluidMoment;
    double denominator;
    int i;
    *RBEError=0.0;
    dummyVect=Grids->DynaVariables.NewLinAccel;
    denominator=sqrt(Dot(dummyVect,dummyVect))+1e-10;

    CF.X=CF.Y=CF.Z=0.0;
    CM.X=CM.Y=CM.Z=0.0;

    ///mahdy

    for(i=0;i<Grids->NRBG;i++)
    {
        FluidForce=Sum(FluidForce,CalcPressureForces(Grids,i));
        FluidForce=Sum(FluidForce,CalcViscousForces(Grids,i));

        FluidMoment=Sum(FluidMoment,CalcViscousMoments(Grids,i));
        FluidMoment=Sum(FluidMoment,CalcPressureMoments(Grids,i));
    }

    SixDOFSolver(CF,CM,FluidForce,FluidMoment,Grids);
    *RBEError+=sqrt(Dot(Sum(MinusVec(dummyVect),Grids->DynaVariables.NewLinAccel),Sum(MinusVec(dummyVect),Grids->DynaVariables.NewLinAccel)))/denominator;
}
