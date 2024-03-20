//void UpdatePressLSEGradient(struct Grid *Grids)
//{
//    struct MatrixCoefficient GMC;
//    int PCell,NCol;
//    GMC.Elem=(double**) malloc(sizeof(double*)*Grids->NC*3);
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NCol=Grids->Cells[PCell].NoOfFace+1;
//
//        GMC.Elem[PCell*3  ]=(double*)malloc(sizeof(double)*NCol*3);
//        GMC.Elem[PCell*3+1]=(double*)malloc(sizeof(double)*NCol*3);
//        GMC.Elem[PCell*3+2]=(double*)malloc(sizeof(double)*NCol*3);
//    }
//
//    QMatrix LaspackMC;
//    Vector ST,Fee;
//    Q_Constr(&LaspackMC,"LaspackMC",Grids->NC*3,False,Rowws,Normal,True);
//    V_Constr(&ST,"ST",Grids->NC*3,Normal,True);
//    V_Constr(&Fee,"Fee",Grids->NC*3,Normal,True);
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
//    int s=0;
//    double *GradComp,*MyST;
//    GradComp = (double *)malloc(sizeof(double)*Grids->NC*3);
//    MyST = (double *)malloc(sizeof(double)*Grids->NC*3);
//
//    SetSVect(MyST,3*Grids->NC,0.0);
//    SetSVect(GradComp,3*Grids->NC,0.0);
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NCol=Grids->Cells[PCell].NoOfFace+1;
//        SetSVect(GMC.Elem[PCell*3  ],3*NCol,0.0);
//        SetSVect(GMC.Elem[PCell*3+1],3*NCol,0.0);
//        SetSVect(GMC.Elem[PCell*3+2],3*NCol,0.0);
//        GradComp[PCell*3  ]=Grids->Variables.PressLSEGrad[PCell].X;
//        GradComp[PCell*3+1]=Grids->Variables.PressLSEGrad[PCell].Y;
//        GradComp[PCell*3+2]=Grids->Variables.PressLSEGrad[PCell].Z;
//    }
//
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        GMC.Elem[PCell*3  ][0] = 1.0;
//        GMC.Elem[PCell*3+1][1] = 1.0;
//        GMC.Elem[PCell*3+2][2] = 1.0;
//
//        for (s=0;s<Grids->Cells[PCell].NoOfFace;s++)
//        {
//            int Ngb;
//            double Vol,FeeP,FeeNgb;
//            struct MyVector rP,rN,Area;
//
//            Ngb = Grids->Cells[PCell].Neighbor[s];
//            Vol = Grids->Cells[PCell].Volume;
//            Area= Grids->Cells[PCell].Area[s];
//            double aP,aN,CST;
//            GetFaceValueParts_Disc(&aP,&aN,&CST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,s,Grids);
//            rP = Minus(Grids->Cells[PCell].ProjectedCenter[s],Grids->Cells[PCell].CellCenter);
//            FeeP  =Grids->Variables.NP[PCell];
//            FeeNgb=(Ngb>=0)?Grids->Variables.NP[Ngb]:0.0;
//
//            GMC.Elem[PCell*3  ][0]+= -Area.X/Vol*aP*rP.X;
//            GMC.Elem[PCell*3+1][0]+= -Area.Y/Vol*aP*rP.X;
//            GMC.Elem[PCell*3+2][0]+= -Area.Z/Vol*aP*rP.X;
//            GMC.Elem[PCell*3  ][1]+= -Area.X/Vol*aP*rP.Y;
//            GMC.Elem[PCell*3+1][1]+= -Area.Y/Vol*aP*rP.Y;
//            GMC.Elem[PCell*3+2][1]+= -Area.Z/Vol*aP*rP.Y;
//            GMC.Elem[PCell*3  ][2]+= -Area.X/Vol*aP*rP.Z;
//            GMC.Elem[PCell*3+1][2]+= -Area.Y/Vol*aP*rP.Z;
//            GMC.Elem[PCell*3+2][2]+= -Area.Z/Vol*aP*rP.Z;
//
//            if (Ngb>=0)
//            {
//                int i,MatchedFaceNum=0;
//                for(i=0;i<Grids->Cells[Ngb].NoOfFace;i++)
//                    if (Mag(Sum(Grids->Cells[PCell].Area[s],
//                                Grids->Cells[Ngb  ].Area[i]))<geps)
//                    {
//                        MatchedFaceNum=i;
//                        break;
//                    }
//                rN = Minus(Grids->Cells[Ngb].ProjectedCenter[MatchedFaceNum],Grids->Cells[Ngb].CellCenter);
//
//                GMC.Elem[PCell*3  ][(s+1)*3  ]+= -Area.X/Vol*aN*rN.X;
//                GMC.Elem[PCell*3+1][(s+1)*3  ]+= -Area.Y/Vol*aN*rN.X;
//                GMC.Elem[PCell*3+2][(s+1)*3  ]+= -Area.Z/Vol*aN*rN.X;
//                GMC.Elem[PCell*3  ][(s+1)*3+1]+= -Area.X/Vol*aN*rN.Y;
//                GMC.Elem[PCell*3+1][(s+1)*3+1]+= -Area.Y/Vol*aN*rN.Y;
//                GMC.Elem[PCell*3+2][(s+1)*3+1]+= -Area.Z/Vol*aN*rN.Y;
//                GMC.Elem[PCell*3  ][(s+1)*3+2]+= -Area.X/Vol*aN*rN.Z;
//                GMC.Elem[PCell*3+1][(s+1)*3+2]+= -Area.Y/Vol*aN*rN.Z;
//                GMC.Elem[PCell*3+2][(s+1)*3+2]+= -Area.Z/Vol*aN*rN.Z;
//
//                MyST[PCell*3  ]+=Area.X/Vol*(aP*FeeP+aN*FeeNgb);
//                MyST[PCell*3+1]+=Area.Y/Vol*(aP*FeeP+aN*FeeNgb);
//                MyST[PCell*3+2]+=Area.Z/Vol*(aP*FeeP+aN*FeeNgb);
//            }
//            else
//            {
//                struct MyVector rF;
//                rF = Minus(Grids->Cells[PCell].FaceCenter[s],Grids->Cells[PCell].ProjectedCenter[s]);
//                double Dens;
//                Dens=DensCal(Grids->Variables.Alpha[PCell]);
//                MyST[PCell*3  ]+=Area.X/Vol*(FeeP+Dens*Dot(Gravity,rF));
//                MyST[PCell*3+1]+=Area.Y/Vol*(FeeP+Dens*Dot(Gravity,rF));
//                MyST[PCell*3+2]+=Area.Z/Vol*(FeeP+Dens*Dot(Gravity,rF));
//            }
//        }
//    }
//
//    double Val;
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        int i,j,k,CellNo,BdrSum;
//        BdrSum=0;
//        int NoOfFace;
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        for (s=0;s<NoOfFace;s++)
//            if (Grids->Cells[PCell].Neighbor[s]<0)
//                BdrSum++;
//        Q_SetLen(&LaspackMC,(PCell*3  )+1,(NoOfFace+1-BdrSum)*3);
//        Q_SetLen(&LaspackMC,(PCell*3+1)+1,(NoOfFace+1-BdrSum)*3);
//        Q_SetLen(&LaspackMC,(PCell*3+2)+1,(NoOfFace+1-BdrSum)*3);
//
//        i=0;j=0;k=0;
//
//        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3  )+1,GMC.Elem[PCell*3  ][0]);
//        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3+1)+1,GMC.Elem[PCell*3  ][1]);
//        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3+2)+1,GMC.Elem[PCell*3  ][2]);
//
//        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3  )+1,GMC.Elem[PCell*3+1][0]);
//        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3+1)+1,GMC.Elem[PCell*3+1][1]);
//        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3+2)+1,GMC.Elem[PCell*3+1][2]);
//
//        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3  )+1,GMC.Elem[PCell*3+2][0]);
//        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3+1)+1,GMC.Elem[PCell*3+2][1]);
//        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3+2)+1,GMC.Elem[PCell*3+2][2]);
//
//        for(s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//            {
//                CellNo = Grids->Cells[PCell].Neighbor[s];
//                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3  )+1,GMC.Elem[PCell*3  ][(s+1)*3  ]);
//                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3+1)+1,GMC.Elem[PCell*3  ][(s+1)*3+1]);
//                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3+2)+1,GMC.Elem[PCell*3  ][(s+1)*3+2]);
//
//                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3  )+1,GMC.Elem[PCell*3+1][(s+1)*3  ]);
//                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3+1)+1,GMC.Elem[PCell*3+1][(s+1)*3+1]);
//                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3+2)+1,GMC.Elem[PCell*3+1][(s+1)*3+2]);
//
//                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3  )+1,GMC.Elem[PCell*3+2][(s+1)*3  ]);
//                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3+1)+1,GMC.Elem[PCell*3+2][(s+1)*3+1]);
//                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3+2)+1,GMC.Elem[PCell*3+2][(s+1)*3+2]);
//            }
//        free(GMC.Elem[PCell*3  ]);
//        free(GMC.Elem[PCell*3+1]);
//        free(GMC.Elem[PCell*3+2]);
//
//        V_SetCmp(&ST ,(PCell*3  )+1,MyST[PCell*3  ]);
//        V_SetCmp(&ST ,(PCell*3+1)+1,MyST[PCell*3+1]);
//        V_SetCmp(&ST ,(PCell*3+2)+1,MyST[PCell*3+2]);
//        V_SetCmp(&Fee,(PCell*3  )+1,GradComp[PCell*3  ]);
//        V_SetCmp(&Fee,(PCell*3+1)+1,GradComp[PCell*3+1]);
//        V_SetCmp(&Fee,(PCell*3+2)+1,GradComp[PCell*3+2]);
//    }
//    free(GMC.Elem);
//
//    SetRTCAccuracy(1e-7);
//    Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,1000,ILUPrecond,1.0);
//
////    int dInt;
////    dInt=GetLastNoIter();
////    printf("########################  Gradient Solver Iterations %d with Error %e\n",dInt,GetLastAccuracy());
//
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        GradComp[3*PCell  ]=V_GetCmp(&Fee,(3*PCell  )+1);
//        GradComp[3*PCell+1]=V_GetCmp(&Fee,(3*PCell+1)+1);
//        GradComp[3*PCell+2]=V_GetCmp(&Fee,(3*PCell+2)+1);
//    }
//
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        Grids->Variables.PressLSEGrad[PCell].X=GradComp[3*PCell  ];
//        Grids->Variables.PressLSEGrad[PCell].Y=GradComp[3*PCell+1];
//        Grids->Variables.PressLSEGrad[PCell].Z=GradComp[3*PCell+2];
//    }
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//    free(GradComp);
//    free(MyST);
//    Q_Destr(&LaspackMC);
//    V_Destr(&ST);
//    V_Destr(&Fee);
//}
//
//
//void FindPressureCoeff(double *a1,struct Grid *Grids,int PCell)
//{
//    int i;
//    int NoOfFace;
//    double DensF,NuF,AMag;
//    double aP,aN,ST;
//
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    for (i=0;i<NoOfFace;i++)
//        if (Grids->Cells[PCell].Neighbor[i]>=0)
//        {
//            AMag=Mag(Grids->Cells[PCell].Area[i]);
//            GetNormalGradParts_Disc_AxNode(&aP,&aN,&ST,Grids->Variables.NP,
//                                           Grids->Variables.PressLSEGrad,PCell,i,Grids);
////            GetNormalGradParts_Disc_OverRelaxed(&aP,&aN,&ST,Grids->Variables.NP,
////                                                Grids->Variables.PressLSEGrad,PCell,i,Grids);
//            a1[0]  += aP*AMag;
//            a1[i+1] = aN*AMag;
//            a1[NoOfFace+1]  += ST*AMag;
//        }
////        else
////        {
////        }
//}
//
///*
//void CompSurfFlux(struct Grid *Grids)
//{
//    struct MyVector VF;
//    double Added2F;
//    struct MyVector Area;
//    int PCell,i,NoOfFace;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        for(i=0;i<NoOfFace;i++)
//            if (Grids->Cells[PCell].Neighbor[i]>=0)
//            {
//                VF.X=FaceValue(Grids->Variables.NU,Grids,PCell,i);
//                VF.Y=FaceValue(Grids->Variables.NV,Grids,PCell,i);
//                VF.Z=FaceValue(Grids->Variables.NW,Grids,PCell,i);
//                double deltaP,W,DensP,DensNgb;
//                int NgbCell;
//                deltaP = sqrt(Dot(Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[PCell].FaceCenter[i]))
//                                 ,Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[PCell].FaceCenter[i]))));
//                W=Weight(Grids,PCell,i);
//                NgbCell=Grids->Cells[PCell].Neighbor[i];
//                DensP=DensCal(Grids->Variables.Alpha[PCell]);
//                DensNgb=DensCal(Grids->Variables.Alpha[NgbCell]);
//                W=(DensP*W)/(DensNgb*(1.0-W)+DensP*W);
//                Area=Grids->Cells[PCell].Area[i];
//                Added2F=1.0/(deltaP*DensP)*W*(Grids->Variables.NP[Grids->Cells[PCell].Neighbor[i]]-Grids->Variables.NP[PCell])*sqrt(Dot(Area,Area));
//                Grids->Variables.NSF[PCell].Face[i]=Dot(VF,Area)-dt*Added2F;
//            }
//    }// PCell for
//}
//*/
//
//
//void CompSurfFlux(struct Grid *Grids)
//{
//    struct MyVector VF,Area;
//    double Added2F,aP,aN,ST;
//    double DensF,NuF;
//    int NoOfFace,PCell,i;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        for(i=0;i<NoOfFace;i++)
//            if (Grids->Cells[PCell].Neighbor[i]>=0)
//            {
//                Area=Grids->Cells[PCell].Area[i];
//                VF.X=FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,i);
//                VF.Y=FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,i);
//                VF.Z=FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,i);
//                GetNormalGradParts_Disc_AxNode(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
////                GetNormalGradParts_Disc_OverRelaxed(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
//                Added2F = (aN*Grids->Variables.NP[Grids->Cells[PCell].Neighbor[i]]
//                          +aP*Grids->Variables.NP[PCell]
//                          +ST)*Mag(Area);
//                Grids->Variables.NSF[PCell].Face[i]=Dot(VF,Area)-dt*Added2F;
//            }// if
//    }// PCell for
//}
//
///*
//void FindDiffCoeff(double *a1,struct Grid *Grids,int PCell,int Equation)
//{
//    int i;
//    int NoOfFace;
//    double ViscF,DensF,NuF;
//    struct MyVector S,K,Delta,D,GradP,GradNgb,GradF;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    switch (Equation)
//    {
//        case 'U': GradP = Gradient_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,PCell,Grids); break;
//        case 'V': GradP = Gradient_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,PCell,Grids); break;
//        case 'W': GradP = Gradient_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,PCell,Grids); break;
//    }
//    for (i=0;i<NoOfFace;i++)
//    {
//        if (Grids->Cells[PCell].Neighbor[i]>=0)
//        {
//            switch (Equation)
//            {
//             case 'U': GradNgb = Gradient_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids->Cells[PCell].Neighbor[i],Grids); break;
//             case 'V': GradNgb = Gradient_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids->Cells[PCell].Neighbor[i],Grids); break;
//             case 'W': GradNgb = Gradient_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids->Cells[PCell].Neighbor[i],Grids); break;
//            }
//            GradF=SVP(0.5,Sum(GradP,GradNgb));
//            ViscF=0.5*(ViscCal(Grids->Variables.Alpha[PCell])+ViscCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            DensF=0.5*(DensCal(Grids->Variables.Alpha[PCell])+DensCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            NuF=ViscF/DensF;
//            D=Sum(Grids->Cells[PCell].CellCenter,MinusVec(Grids->Cells[Grids->Cells[PCell].Neighbor[i]].CellCenter));
//            S=Grids->Cells[PCell].Area[i];
////Minimum Correction Approach
////      Delta=SVP(Dot(D,S)/Dot(D,D),D);
////Orthogonal Correction Approach
////      Delta=SVP(sqrt(Dot(S,S)/Dot(D,D)),D);
////Over-relaxed Approach
//            Delta=SVP(Dot(S,S)/Dot(D,S),D);
//
//            K=Sum(S,MinusVec(Delta));
//            a1[0]  +=-NuF*sqrt(Dot(Delta,Delta)/Dot(D,D));
//            a1[i+1] = NuF*sqrt(Dot(Delta,Delta)/Dot(D,D));
//            a1[NoOfFace+1]  += NuF*Dot(K,GradF);
//        }
//    }
//}
//*/
//
//void FindDiffCoeff_AxNode(double *a1,struct Grid *Grids,int PCell,int Equation)
//{
//    int i;
//    int NoOfFace;
//    double ViscF,DensF,NuF,AMag,aP,aN,ST;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    for (i=0;i<NoOfFace;i++)
//        if (Grids->Cells[PCell].Neighbor[i]>=0)
//        {
//            switch (Equation)
//            {
//                case 'U': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NU
//                                            ,Grids->Variables.ULSEGrad,PCell,i,Grids); break;
//                case 'V': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NV
//                                            ,Grids->Variables.VLSEGrad,PCell,i,Grids); break;
//                case 'W': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NW
//                                            ,Grids->Variables.WLSEGrad,PCell,i,Grids); break;
//            }
//            AMag=Mag(Grids->Cells[PCell].Area[i]);
//            ViscF=0.5*(ViscCal(Grids->Variables.Alpha[PCell])+ViscCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            DensF=0.5*(DensCal(Grids->Variables.Alpha[PCell])+DensCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            NuF=ViscF/DensF;
//
//            a1[0]  += NuF*aP*AMag;
//            a1[i+1] = NuF*aN*AMag;
//            a1[NoOfFace+1]  += NuF*ST*AMag;
//        }
//}
//
//
//
//void FindAllDiffCoeff_AxNode(double *ad,
//                             double *ST_Diff_U,
//                             double *ST_Diff_V,
//                             double *ST_Diff_W,
//                             struct Grid *Grids,int PCell)
//{
//    int i;
//    int NoOfFace;
//    double ViscF,DensF,NuF,AMag,aP,aN,UST,VST,WST;
//
//    *ST_Diff_U=*ST_Diff_V=*ST_Diff_W=0.0;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    for (i=0;i<NoOfFace;i++)
//        if (Grids->Cells[PCell].Neighbor[i]>=0)
//        {
//            GetAllNormalGradParts(&aP,&aN,
//                                  &UST,&VST,&WST,
//                                  Grids->Variables.NU,Grids->Variables.ULSEGrad,
//                                  Grids->Variables.NV,Grids->Variables.VLSEGrad,
//                                  Grids->Variables.NW,Grids->Variables.WLSEGrad,
//                                  PCell,i,Grids);
//
//            AMag=Mag(Grids->Cells[PCell].Area[i]);
//            ViscF=0.5*(ViscCal(Grids->Variables.Alpha[PCell])+ViscCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            DensF=0.5*(DensCal(Grids->Variables.Alpha[PCell])+DensCal(Grids->Variables.Alpha[Grids->Cells[PCell].Neighbor[i]]));
//            NuF=ViscF/DensF;
//
//            ad[0]  += NuF*aP*AMag;
//            ad[i+1] = NuF*aN*AMag;
//            *ST_Diff_U += -NuF*UST*AMag;
//            *ST_Diff_V += -NuF*VST*AMag;
//            *ST_Diff_W += -NuF*WST*AMag;
//        }
//}
//
//
//void FindAllConvCoeff_AxNode(double *ac,double *ST_Conv_U,double *ST_Conv_V,double *ST_Conv_W,
//                             struct Grid *Grids,int PCell,
//                             double *U,struct MyVector *LSEGradU,
//                             double *V,struct MyVector *LSEGradV,
//                             double *W,struct MyVector *LSEGradW,char Method)
//{
//    double Flux,G=.7,aP,aN,UST,VST,WST;
//    struct MyVector Area,FaceCenter,C;
//    int s,NoOfFace;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    (*ST_Conv_U)=(*ST_Conv_V)=(*ST_Conv_W)=0.0;
//    for (s=0;s<NoOfFace;s++)
//    {
//        if (Grids->Cells[PCell].Neighbor[s]>=0)
//        {
//            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
//            Area=Grids->Cells[PCell].Area[s];
//            C=RigidBodyVel(Grids,FaceCenter);
//            Flux=Grids->Variables.NSF[PCell].Face[s]-Dot(C,Area);
////  Cell Matrix Coefficients Calculation with different interpolation schemes
//            switch (Method)
//            {
//                case 'C': // Linear Interpolation (CDS)
//                {
//                    GetAllFaceValueParts(&aP,&aN,
//                                         &UST,&VST,&WST,
//                                         U,LSEGradU,
//                                         V,LSEGradV,
//                                         W,LSEGradW,
//                                         PCell,s,Grids);
//                    ac[0]+=aP*Flux;
//                    ac[s+1]=aN*Flux;
//                    (*ST_Conv_U)+=-UST*Flux;
//                    (*ST_Conv_V)+=-VST*Flux;
//                    (*ST_Conv_W)+=-WST*Flux;
//                }break;
//                case 'T': // Interpolation (UDS+CDS)
//                {
//                    ac[0]+=G*0.5*Flux+(1.0-G)*max(Flux,0.0);
//                    ac[s+1]=G*0.5*Flux+(1.0-G)*(-max(-Flux,0.0));
//                    (*ST_Conv_U)+=0.0;
//                    (*ST_Conv_V)+=0.0;
//                    (*ST_Conv_W)+=0.0;
//                }break;
//
//                case 'U': // Upwind Interpolation (UDS)
//                {
//                    GetAllUpwindFaceValueParts(&aP,&aN,
//                                               &UST,&VST,&WST,
//                                               U,LSEGradU,
//                                               V,LSEGradV,
//                                               W,LSEGradW,
//                                               PCell,s,Grids);
//                    ac[0]+=aP*Flux;
//                    ac[s+1]=aN*Flux;
//                    (*ST_Conv_U)+=-UST*Flux;
//                    (*ST_Conv_V)+=-VST*Flux;
//                    (*ST_Conv_W)+=-WST*Flux;
//                }break;
//
//                case 'D': // Deffered CDS and UDS Interpolation
//                {
//                    GetAllFaceValueParts(&aP,&aN,
//                                         &UST,&VST,&WST,
//                                         U,LSEGradU,
//                                         V,LSEGradV,
//                                         W,LSEGradW,
//                                         PCell,s,Grids);
//                    ac[0]+=max(Flux,0.0);
//                    ac[s+1]=-max(-Flux,0.0);
//                    (*ST_Conv_U)+=(max(Flux,0.0)*U[PCell]-max(-Flux,0.0)*U[Grids->Cells[PCell].Neighbor[s]])
//                                 -(aP*U[PCell]+aN*U[Grids->Cells[PCell].Neighbor[s]]+UST)*Flux;
//                    (*ST_Conv_V)+=(max(Flux,0.0)*V[PCell]-max(-Flux,0.0)*V[Grids->Cells[PCell].Neighbor[s]])
//                                 -(aP*V[PCell]+aN*V[Grids->Cells[PCell].Neighbor[s]]+VST)*Flux;
//                    (*ST_Conv_W)+=(max(Flux,0.0)*W[PCell]-max(-Flux,0.0)*W[Grids->Cells[PCell].Neighbor[s]])
//                                 -(aP*W[PCell]+aN*W[Grids->Cells[PCell].Neighbor[s]]+WST)*Flux;
//                }break;
//            }
//        }
//    }
//}
//
//void FindConvCoeff_AxNode(double *a2,struct Grid *Grids,int PCell,double *H,struct MyVector *LSEGrad,char Method)
//{
//    double Flux,G=.7,aP,aN,ST;
//    struct MyVector Area,FaceCenter,W;
//    int s,NoOfFace;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//    a2[NoOfFace+1]=0.0;
//    for (s=0;s<NoOfFace;s++)
//    {
//        if (Grids->Cells[PCell].Neighbor[s]>=0)
//        {
//            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
//            Area=Grids->Cells[PCell].Area[s];
//            W=RigidBodyVel(Grids,FaceCenter);
//            Flux=Grids->Variables.NSF[PCell].Face[s]-Dot(W,Area);
//      //Cell Matrix Coefficients Calculation with different interpolation schemes
//            switch (Method)
//            {
//                case 'C': // Linear Interpolation (CDS)
//                {
//                    GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,s,Grids);
//                    a2[0]=a2[0]+aP*Flux;
//                    a2[s+1]=aN*Flux;
//                    a2[NoOfFace+1]+=-ST*Flux;
//                }break;
//                case 'T': // Interpolation (UDS+CDS)
//                {
//                    a2[0]=a2[0]+G*0.5*Flux+(1.0-G)*max(Flux,0.0);
//                    a2[s+1]=G*0.5*Flux+(1.0-G)*(-max(-Flux,0.0));
//                    a2[NoOfFace+1]=0.0;
//                }break;
//
//                case 'U': // Upwind Interpolation (UDS)
//                {
//                    a2[0]=a2[0]+max(Flux,0.0);
//                    a2[s+1]=-max(-Flux,0.0);
//                    a2[NoOfFace+1]=0.0;
//                }break;
//
//                case 'D': // Deffered CDS and UDS Interpolation
//                {
//                    GetFaceValueParts(&aP,&aN,&ST,LSEGrad,PCell,s,Grids);
//                    a2[0]=a2[0]+max(Flux,0.0);
//                    a2[s+1]=-max(-Flux,0.0);
//                    a2[NoOfFace+1]=+(max(Flux,0.0)*H[PCell]-max(-Flux,0.0)*H[Grids->Cells[PCell].Neighbor[s]])
//                                   -(aP*H[PCell]+aN*H[Grids->Cells[PCell].Neighbor[s]]+ST)*Flux;
//                }break;
//            }
//        }
//    }
//}
//
//
//void FindConvCoeffLagged(double *a2,struct Grid *Grids,int PCell,double *H,char Method)
//{
//    double Flux,G=.7;
//    struct MyVector Area,FaceCenter,W;
//    int s,NoOfFace;
//    NoOfFace=Grids->Cells[PCell].NoOfFace;
//
//    for (s=0;s<NoOfFace;s++)
//    {
//        if (Grids->Cells[PCell].Neighbor[s]>=0)
//        {
//            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
//            Area=Grids->Cells[PCell].Area[s];
//            W=LaggedRigidBodyVel(Grids,FaceCenter);
//            Flux=Grids->Variables.SF[PCell].Face[s]-Dot(W,Area);
//            //Cell Matrix Coefficients Calculation with different interpolation schemes
//            switch (Method)
//            {
//                case 'C': // Linear Interpolation (CDS)
//                {
//                    a2[0]=a2[0]+0.5*Flux;
//                    a2[s+1]=0.5*Flux;
//                    a2[NoOfFace+1]=0.0;
//                }break;
//
//                case 'T': // Interpolation (UDS+CDS)
//                {
//                    a2[0]=a2[0]+G*0.5*Flux+(1.0-G)*max(Flux,0.0);
//                    a2[s+1]=G*0.5*Flux+(1.0-G)*(-max(-Flux,0.0));
//                    a2[NoOfFace+1]=0.0;
//                }break;
//
//                case 'U': // Upwind Interpolation (UDS)
//                {
//                    a2[0]=a2[0]+max(Flux,0.0);
//                    a2[s+1]=-max(-Flux,0.0);
//                    a2[NoOfFace+1]=0.0;
//                }break;
//
//                case 'D': // Deffered CDS and UDS Interpolation
//                {
//                    a2[0]=a2[0]+max(Flux,0.0);
//                    a2[s+1]=-max(-Flux,0.0);
//                    a2[NoOfFace+1]=+(max(Flux,0.0)*H[PCell]-max(-Flux,0.0)*H[Grids->Cells[PCell].Neighbor[s]])
//                                   -(H[PCell]+H[Grids->Cells[PCell].Neighbor[s]])*0.5*Flux;
//                }break;
//            }
//        }
//    }
//}
//
///*
//void MakeMomentumMatrix(struct MatrixCoefficient *SMC,struct Grid *Grids,char Equation
//                       ,char ConvectionMethod
//                       ,char TemporalTermDiscretizationMethodForMomentum
//                       ,char TimeDiscretizationMethodForConvection
//                       ,char TimeDiscretizationMethodForDiffusion)
//{
//    int PCell,s;
//    double *NH,*H
//          ,Omega0=1.0
//          ,DGamma0=1.0
//          ,CGamma0=1.0;
//    double ad[8]
//          ,ac[8];
//    struct MyVector *LSEGrad;
//
//    switch (TemporalTermDiscretizationMethodForMomentum)
//    {
//        case 'A' : Omega0 =     1.0;break;
//        case 'B' : Omega0 = 3.0/2.0;break;
//    }
//
//    switch (TimeDiscretizationMethodForConvection)
//    {
//        case 'I' : CGamma0 = 0.0; break;
//        case 'C' : CGamma0 = 0.5; break;
//        case 'E' : CGamma0 = 1.0; break;
//    }
//
//    switch (TimeDiscretizationMethodForDiffusion)
//    {
//        case 'I' : DGamma0 = 0.0; break;
//        case 'C' : DGamma0 = 0.5; break;
//        case 'E' : DGamma0 = 1.0; break;
//    }
//
//    switch (Equation)
//    {
//        case 'U': {NH=Grids->Variables.NU;H=Grids->Variables.U;LSEGrad=Grids->Variables.ULSEGrad;} break;
//        case 'V': {NH=Grids->Variables.NV;H=Grids->Variables.V;LSEGrad=Grids->Variables.VLSEGrad;} break;
//        case 'W': {NH=Grids->Variables.NW;H=Grids->Variables.W;LSEGrad=Grids->Variables.WLSEGrad;} break;
//    }
//
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        int NoOfFace;
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        SetSVect(ad,8,0.0);
//        SetSVect(ac,8,0.0);
//        if (TimeDiscretizationMethodForConvection=='E')
//            FindConvCoeffLagged(ac,Grids,PCell,H,ConvectionMethod);
//        else
//            FindConvCoeff_AxNode(ac,Grids,PCell,NH,LSEGrad,ConvectionMethod);
////            FindConvCoeff(ac,Grids,PCell,NH,ConvectionMethod);
//
//        FindDiffCoeff_AxNode(ad,Grids,PCell,Equation);
//
//        SMC->Elem[PCell][0]=(Grids->Cells[PCell].Volume/dt*Omega0-ad[0]*(1.0-DGamma0)+ac[0]*(1.0-CGamma0));
//
//        for (s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//                SMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0)-ad[s+1]*(1.0-DGamma0);
//
//        SMC->ST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*H[PCell])+(ac[NoOfFace+1]-ad[NoOfFace+1]);
//        for (s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//                SMC->ST[PCell]+=-((-ad[s+1]*DGamma0+ac[s+1]*CGamma0)*H[Grids->Cells[PCell].Neighbor[s]]);
//
//    }
//}
//*/
//
//void MakeAllMomentumMatrix(struct MatrixCoefficient *SMC,
//                           double *UST,double *VST,double *WST,
//                           struct Grid *Grids,
//                           char ConvectionMethod,
//                           char TemporalTermDiscretizationMethodForMomentum,
//                           char TimeDiscretizationMethodForConvection,
//                           char TimeDiscretizationMethodForDiffusion)
//{
//    int PCell,s;
//    double *NU,*U
//          ,*NV,*V
//          ,*NW,*W
//          ,Omega0=1.0
//          ,DGamma0=1.0
//          ,CGamma0=1.0
//          ,ST_Conv_U=0.0
//          ,ST_Conv_V=0.0
//          ,ST_Conv_W=0.0
//          ,ST_Diff_U=0.0
//          ,ST_Diff_V=0.0
//          ,ST_Diff_W=0.0;
//    double ad[8]
//          ,ac[8];
//    struct MyVector *LSEGradU
//                   ,*LSEGradV
//                   ,*LSEGradW;
//
//    switch (TemporalTermDiscretizationMethodForMomentum)
//    {
//        case 'A' : Omega0 =     1.0;break;
//        case 'B' : Omega0 = 3.0/2.0;break;
//    }
//
//    switch (TimeDiscretizationMethodForConvection)
//    {
//        case 'I' : CGamma0 = 0.0; break;
//        case 'C' : CGamma0 = 0.5; break;
//        case 'E' : CGamma0 = 1.0; break;
//    }
//
//    switch (TimeDiscretizationMethodForDiffusion)
//    {
//        case 'I' : DGamma0 = 0.0; break;
//        case 'C' : DGamma0 = 0.5; break;
//        case 'E' : DGamma0 = 1.0; break;
//    }
//
//    NU=Grids->Variables.NU;U=Grids->Variables.U;LSEGradU=Grids->Variables.ULSEGrad;
//    NV=Grids->Variables.NV;V=Grids->Variables.V;LSEGradV=Grids->Variables.VLSEGrad;
//    NW=Grids->Variables.NW;W=Grids->Variables.W;LSEGradW=Grids->Variables.WLSEGrad;
//
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        int NoOfFace;
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        SetSVect(ad,8,0.0);
//        SetSVect(ac,8,0.0);
//
//        FindAllConvCoeff_AxNode(ac,&ST_Conv_U,&ST_Conv_V,&ST_Conv_W,
//                                Grids,PCell,
//                                U,LSEGradU,
//                                V,LSEGradV,
//                                W,LSEGradW,ConvectionMethod);
//        FindAllDiffCoeff_AxNode(ad,&ST_Diff_U,&ST_Diff_V,&ST_Diff_W,Grids,PCell);
//
//        SMC->Elem[PCell][0]=(Grids->Cells[PCell].Volume/dt*Omega0-ad[0]*(1.0-DGamma0)+ac[0]*(1.0-CGamma0));
//
//        for (s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//                SMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0)-ad[s+1]*(1.0-DGamma0);
//
//
//
//        UST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*U[PCell])+(ST_Conv_U-ST_Diff_U);
//        VST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*V[PCell])+(ST_Conv_V-ST_Diff_V);
//        WST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*W[PCell])+(ST_Conv_W-ST_Diff_W);
//        for (s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//            {
//                UST[PCell]+=-((-ad[s+1]*DGamma0+ac[s+1]*CGamma0)*U[Grids->Cells[PCell].Neighbor[s]]);
//                VST[PCell]+=-((-ad[s+1]*DGamma0+ac[s+1]*CGamma0)*V[Grids->Cells[PCell].Neighbor[s]]);
//                WST[PCell]+=-((-ad[s+1]*DGamma0+ac[s+1]*CGamma0)*W[Grids->Cells[PCell].Neighbor[s]]);
//            }
//
//    }
//}
//
//void SetSourceTerms(double *UST,double *VST,double *WST,
//                    struct Grid *Grids,char TemporalTermDiscretizationMethodForMomentum)
//{
//    int PCell;
//    double Vol,Dens;
//    struct MyVector GradP,Vel,Added;
//    double Omega1,Omega2;
//    switch (TemporalTermDiscretizationMethodForMomentum)
//    {
//        case 'A' : {Omega1 =     1.0 ; Omega2 =     0.0;}break;
//        case 'B' : {Omega1 = 4.0/2.0 ; Omega2 = 1.0/2.0;}break;
//    }
//
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
//        Dens=DensCal(Grids->Variables.Alpha[PCell]);
//        GradP=Gradient_Disc_AxNode(Grids->Variables.P,Grids->Variables.PressLSEGrad,PCell,Grids);
//        Added=SVP(dt/Dens,GradP);
//        Dens=SelectDens(Vel,Added,Dens);
//
//        Vol=Grids->Cells[PCell].Volume;
//        UST[PCell]+=(Grids->Variables.U[PCell]*Omega1-Grids->Variables.PU[PCell]*Omega2)*Vol/dt+Gravity.X*Vol-GradP.X*Vol/Dens;
//        VST[PCell]+=(Grids->Variables.V[PCell]*Omega1-Grids->Variables.PV[PCell]*Omega2)*Vol/dt+Gravity.Y*Vol-GradP.Y*Vol/Dens;
//        WST[PCell]+=(Grids->Variables.W[PCell]*Omega1-Grids->Variables.PW[PCell]*Omega2)*Vol/dt+Gravity.Z*Vol-GradP.Z*Vol/Dens;
//    }
//}
//
//void SetOutFlowFlux(struct Grid *Grids)
//{
//    int PCell,k,Type,FaceNum;
//    double InflowFlux=0.0,OutflowFlux=0.0,Landa=1.0;
//    for(k=0;k<Grids->NB;k++)
//    {
//        PCell   = Grids->BoundaryCells[k].CellNum;
//        Type    = Grids->BoundaryCells[k].Type;
//        FaceNum = Grids->BoundaryCells[k].FaceNum;
//        switch (Type)
//        {
//            case 3:// Dirichlet Flow Function;
//                OutflowFlux += Grids->Variables.NSF[PCell].Face[FaceNum];
//            break;
//            case 5:// Dirichlet Flow Function;
//                InflowFlux += Grids->Variables.NSF[PCell].Face[FaceNum];
//            break;
//        }
//    }
//    Landa=-InflowFlux/(OutflowFlux+1e-6);
//    for(k=0;k<Grids->NB;k++)
//    {
//        PCell   = Grids->BoundaryCells[k].CellNum;
//        Type    = Grids->BoundaryCells[k].Type;
//        FaceNum = Grids->BoundaryCells[k].FaceNum;
//        switch (Type)
//        {
//            case 3:// Dirichlet Flow Function;
//                Grids->Variables.SF[PCell].Face[FaceNum]=(Grids->Variables.NSF[PCell].Face[FaceNum]*=Landa);
//            break;
//        }
//    }
//}
//
//void SetSurfFluxBoundary(struct Grid *Grids)
//{
//    int PCell,k,Type,FaceNum;
//    struct MyVector Area,FaceCenter,V,W;
//    double Value;
//    for(k=0;k<Grids->NB;k++)
//    {
//        PCell   = Grids->BoundaryCells[k].CellNum;
//        Type    = Grids->BoundaryCells[k].Type;
//        FaceNum = Grids->BoundaryCells[k].FaceNum;
//        Value   = Grids->BoundaryCells[k].Value;
//        switch (Type)
//        {
//            case 1:// Rigid Body Wall;
//            {
//                FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
//                Area=Grids->Cells[PCell].Area[FaceNum];
//                W=RigidBodyVel(Grids,FaceCenter);
//                Grids->Variables.NSF[PCell].Face[FaceNum]=Dot(W,Area);
//            }break;
//            case 3:// Outflow;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                V.X=FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,FaceNum);
//                V.Y=FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,FaceNum);
//                V.Z=FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,FaceNum);
//                Grids->Variables.NSF[PCell].Face[FaceNum]=Dot(Area,V);
//            }break;
//            case 4:// Symmetry;
//            {
//                Grids->Variables.NSF[PCell].Face[FaceNum]=0.0;;
//            }break;
//            case 5:// Dirichlet Flow Function;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                V=FlowFunction(FaceCenter,time);
//                Grids->Variables.NSF[PCell].Face[FaceNum]=Dot(Area,V);
//            }break;
//            case 6:// Static Wall;
//            {
//                Grids->Variables.NSF[PCell].Face[FaceNum]=0.0;;
//            }break;
//        }
//    }
//}
//
//void SetMomentumBC(struct MatrixCoefficient *MC
//                  ,double *UST
//                  ,double *VST
//                  ,double *WST
//                  ,struct Grid *Grids)
//{
//    int PCell,k,Type,FaceNum;
//    struct MyVector n,Area,FaceCenter,V,W,rC;
//    double Amag,Flux,Value,dn,nu;
//    for(k=0;k<Grids->NB;k++)
//    {
//        PCell   = Grids->BoundaryCells[k].CellNum;
//        Type    = Grids->BoundaryCells[k].Type;
//        FaceNum = Grids->BoundaryCells[k].FaceNum;
//        Value   = Grids->BoundaryCells[k].Value;
//        switch (Type)
//        {
//            case 1:// Rigid Body Wall Flow Function;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                Amag=Mag(Area);
//                n=Normalized(Area);
//                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
//                nu=ViscCal(Grids->Variables.Alpha[PCell])/DensCal(Grids->Variables.Alpha[PCell]);
//                W=RigidBodyVel(Grids,FaceCenter);
//                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
//
//            //  Set Diffusive Flux Boundary
//                MC->Elem[PCell][0]+=nu/dn*Amag;
//
//                UST[PCell]+=nu/dn*Amag*(W.X-Dot(Grids->Variables.ULSEGrad[PCell],rC));
//                VST[PCell]+=nu/dn*Amag*(W.Y-Dot(Grids->Variables.VLSEGrad[PCell],rC));
//                WST[PCell]+=nu/dn*Amag*(W.Z-Dot(Grids->Variables.WLSEGrad[PCell],rC));
//            }break;
//            case 3:// Out-Flow;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                Amag=sqrt(Dot(Area,Area));
//                n=SVP((1./Amag),Area);
//                V.X=Grids->Variables.U[PCell];
//                V.Y=Grids->Variables.V[PCell];
//                V.Z=Grids->Variables.W[PCell];
//                dn=fabs((Dot(Sum(Grids->Cells[PCell].CellCenter,MinusVec(FaceCenter)),n)));
//                nu=ViscCal(Grids->Variables.Alpha[PCell])/DensCal(Grids->Variables.Alpha[PCell]);
//                Flux=Grids->Variables.NSF[PCell].Face[FaceNum];
//                //  Set Diffusive Flux Boundary
//                MC->Elem[PCell][0]+=nu/dn*Amag;
//
//                UST[PCell]+=nu/dn*V.X*Amag;
//                VST[PCell]+=nu/dn*V.Y*Amag;
//                WST[PCell]+=nu/dn*V.Z*Amag;
//                //  Set Convective Flux Boundary
//                UST[PCell]+=-Flux*V.X;
//                VST[PCell]+=-Flux*V.Y;
//                WST[PCell]+=-Flux*V.Z;
//            }break;
//            case 4:// Symmetry;
//            {
//            }break;
//            case 5:// Dirichlet Flow Function;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                Amag=Mag(Area);
//                n=Normalized(Area);
//                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
//                nu=ViscCal(Grids->Variables.Alpha[PCell])/DensCal(Grids->Variables.Alpha[PCell]);
//                W=RigidBodyVel(Grids,FaceCenter);
//                V=FlowFunction(FaceCenter,time);
//                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
//
//                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);
//                //  Set Diffusive Flux Boundary
//                MC->Elem[PCell][0]+=nu/dn*Amag;
//
//                UST[PCell]+=nu/dn*Amag*(V.X-Dot(Grids->Variables.ULSEGrad[PCell],rC));
//                VST[PCell]+=nu/dn*Amag*(V.Y-Dot(Grids->Variables.VLSEGrad[PCell],rC));
//                WST[PCell]+=nu/dn*Amag*(V.Z-Dot(Grids->Variables.WLSEGrad[PCell],rC));
//                //  Set Convective Flux Boundary
//                UST[PCell]+=-Flux*V.X;
//                VST[PCell]+=-Flux*V.Y;
//                WST[PCell]+=-Flux*V.Z;
//            }break;
//            case 6:// Static Wall;
//            {
//                Area =  Grids->Cells[PCell].Area[FaceNum];
//                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
//                Amag=sqrt(Dot(Area,Area));
//                n=SVP((1./Amag),Area);
//                dn=fabs((Dot(Sum(Grids->Cells[PCell].CellCenter,MinusVec(FaceCenter)),n)));
//                nu=ViscCal(Grids->Variables.Alpha[PCell])/DensCal(Grids->Variables.Alpha[PCell]);
//            //  Set Diffusive Flux Boundary
//                MC->Elem[PCell][0]+=nu/dn*Amag;
//
//                UST[PCell]+=0.0;
//                VST[PCell]+=0.0;
//                WST[PCell]+=0.0;
//            }break;
//        }
//    }
//}
//
//void MallocMatrixes(struct Grid *Grids
//                   ,QMatrix *LaspackMC
//                   ,struct MatrixCoefficient *MC
//                   ,double **UST
//                   ,double **VST
//                   ,double **WST
//                   ,double **PST
//                   ,Vector *ST
//                   ,Vector *Fee)
//{
//    int NC=0,NCol=0,PCell;
//    NC=Grids->NC;
//
//    (*MC).Elem =(double**) malloc(sizeof(double*)*NC);
//    for(PCell=0;PCell<NC;PCell++)
//    {
//        NCol=Grids->Cells[PCell].NoOfFace+1;
//        (*MC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
//    }
//
//    *UST=(double*) malloc(sizeof(double)*NC);
//    *VST=(double*) malloc(sizeof(double)*NC);
//    *WST=(double*) malloc(sizeof(double)*NC);
//    *PST=(double*) malloc(sizeof(double)*NC);
//
//    Q_Constr(LaspackMC,"LaspackMC",NC,False,Rowws,Normal,True);
//
//    V_Constr(ST,"ST",NC,Normal,True);
//
//    V_Constr(Fee,"Fee",NC,Normal,True);
//}
//
//
//void ResetLaspackMatrix(QMatrix *LaspackMC,struct Grid *Grids)
//{
//    int NC;
//    NC=Grids->NC;
//    Q_Destr(LaspackMC);
//    Q_Constr(LaspackMC,"GMC",NC,False,Rowws,Normal,True);
//}
//
//void FreeMatrixes(struct Grid *Grids
//                 ,QMatrix *LaspackMC
//                 ,struct MatrixCoefficient *MC
//                 ,double **UST
//                 ,double **VST
//                 ,double **WST
//                 ,double **PST
//                 ,Vector *ST
//                 ,Vector *Fee)
//{
//    int NC=0,PCell;
//    NC=Grids->NC;
//    for(PCell=0;PCell<NC;PCell++)
//        free((*MC).Elem[PCell]);
//    free((*MC).Elem);
//
//    free(*UST);
//    free(*VST);
//    free(*WST);
//    free(*PST);
//
//    Q_Destr(LaspackMC);
//
//    V_Destr(ST);
//
//    V_Destr(Fee);
//}
//
//void VelsSurfIntegral(double *VelSI,struct Grid *Grids)
//{
//    struct MyVector Hf,Area;
//    int i,PCell;
//    int NoOfFace;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        VelSI[PCell]=0.0;
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        for (i=0;i<NoOfFace;i++)
//        {
//            Hf.X=FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,i);
//            Hf.Y=FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,i);
//            Hf.Z=FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,i);
//            Area=Grids->Cells[PCell].Area[i];
//            VelSI[PCell] += ((Grids->Cells[PCell].Neighbor[i]>=0)?(Dot(Hf,Area)):(Grids->Variables.NSF[PCell].Face[i]));
//        }
//    }
//}
//
//void MakePressureMatrix(struct MatrixCoefficient *SMC,double *PST,struct Grid *Grids,double *VelSI,double alfa)
//{
//    int PCell;
//    int NoOfFace,s;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        double ad[8];
//        SetSVect(ad,8,0.0);
//        FindPressureCoeff(ad,Grids,PCell);
//        SMC->Elem[PCell][0]=ad[0]/alfa;
//        for (s=0;s<NoOfFace;s++)
//            if(Grids->Cells[PCell].Neighbor[s]>=0)
//                SMC->Elem[PCell][s+1]=ad[s+1];
//        PST[PCell]=-ad[NoOfFace+1]+VelSI[PCell]/dt+((1.-alfa)/alfa)*(ad[0]*Grids->Variables.NP[PCell]);
//    }
//}
//
//void UpdateVelocities_UStar(struct Grid *Grids)
//{
//    int PCell;
//    double PDens;
//    struct MyVector GradP,Vel,Added;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
//        PDens=DensCal(Grids->Variables.Alpha[PCell]);
//        GradP=Gradient_Disc_AxNode(Grids->Variables.P,Grids->Variables.PressLSEGrad,PCell,Grids);
//        Added=SVP(dt/PDens,GradP);
//        PDens=SelectDens(Vel,Added,PDens);
//        Grids->Variables.NU[PCell] += GradP.X*dt/PDens;
//        Grids->Variables.NV[PCell] += GradP.Y*dt/PDens;
//        Grids->Variables.NW[PCell] += GradP.Z*dt/PDens;
//    }
//}
//
//void UpdateVelocities_UN(struct Grid *Grids)
//{
//    int PCell;
//    double PDens;
//    struct MyVector GradP,Added,Vel;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
//        PDens=DensCal(Grids->Variables.Alpha[PCell]);
//        GradP=Gradient_Disc_AxNode(Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,Grids);
//        Added=SVP(dt/PDens,GradP);
//        PDens=SelectDens(Vel,Added,PDens);
//        Grids->Variables.NU[PCell] += -GradP.X*dt/PDens;
//        Grids->Variables.NV[PCell] += -GradP.Y*dt/PDens;
//        Grids->Variables.NW[PCell] += -GradP.Z*dt/PDens;
//    }
//}
//
//double MassErrorCalc(struct Grid *Grids)
//{
//    double CellFlux=0.0
//          ,TotalFlux=0.0;
//    int PCell,j,NoOfFace;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        CellFlux=0.0;
//        for (j=0;j<NoOfFace;j++)
//            CellFlux+=Grids->Variables.NSF[PCell].Face[j];
//        TotalFlux+=fabs(CellFlux);
//    }
//    return TotalFlux;
//
//}
//
//void InitialConditionForFlowVariables(struct Grid *Grids)
//{
//    int PCell,FaceNum;
//    struct MyVector V,Area;
//    double A;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
////        V=FlowFunction(Grids->Cells[PCell].CellCenter,time);
////        Grids->Variables.PU[PCell] = Grids->Variables.NU[PCell] = Grids->Variables.U[PCell] = V.X;
////        Grids->Variables.PV[PCell] = Grids->Variables.NV[PCell] = Grids->Variables.V[PCell] = V.Y;
////        Grids->Variables.PW[PCell] = Grids->Variables.NW[PCell] = Grids->Variables.W[PCell] = V.Z;
//
////        Grids->Variables.PPhi[PCell] =
////        Grids->Variables.NPhi[PCell] =
////        Grids->Variables.Phi[PCell]  = Grids->Cells[PCell].CellCenter.Z-1.615;
//
//        Grids->Variables.PAlpha[PCell] =
//        Grids->Variables.NAlpha[PCell] =
//        Grids->Variables.Alpha[PCell] =1.0;// (Grids->Cells[PCell].CellCenter.Z<-1.3307)?1.0:0.0;
//
////        Grids->Variables.NP[PCell] = Grids->Variables.P[PCell]  = DensCal(Grids->Variables.Alpha[PCell])*
////                                                                  Gravity.Z*
////                                                                  (Grids->Cells[PCell].CellCenter.Z-1.615);
//
//    }
//    UpdatePressLSEGradient(Grids);
////    UpdatePhiLSEGradient(Grids);
////    UpdateVelsLSEGradient(Grids);
////    for(PCell=0;PCell<Grids->NC;PCell++)
////    {
////        int NoOfFace;
////        NoOfFace=Grids->Cells[PCell].NoOfFace;
////        for(FaceNum=0;FaceNum<NoOfFace;FaceNum++)
////        {
////            Area=Grids->Cells[PCell].Area[FaceNum];
////            V.X = FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,FaceNum);
////            V.Y = FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,FaceNum);
////            V.Z = FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,FaceNum);
////            Grids->Variables.SF[PCell].Face[FaceNum] = Grids->Variables.NSF[PCell].Face[FaceNum] = Dot(Area,V);
////            Grids->Variables.AF[PCell].Face[FaceNum] = Grids->Variables.NAF[PCell].Face[FaceNum] =
////                                                       FaceValue(Grids->Variables.NAlpha,Grids,PCell,FaceNum);
////        }
////    }
////    SetSurfFluxBoundary(Grids);
//}
//
//
//void SetDamping(struct Grid *Grids,double *UST,double *VST,double *WST,
//                                   double *u  ,double *v  ,double *w)
//{
//    double DampTermu,DampTermv,DampTermw;
//    int PCell;
//
//    for (PCell=0;PCell<Grids->NC;PCell++)
//    {
//        double R1,R2,R3,DampTerm=0.0,alpha=-50.0,Vol;
///*
//  R1=Position.X-MassCenter.X;
//  R2=-10.0;
//  R3=-30.0;
//  alpha=-50.0;
//  DampTerm=DampTerm+alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//*/
//
//        R1=Grids->Cells[PCell].CellCenter.X-Grids->DynaVariables.MassCenter.X;
//        R2=8.0;
//        R3=17.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.X-Grids->DynaVariables.MassCenter.X;
//        R2=-35.0;
//        R3=-56.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Y-Grids->DynaVariables.MassCenter.Y;
//        R2=8.0;
//        R3=19.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Y-Grids->DynaVariables.MassCenter.Y;
//        R2=-8.0;
//        R3=-19.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Z-Grids->DynaVariables.MassCenter.Z;
//        R2=3.0;
//        R3=11.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Z-Grids->DynaVariables.MassCenter.Z;
//        R2=-7.0;
//        R3=-13.0;
//        alpha=-25.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        Vol=Grids->Cells[PCell].Volume;
//        DampTermu=DampTermv=DampTermw=DampTerm=0.0;
//
//        UST[PCell]+=DampTermu*u[PCell]*Vol;
//        VST[PCell]+=DampTermv*v[PCell]*Vol;
//        WST[PCell]+=DampTermw*w[PCell]*Vol;
//
//    }
//}
//
//
//void UpdateVelsLSEGradient(struct Grid *Grids)
//{
//    int PCell;
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        Grids->Variables.ULSEGrad[PCell]=LSEGradient(Grids->Variables.NU,PCell,Grids);
//        Grids->Variables.VLSEGrad[PCell]=LSEGradient(Grids->Variables.NV,PCell,Grids);
//        Grids->Variables.WLSEGrad[PCell]=LSEGradient(Grids->Variables.NW,PCell,Grids);
//    }
//}
//
//void SolveNSEquation(struct Grid *Grids
//                    ,char ConvectionMethod
//                    ,char TemporalTermDiscretizationMethodForMomentum
//                    ,char TimeDiscretizationMethodForConvection
//                    ,char TimeDiscretizationMethodForDiffusion)
//{
//    struct MatrixCoefficient MC;
//
//    double *UST,*VST,*WST,*PST;
//    int i,j,dummyCNT,dInt;
//    int PCell;
//
//    QMatrix LaspackMC;
//    Vector ST;
//    Vector Fee;
//
//    MallocMatrixes(Grids
//                  ,&LaspackMC,&MC
//                  ,&UST,&VST,&WST,&PST,&ST,&Fee);
//
//    double *VelSI;
//    VelSI=(double *)malloc(sizeof(double)*Grids->NC);
//
////    printf("NSE started ...\n");
//
////#############################################################################
////#############################################################################
////#############################################################################
//    do
//    {
//        dInt=0;
//
//        MakeAllMomentumMatrix(&MC,
//                              UST,VST,WST,
//                              Grids,
//                              ConvectionMethod,
//                              TemporalTermDiscretizationMethodForMomentum,
//                              TimeDiscretizationMethodForConvection,
//                              TimeDiscretizationMethodForDiffusion);
//
//        SetSourceTerms(UST,VST,WST,Grids,TemporalTermDiscretizationMethodForMomentum);
//        SetMomentumBC(&MC,UST,VST,WST,Grids);
//        SetDamping(Grids,UST,VST,WST,Grids->Variables.U,Grids->Variables.V,Grids->Variables.W);
//
//        SetRTCAccuracy(1e-6);
//        ResetLaspackMatrix(&LaspackMC,Grids);
//        CopyToLaspackMatrix(&LaspackMC,&MC,Grids);
//
//        CopyToLaspackVector(&ST ,UST,Grids);
//        CopyToLaspackVector(&Fee,Grids->Variables.NU,Grids);
//        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
//        dInt+=GetLastNoIter();
//        CopyFromLaspackVector(Grids->Variables.NU,&Fee,Grids);
//
//        CopyToLaspackVector(&ST ,VST,Grids);
//        CopyToLaspackVector(&Fee,Grids->Variables.NV,Grids);
//        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
//        dInt+=GetLastNoIter();
//        CopyFromLaspackVector(Grids->Variables.NV,&Fee,Grids);
//
//        CopyToLaspackVector(&ST ,WST,Grids);
//        CopyToLaspackVector(&Fee,Grids->Variables.NW,Grids);
//        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
//        dInt+=GetLastNoIter();
//        CopyFromLaspackVector(Grids->Variables.NW,&Fee,Grids);
//
//        UpdateVelsLSEGradient(Grids);
//        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        printf("########################  Momentum Solver Iterations %d\n",dInt);
//        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//    }while(dInt>0);
//
//
////&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//    UpdateVelocities_UStar(Grids);
//
//    UpdateVelsLSEGradient(Grids);
//    VelsSurfIntegral(VelSI,Grids);
//
//    i=0;
//    do
//    {
//        i++;
//        ResetLaspackMatrix(&LaspackMC,Grids);
//        MakePressureMatrix(&MC,PST,Grids,VelSI,1.0-1e-10);
//        SetRTCAccuracy(1e-4/pow(10.0,(i<3)?i:3));
//        CopyToLaspackMatrix(&LaspackMC,&MC,Grids);
//        CopyToLaspackVector(&ST,PST,Grids);
//        CopyToLaspackVector(&Fee,Grids->Variables.NP,Grids);
////        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,1500,ILUPrecond,1.0);
//        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,1500,ILUPrecond,1.5);
//////        P=*SSORIter(&LaspackPMC,&P,&PST,500,ILUPrecond,1.5);
//        dInt=GetLastNoIter();
//        printf("####################  Pressure Solver Iterations %d with Error %e\n",dInt,GetLastAccuracy());
//        double dummy=0.0;
//        dummy = V_GetCmp(&Fee,1);
//        CopyFromLaspackVectorByOffset(Grids->Variables.NP,&Fee,dummy,Grids);
//
//
//        UpdatePressLSEGradient(Grids);
////        for(PCell=0;PCell<Grids->NC;PCell++)
////            Grids->Variables.PressLSEGrad[PCell]=Gradient(Grids->Variables.NP,PCell,Grids);
//
//
//
//    }while((dInt>0)&&(i<15));
////    printf("Continuity solved ...\n");
//
//    CompSurfFlux(Grids);
//    SetSurfFluxBoundary(Grids);
//
//    double TotalFlux=0.0,CellFlux;
//
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        CellFlux=0.0;
//        for(i=0;i<Grids->Cells[PCell].NoOfFace;i++)
//                CellFlux+=Grids->Variables.NSF[PCell].Face[i];
//        TotalFlux+=fabs(CellFlux);
//    }
//
//    printf("Total Flux = %e\n",TotalFlux);
////    SetOutFlowFlux(Grids);
//    UpdateVelocities_UN(Grids);
//    UpdateVelsLSEGradient(Grids);
//
////#############################################################################
////#############################################################################
////#############################################################################
//    free(VelSI);
//
//    FreeMatrixes(Grids
//                ,&LaspackMC
//                ,&MC,&UST,&VST,&WST,&PST
//                ,&ST
//                ,&Fee);
////    printf("NSE ended ...\n");
//}


/// mahdy
//Determination of the donner and acceptor cells of all faces for current cell
void DonnAccCell(struct DonnerAcceptor *DA,struct Grid *Grids,int PCell)
{
    double Flux;
    int FaceNum,NgbCell;
    struct MyVector Area,FaceCenter,W;
    int NoOfFace;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (FaceNum=0;FaceNum<NoOfFace;FaceNum++)
    {
        Area=Grids->Cells[PCell].Area[FaceNum];
        FaceCenter=Grids->Cells[PCell].FaceCenter[FaceNum];
        W=RigidBodyVel(Grids,FaceCenter);
        Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(W,Area);
        NgbCell=(Grids->Cells[PCell].Neighbor[FaceNum]>=0)?Grids->Cells[PCell].Neighbor[FaceNum]:PCell;
        if (fabs(Flux)<=geps)
        {
            DA->Face[FaceNum].d=PCell;
            DA->Face[FaceNum].a=NgbCell;
        }
        else if (Flux>geps)
        {
            DA->Face[FaceNum].d=PCell;
            DA->Face[FaceNum].a=NgbCell;
        }
        else
        {
            DA->Face[FaceNum].d=NgbCell;
            DA->Face[FaceNum].a=PCell;
        }
    }
}

void UpdatePressLSEGradient(struct Grid *Grids)
{
    struct MatrixCoefficient GMC;
    int PCell,NCol;
    GMC.Elem=(double**) malloc(sizeof(double*)*Grids->NC*3);
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        NCol=Grids->Cells[PCell].NoOfFace+1;

        GMC.Elem[PCell*3  ]=(double*)malloc(sizeof(double)*NCol*3);
        GMC.Elem[PCell*3+1]=(double*)malloc(sizeof(double)*NCol*3);
        GMC.Elem[PCell*3+2]=(double*)malloc(sizeof(double)*NCol*3);
    }

    QMatrix LaspackMC;
    Vector ST,Fee;
    Q_Constr(&LaspackMC,"LaspackMC",Grids->NC*3,False,Rowws,Normal,True);
    V_Constr(&ST,"ST",Grids->NC*3,Normal,True);
    V_Constr(&Fee,"Fee",Grids->NC*3,Normal,True);
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

    int s=0;
    double *GradComp,*MyST;
    GradComp = (double *)malloc(sizeof(double)*Grids->NC*3);
    MyST = (double *)malloc(sizeof(double)*Grids->NC*3);

    SetSVect(MyST,3*Grids->NC,0.0);
    SetSVect(GradComp,3*Grids->NC,0.0);
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        NCol=Grids->Cells[PCell].NoOfFace+1;
        SetSVect(GMC.Elem[PCell*3  ],3*NCol,0.0);
        SetSVect(GMC.Elem[PCell*3+1],3*NCol,0.0);
        SetSVect(GMC.Elem[PCell*3+2],3*NCol,0.0);
        GradComp[PCell*3  ]=Grids->Variables.PressLSEGrad[PCell].X;
        GradComp[PCell*3+1]=Grids->Variables.PressLSEGrad[PCell].Y;
        GradComp[PCell*3+2]=Grids->Variables.PressLSEGrad[PCell].Z;
    }

    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        GMC.Elem[PCell*3  ][0] = 1.0;
        GMC.Elem[PCell*3+1][1] = 1.0;
        GMC.Elem[PCell*3+2][2] = 1.0;

        for (s=0;s<Grids->Cells[PCell].NoOfFace;s++)
        {
            int Ngb;
            double Vol,FeeP,FeeNgb;
            struct MyVector rP,rN,Area;

            Ngb = Grids->Cells[PCell].Neighbor[s];
            Vol = Grids->Cells[PCell].Volume;
            Area= Grids->Cells[PCell].Area[s];
            double aP,aN,CST;
            GetFaceValueParts_Disc(&aP,&aN,&CST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,s,Grids);
            rP = Minus(Grids->Cells[PCell].ProjectedCenter[s],Grids->Cells[PCell].CellCenter);
            FeeP  =Grids->Variables.NP[PCell];
            FeeNgb=(Ngb>=0)?Grids->Variables.NP[Ngb]:0.0;

            GMC.Elem[PCell*3  ][0]+= -Area.X/Vol*aP*rP.X;
            GMC.Elem[PCell*3+1][0]+= -Area.Y/Vol*aP*rP.X;
            GMC.Elem[PCell*3+2][0]+= -Area.Z/Vol*aP*rP.X;
            GMC.Elem[PCell*3  ][1]+= -Area.X/Vol*aP*rP.Y;
            GMC.Elem[PCell*3+1][1]+= -Area.Y/Vol*aP*rP.Y;
            GMC.Elem[PCell*3+2][1]+= -Area.Z/Vol*aP*rP.Y;
            GMC.Elem[PCell*3  ][2]+= -Area.X/Vol*aP*rP.Z;
            GMC.Elem[PCell*3+1][2]+= -Area.Y/Vol*aP*rP.Z;
            GMC.Elem[PCell*3+2][2]+= -Area.Z/Vol*aP*rP.Z;

            if (Ngb>=0)
            {
                int i,MatchedFaceNum=0;
                for(i=0;i<Grids->Cells[Ngb].NoOfFace;i++)
                    if (Mag(Sum(Grids->Cells[PCell].Area[s],
                                Grids->Cells[Ngb  ].Area[i]))<geps)
                    {
                        MatchedFaceNum=i;
                        break;
                    }
                rN = Minus(Grids->Cells[Ngb].ProjectedCenter[MatchedFaceNum],Grids->Cells[Ngb].CellCenter);

                GMC.Elem[PCell*3  ][(s+1)*3  ]+= -Area.X/Vol*aN*rN.X;
                GMC.Elem[PCell*3+1][(s+1)*3  ]+= -Area.Y/Vol*aN*rN.X;
                GMC.Elem[PCell*3+2][(s+1)*3  ]+= -Area.Z/Vol*aN*rN.X;
                GMC.Elem[PCell*3  ][(s+1)*3+1]+= -Area.X/Vol*aN*rN.Y;
                GMC.Elem[PCell*3+1][(s+1)*3+1]+= -Area.Y/Vol*aN*rN.Y;
                GMC.Elem[PCell*3+2][(s+1)*3+1]+= -Area.Z/Vol*aN*rN.Y;
                GMC.Elem[PCell*3  ][(s+1)*3+2]+= -Area.X/Vol*aN*rN.Z;
                GMC.Elem[PCell*3+1][(s+1)*3+2]+= -Area.Y/Vol*aN*rN.Z;
                GMC.Elem[PCell*3+2][(s+1)*3+2]+= -Area.Z/Vol*aN*rN.Z;

                MyST[PCell*3  ]+=Area.X/Vol*(aP*FeeP+aN*FeeNgb);
                MyST[PCell*3+1]+=Area.Y/Vol*(aP*FeeP+aN*FeeNgb);
                MyST[PCell*3+2]+=Area.Z/Vol*(aP*FeeP+aN*FeeNgb);
            }
            else
            {
                struct MyVector rF;
                rF = Minus(Grids->Cells[PCell].FaceCenter[s],Grids->Cells[PCell].ProjectedCenter[s]);
                double Dens;
                Dens=DensCal(Grids->Variables.NAlpha[PCell]);
                MyST[PCell*3  ]+=Area.X/Vol*(FeeP+Dens*Dot(Gravity,rF));
                MyST[PCell*3+1]+=Area.Y/Vol*(FeeP+Dens*Dot(Gravity,rF));
                MyST[PCell*3+2]+=Area.Z/Vol*(FeeP+Dens*Dot(Gravity,rF));
            }
        }
    }

    double Val;
    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        int i,j,k,CellNo,BdrSum;
        BdrSum=0;
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for (s=0;s<NoOfFace;s++)
            if (Grids->Cells[PCell].Neighbor[s]<0)
                BdrSum++;
        Q_SetLen(&LaspackMC,(PCell*3  )+1,(NoOfFace+1-BdrSum)*3);
        Q_SetLen(&LaspackMC,(PCell*3+1)+1,(NoOfFace+1-BdrSum)*3);
        Q_SetLen(&LaspackMC,(PCell*3+2)+1,(NoOfFace+1-BdrSum)*3);

        i=0;j=0;k=0;

        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3  )+1,GMC.Elem[PCell*3  ][0]);
        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3+1)+1,GMC.Elem[PCell*3  ][1]);
        Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(PCell*3+2)+1,GMC.Elem[PCell*3  ][2]);

        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3  )+1,GMC.Elem[PCell*3+1][0]);
        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3+1)+1,GMC.Elem[PCell*3+1][1]);
        Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(PCell*3+2)+1,GMC.Elem[PCell*3+1][2]);

        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3  )+1,GMC.Elem[PCell*3+2][0]);
        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3+1)+1,GMC.Elem[PCell*3+2][1]);
        Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(PCell*3+2)+1,GMC.Elem[PCell*3+2][2]);

        for(s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
            {
                CellNo = Grids->Cells[PCell].Neighbor[s];
                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3  )+1,GMC.Elem[PCell*3  ][(s+1)*3  ]);
                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3+1)+1,GMC.Elem[PCell*3  ][(s+1)*3+1]);
                Q_SetEntry(&LaspackMC,(PCell*3  )+1,i++,(CellNo*3+2)+1,GMC.Elem[PCell*3  ][(s+1)*3+2]);

                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3  )+1,GMC.Elem[PCell*3+1][(s+1)*3  ]);
                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3+1)+1,GMC.Elem[PCell*3+1][(s+1)*3+1]);
                Q_SetEntry(&LaspackMC,(PCell*3+1)+1,j++,(CellNo*3+2)+1,GMC.Elem[PCell*3+1][(s+1)*3+2]);

                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3  )+1,GMC.Elem[PCell*3+2][(s+1)*3  ]);
                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3+1)+1,GMC.Elem[PCell*3+2][(s+1)*3+1]);
                Q_SetEntry(&LaspackMC,(PCell*3+2)+1,k++,(CellNo*3+2)+1,GMC.Elem[PCell*3+2][(s+1)*3+2]);
            }
        free(GMC.Elem[PCell*3  ]);
        free(GMC.Elem[PCell*3+1]);
        free(GMC.Elem[PCell*3+2]);

        V_SetCmp(&ST ,(PCell*3  )+1,MyST[PCell*3  ]);
        V_SetCmp(&ST ,(PCell*3+1)+1,MyST[PCell*3+1]);
        V_SetCmp(&ST ,(PCell*3+2)+1,MyST[PCell*3+2]);
        V_SetCmp(&Fee,(PCell*3  )+1,GradComp[PCell*3  ]);
        V_SetCmp(&Fee,(PCell*3+1)+1,GradComp[PCell*3+1]);
        V_SetCmp(&Fee,(PCell*3+2)+1,GradComp[PCell*3+2]);
    }
    free(GMC.Elem);

    SetRTCAccuracy(1e-7);
    Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,1000,ILUPrecond,1.0);

//    int dInt;
//    dInt=GetLastNoIter();
//    printf("########################  Gradient Solver Iterations %d with Error %e\n",dInt,GetLastAccuracy());

    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        GradComp[3*PCell  ]=V_GetCmp(&Fee,(3*PCell  )+1);
        GradComp[3*PCell+1]=V_GetCmp(&Fee,(3*PCell+1)+1);
        GradComp[3*PCell+2]=V_GetCmp(&Fee,(3*PCell+2)+1);
    }

    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.PressLSEGrad[PCell].X=GradComp[3*PCell  ];
        Grids->Variables.PressLSEGrad[PCell].Y=GradComp[3*PCell+1];
        Grids->Variables.PressLSEGrad[PCell].Z=GradComp[3*PCell+2];
    }
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
    free(GradComp);
    free(MyST);
    Q_Destr(&LaspackMC);
    V_Destr(&ST);
    V_Destr(&Fee);
}

void UpdateVelsLSEGradient(struct Grid *Grids)
{
    int PCell;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Grids->Variables.ULSEGrad[PCell]=LSEGradient(Grids->Variables.NU,PCell,Grids);
        Grids->Variables.VLSEGrad[PCell]=LSEGradient(Grids->Variables.NV,PCell,Grids);
        Grids->Variables.WLSEGrad[PCell]=LSEGradient(Grids->Variables.NW,PCell,Grids);
    }
}

void FindPressureCoeff(double *a1,struct Grid *Grids,int PCell)
{
    int i;
    int NoOfFace;
    double DensF,NuF,AMag;
    double aP,aN,ST;

    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
        if (Grids->Cells[PCell].Neighbor[i]>=0)
        {
            AMag=Mag(Grids->Cells[PCell].Area[i]);
            GetNormalGradParts_Disc_AxNode(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
//          GetNormalGradParts_Disc_OverRelaxed(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
            a1[0]  += aP*AMag;
            a1[i+1] = aN*AMag;
            a1[NoOfFace+1]  += ST*AMag;
        }

}

void CompSurfFlux(struct Grid *Grids)
{
    struct MyVector VF,Area;
    double Added2F,aP,aN,ST;
    double DensF,NuF;
    int NoOfFace,PCell,i;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for(i=0;i<NoOfFace;i++)
            if (Grids->Cells[PCell].Neighbor[i]>=0)
            {
                Area=Grids->Cells[PCell].Area[i];
                VF.X=FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,i);
                VF.Y=FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,i);
                VF.Z=FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,i);
                GetNormalGradParts_Disc_AxNode(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
//                GetNormalGradParts_Disc_OverRelaxed(&aP,&aN,&ST,Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,i,Grids);
                Added2F = (aN*Grids->Variables.NP[Grids->Cells[PCell].Neighbor[i]]
                          +aP*Grids->Variables.NP[PCell]
                          +ST)*Mag(Area);
                Grids->Variables.NSF[PCell].Face[i]=Dot(VF,Area)-dt*Added2F;
            }// if
    }// PCell for
}

void FindDiffCoeff_AxNode(double *a1,struct Grid *Grids,int PCell,int Equation)
{
    int i;
    int NoOfFace;
    double ViscF,DensF,NuF,AMag,aP,aN,ST,Nut;
    double w1,w2;
    NoOfFace=Grids->Cells[PCell].NoOfFace;
    for (i=0;i<NoOfFace;i++)
        if (Grids->Cells[PCell].Neighbor[i]>=0)
        {
            switch (Equation)
            {
                case 'U': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NU
                                            ,Grids->Variables.ULSEGrad,PCell,i,Grids); break;
                case 'V': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NV
                                            ,Grids->Variables.VLSEGrad,PCell,i,Grids); break;
                case 'W': GetNormalGradParts(&aP,&aN,&ST,Grids->Variables.NW
                                            ,Grids->Variables.WLSEGrad,PCell,i,Grids); break;
            }
            AMag=Mag(Grids->Cells[PCell].Area[i]);
            ViscF=0.5*(ViscCal(Grids->Variables.PAlpha[PCell])+ViscCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));
            DensF=0.5*(DensCal(Grids->Variables.PAlpha[PCell])+DensCal(Grids->Variables.PAlpha[Grids->Cells[PCell].Neighbor[i]]));

            ///Test
            w1=Mag(Minus(Grids->Cells[PCell].CellCenter
                 ,Grids->Cells[PCell].FaceCenter[i]));
            w2=Mag(Minus(Grids->Cells[Grids->Cells[PCell].Neighbor[i]].CellCenter
                 ,Grids->Cells[PCell].FaceCenter[i]));
            Nut=(w2*(Grids->Variables.NuTurb[PCell])
                     +w1*(Grids->Variables.NuTurb[Grids->Cells[PCell].Neighbor[i]]))/(w1+w2);

            ///saman
//            Nut=0.5*((Grids->Variables.NuTurb[PCell])+(Grids->Variables.NuTurb[Grids->Cells[PCell].Neighbor[i]]));
            NuF=(ViscF/DensF)+Nut;

            a1[0]  += NuF*aP*AMag;
            a1[i+1] = NuF*aN*AMag;
            a1[NoOfFace+1]-=NuF*ST*AMag;
        }
}

void FindConvCoeff_AxNode(double *a2,struct Grid *Grids,int PCell,double *H,struct MyVector *LSEGrad,char Method)
{
    double Flux,G=0.7,aP,aN,ST,kapa=0.5,betaf,AU,AN;
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

void FindConvCoeffLagged(double *a2,struct Grid *Grids,int PCell,double *H,char Method)
{
    double Flux,G=.7;
    struct MyVector Area,FaceCenter,W;
    int s,NoOfFace;
    NoOfFace=Grids->Cells[PCell].NoOfFace;

    for (s=0;s<NoOfFace;s++)
    {
        if (Grids->Cells[PCell].Neighbor[s]>=0)
        {
            FaceCenter=Grids->Cells[PCell].FaceCenter[s];
            Area=Grids->Cells[PCell].Area[s];
            W=LaggedRigidBodyVel(Grids,FaceCenter);
            Flux=Grids->Variables.SF[PCell].Face[s]-Dot(W,Area);
            //Cell Matrix Coefficients Calculation with different interpolation schemes
            switch (Method)
            {
                case 'C': // Linear Interpolation (CDS)
                {
                    a2[0]=a2[0]+0.5*Flux;
                    a2[s+1]=0.5*Flux;
                    a2[NoOfFace+1]=0.0;
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
                    a2[0]=a2[0]+max(Flux,0.0);
                    a2[s+1]=-max(-Flux,0.0);
                    a2[NoOfFace+1]=+(max(Flux,0.0)*H[PCell]-max(-Flux,0.0)*H[Grids->Cells[PCell].Neighbor[s]])
                                   -(H[PCell]+H[Grids->Cells[PCell].Neighbor[s]])*0.5*Flux;
                }break;
            }
        }
    }
}

void MakeMomentumMatrix(struct MatrixCoefficient *SMC,double *ST,struct Grid *Grids,char Equation
                       ,char ConvectionMethod
                       ,char TemporalTermDiscretizationMethodForMomentum
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

    switch (TemporalTermDiscretizationMethodForMomentum)
    {
        case 'A' : Omega0 =     1.0;break;
        case 'B' : Omega0 = 3.0/2.0;break;
    }

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

    switch (Equation)
    {
        case 'U': {NH=Grids->Variables.NU;H=Grids->Variables.U;LSEGrad=Grids->Variables.ULSEGrad;} break;
        case 'V': {NH=Grids->Variables.NV;H=Grids->Variables.V;LSEGrad=Grids->Variables.VLSEGrad;} break;
        case 'W': {NH=Grids->Variables.NW;H=Grids->Variables.W;LSEGrad=Grids->Variables.WLSEGrad;} break;
    }

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        int NoOfFace;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        SetSVect(ad,8,0.0);
        SetSVect(ac,8,0.0);

        if (TimeDiscretizationMethodForConvection=='E')
            FindConvCoeffLagged(ac,Grids,PCell,H,ConvectionMethod);
        else
            FindConvCoeff_AxNode(ac,Grids,PCell,NH,LSEGrad,ConvectionMethod);

        FindDiffCoeff_AxNode(ad,Grids,PCell,Equation);


        SMC->Elem[PCell][0]=(Grids->Cells[PCell].Volume/dt*Omega0-ad[0]*(1.0-DGamma0)+ac[0]*(1.0-CGamma0));
        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                SMC->Elem[PCell][s+1]=ac[s+1]*(1.0-CGamma0)-ad[s+1]*(1.0-DGamma0);

        ST[PCell]=-((ac[0]*CGamma0-ad[0]*DGamma0)*H[PCell])+(ac[NoOfFace+1]-ad[NoOfFace+1]);

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                ST[PCell]+= -((ac[s+1]*CGamma0-ad[s+1]*DGamma0)*H[Grids->Cells[PCell].Neighbor[s]]);

    }
}

void SetSourceTerms(double *UST,double *VST,double *WST,
                    struct Grid *Grids,char TemporalTermDiscretizationMethodForMomentum)
{
    int PCell;
    double Vol,Dens;
    struct MyVector GradP,Vel,Added,GradK;
    double Omega1,Omega2;
    switch (TemporalTermDiscretizationMethodForMomentum)
    {
        case 'A' : {Omega1 =     1.0 ; Omega2 =     0.0;}break;
        case 'B' : {Omega1 = 4.0/2.0 ; Omega2 = 1.0/2.0;}break;
    }

    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
        Dens=DensCal(Grids->Variables.PAlpha[PCell]);
        GradP=Gradient_Disc_AxNode(Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,Grids);///OK

        ///saman
        GradK=Gradient_Disc_AxNode(Grids->Variables.NK,Grids->Variables.KLSEGrad,PCell,Grids);
//        GradP=Gradient(Grids->Variables.NP,PCell,Grids);///OK

//        Added=SVP(dt/Dens,GradP);
//        Dens=SelectDens(Vel,Added,Dens);      ///   Only In The Case of Two Phase Flow...

        Vol=Grids->Cells[PCell].Volume;
        UST[PCell]+=(Grids->Variables.U[PCell]*Omega1-Grids->Variables.PU[PCell]*Omega2)*Vol/dt+Gravity.X*Vol-GradP.X*Vol/Dens-(2.0/3.0)*GradK.X*Vol;
        VST[PCell]+=(Grids->Variables.V[PCell]*Omega1-Grids->Variables.PV[PCell]*Omega2)*Vol/dt+Gravity.Y*Vol-GradP.Y*Vol/Dens-(2.0/3.0)*GradK.Y*Vol;
        WST[PCell]+=(Grids->Variables.W[PCell]*Omega1-Grids->Variables.PW[PCell]*Omega2)*Vol/dt+Gravity.Z*Vol-GradP.Z*Vol/Dens-(2.0/3.0)*GradK.Z*Vol;
    }
}

void SetOutFlowFlux(struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    double InflowFlux=0.0,OutflowFlux=0.0,Landa=1.0;
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        switch (Type)
        {
            case 3:// Outlet;
                OutflowFlux += Grids->Variables.NSF[PCell].Face[FaceNum];
            break;
            case 5:// Inlet;
                InflowFlux += Grids->Variables.NSF[PCell].Face[FaceNum];
            break;
        }
    }
    Landa=-InflowFlux/(OutflowFlux+1e-6);
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        switch (Type)
        {
            case 3:// Dirichlet Flow Function;
//                Grids->Variables.SF[PCell].Face[FaceNum]=(Grids->Variables.NSF[PCell].Face[FaceNum]*=Landa);
                Grids->Variables.NSF[PCell].Face[FaceNum]=(Grids->Variables.NSF[PCell].Face[FaceNum])*Landa;
            break;
        }
    }
}

void SetSurfFluxBoundary(struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector Area,FaceCenter,V,W,PC,PCF;
    double Value,dn,Densf,Amag;
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        Value   = Grids->BoundaryCells[k].Value;
        switch (Type)
        {
            case 1:// Rigid Body Moving Wall;
            {
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Area = Grids->Cells[PCell].Area[FaceNum];
                W = RigidBodyVel(Grids,FaceCenter);
                Grids->Variables.NSF[PCell].Face[FaceNum] = Dot(W,Area);
            }break;

            case 3:// Outflow;
            {
                if(Grids->RCells[PCell].Type==4)   ///Hex
                          {
                            if(FaceNum==0)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[1]+
                                                                                Grids->Variables.NSF[PCell].Face[2]+
                                                                                Grids->Variables.NSF[PCell].Face[3]+
                                                                                Grids->Variables.NSF[PCell].Face[4]+
                                                                                Grids->Variables.NSF[PCell].Face[5]);
                            if(FaceNum==1)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[0]+
                                                                                Grids->Variables.NSF[PCell].Face[2]+
                                                                                Grids->Variables.NSF[PCell].Face[3]+
                                                                                Grids->Variables.NSF[PCell].Face[4]+
                                                                                Grids->Variables.NSF[PCell].Face[5]);
                            if(FaceNum==2)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[1]+
                                                                                Grids->Variables.NSF[PCell].Face[0]+
                                                                                Grids->Variables.NSF[PCell].Face[3]+
                                                                                Grids->Variables.NSF[PCell].Face[4]+
                                                                                Grids->Variables.NSF[PCell].Face[5]);
                            if(FaceNum==3)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[1]+
                                                                                Grids->Variables.NSF[PCell].Face[2]+
                                                                                Grids->Variables.NSF[PCell].Face[0]+
                                                                                Grids->Variables.NSF[PCell].Face[4]+
                                                                                Grids->Variables.NSF[PCell].Face[5]);
                            if(FaceNum==4)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[1]+
                                                                                Grids->Variables.NSF[PCell].Face[2]+
                                                                                Grids->Variables.NSF[PCell].Face[3]+
                                                                                Grids->Variables.NSF[PCell].Face[0]+
                                                                                Grids->Variables.NSF[PCell].Face[5]);
                            if(FaceNum==5)
                            Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*  (Grids->Variables.NSF[PCell].Face[1]+
                                                                                Grids->Variables.NSF[PCell].Face[2]+
                                                                                Grids->Variables.NSF[PCell].Face[3]+
                                                                                Grids->Variables.NSF[PCell].Face[4]+
                                                                                Grids->Variables.NSF[PCell].Face[0]);
                          }
//
//                          if(Grids->RCells[PCell].Type==5)   ///Prism
//                          {
//                            if(FaceNum==3)
//                                Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*(Grids->Variables.NSF[PCell].Face[4]);
//                            if(FaceNum==4)
//                                Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*(Grids->Variables.NSF[PCell].Face[3]);
//                            if(FaceNum==0)
//                                Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*(Grids->Variables.NSF[PCell].Face[1]+
//                                                                                  Grids->Variables.NSF[PCell].Face[2]);
//                            if(FaceNum==1)
//                                Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*(Grids->Variables.NSF[PCell].Face[0]+
//                                                                                  Grids->Variables.NSF[PCell].Face[2]);
//                            if(FaceNum==2)
//                                Grids->Variables.NSF[PCell].Face[FaceNum]=(-1.0)*(Grids->Variables.NSF[PCell].Face[0]+
//                                                                                  Grids->Variables.NSF[PCell].Face[1]);
//                          }
//                        double CellFlux=0.0;
//                        int i;
//                        for(i=0;i<Grids->Cells[PCell].NoOfFace;i++)
//                                    CellFlux+=Grids->Variables.NSF[PCell].Face[i];
//                        Grids->Variables.NSF[PCell].Face[FaceNum]=-CellFlux;
//
//                          if(G->Cells[PCell].CType==6)   ///Tet
//                          {
//
//                          }
//
//                          if(G->Cells[PCell].CType==7)   ///Pyramid
//                          {
//
//                          }
            }break;

            case 4:// Symmetry;
            {
                Grids->Variables.NSF[PCell].Face[FaceNum] = 0.0;
            }break;

            case 5:// Dirichlet Flow Function;
            {
                Area = Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                V = FlowFunction(FaceCenter,time);
                Grids->Variables.NSF[PCell].Face[FaceNum] = Dot(Area,V);
            }break;

            case 6:// Static Wall;
            {
                Grids->Variables.NSF[PCell].Face[FaceNum] = 0.0;;
            }break;
        }
    }
}

void SetOutletfFluxBoundary(struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector Area,FaceCenter,V,W,PC,PCF;
    double Value,dn,Densf,Amag;
    for(k=0;k<Grids->NB;k++)
    {
        PCell   = Grids->BoundaryCells[k].CellNum;
        Type    = Grids->BoundaryCells[k].Type;
        FaceNum = Grids->BoundaryCells[k].FaceNum;
        Value   = Grids->BoundaryCells[k].Value;
        switch (Type)
        {
            case 3:
            {
                double CellFlux=0.0;
                int i;
                for(i=0;i<Grids->Cells[PCell].NoOfFace;i++)
                    if (Grids->Cells[PCell].Neighbor[i]>=0)
                        CellFlux+=Grids->Variables.NSF[PCell].Face[i];
                Grids->Variables.NSF[PCell].Face[FaceNum]=-CellFlux;
            }break;

        }
    }
}

void SetMomentumBC(struct MatrixCoefficient *USMC
                  ,struct MatrixCoefficient *VSMC
                  ,struct MatrixCoefficient *WSMC
                  ,double *UST
                  ,double *VST
                  ,double *WST
                  ,struct Grid *Grids)
{
    int PCell,k,Type,FaceNum;
    struct MyVector n,Area,FaceCenter,V,W,rC,Pf,PC;
    struct MyVector Projected,CellVel;
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
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);

                ///saman
                nu=nu+(Grids->Variables.NuTurb[PCell]);

                W=RigidBodyVel(Grids,FaceCenter);
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);

            //  Set Diffusive Flux Boundary
                USMC->Elem[PCell][0]+=nu/dn*Amag;
                VSMC->Elem[PCell][0]+=nu/dn*Amag;
                WSMC->Elem[PCell][0]+=nu/dn*Amag;

                UST[PCell]+=nu/dn*Amag*(W.X+Dot(Grids->Variables.ULSEGrad[PCell],rC));
                VST[PCell]+=nu/dn*Amag*(W.Y+Dot(Grids->Variables.VLSEGrad[PCell],rC));
                WST[PCell]+=nu/dn*Amag*(W.Z+Dot(Grids->Variables.WLSEGrad[PCell],rC));
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
                PC=Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
                W=RigidBodyVel(Grids,FaceCenter);
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);

                ///saman
                nu=nu+(Grids->Variables.NuTurb[PCell]);

                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);;

                //  Set Diffusive Flux Boundary
                UST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.ULSEGrad[PCell],PC);///OK
                VST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.VLSEGrad[PCell],PC);
                WST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.WLSEGrad[PCell],PC);

                //  Set Convective Flux Boundary
                USMC->Elem[PCell][0]+=Flux;///OK
                VSMC->Elem[PCell][0]+=Flux;
                WSMC->Elem[PCell][0]+=Flux;

                UST[PCell]+=-Flux*Dot(Grids->Variables.ULSEGrad[PCell],Pf);///Test
                VST[PCell]+=-Flux*Dot(Grids->Variables.VLSEGrad[PCell],Pf);
                WST[PCell]+=-Flux*Dot(Grids->Variables.WLSEGrad[PCell],Pf);
            }break;
            case 4:// Symmetry;
            {
            }break;
            case 5:// Dirichlet Flow Function;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);

                ///saman
                nu=nu+(Grids->Variables.NuTurb[PCell]);

                W=RigidBodyVel(Grids,FaceCenter);
                V=FlowFunction(FaceCenter,time);
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);

                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);
                //  Set Diffusive Flux Boundary
                USMC->Elem[PCell][0]+=nu/dn*Amag;///OK
                VSMC->Elem[PCell][0]+=nu/dn*Amag;
                WSMC->Elem[PCell][0]+=nu/dn*Amag;

                UST[PCell]+=nu/dn*Amag*(V.X-Dot(Grids->Variables.ULSEGrad[PCell],rC));///OK
                VST[PCell]+=nu/dn*Amag*(V.Y-Dot(Grids->Variables.VLSEGrad[PCell],rC));
                WST[PCell]+=nu/dn*Amag*(V.Z-Dot(Grids->Variables.WLSEGrad[PCell],rC));
                //  Set Convective Flux Boundary
                UST[PCell]+=-Flux*V.X;///OK
                VST[PCell]+=-Flux*V.Y;
                WST[PCell]+=-Flux*V.Z;
            }break;
            case 6:// Static Wall;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);

                ///saman
                nu=nu+(Grids->Variables.NuTurb[PCell]);
            //  Set Diffusive Flux Boundary
                USMC->Elem[PCell][0]+=nu/dn*Amag;
                VSMC->Elem[PCell][0]+=nu/dn*Amag;
                WSMC->Elem[PCell][0]+=nu/dn*Amag;

                UST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.ULSEGrad[PCell],rC);///OK (I Changed W.X- to W.X+)
                VST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.VLSEGrad[PCell],rC);
                WST[PCell]+=nu/dn*Amag*Dot(Grids->Variables.WLSEGrad[PCell],rC);

            }break;
            case 7:// Free Slip wall;
            {
                Area =  Grids->Cells[PCell].Area[FaceNum];
                FaceCenter = Grids->Cells[PCell].FaceCenter[FaceNum];
                Amag=Mag(Area);
                n=Normalized(Area);
                dn=Mag(Minus(FaceCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]));
                nu=ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);

                ///saman
                nu=nu+(Grids->Variables.NuTurb[PCell]);

                W=RigidBodyVel(Grids,FaceCenter);

                CellVel.X=Grids->Variables.NU[PCell];
                CellVel.Y=Grids->Variables.NV[PCell];
                CellVel.Z=Grids->Variables.NW[PCell];
                Projected=SVP(Dot(CellVel,n),n);
                V=Minus(CellVel,Projected);
                rC=Minus(Grids->Cells[PCell].CellCenter,Grids->Cells[PCell].ProjectedCenter[FaceNum]);
                Flux=Grids->Variables.NSF[PCell].Face[FaceNum]-Dot(Area,W);
                //  Set Diffusive Flux Boundary
                USMC->Elem[PCell][0]+=nu/dn*Amag;///OK
                VSMC->Elem[PCell][0]+=nu/dn*Amag;
                WSMC->Elem[PCell][0]+=nu/dn*Amag;

                UST[PCell]+=nu/dn*Amag*(V.X-Dot(Grids->Variables.ULSEGrad[PCell],rC));///OK
                VST[PCell]+=nu/dn*Amag*(V.Y-Dot(Grids->Variables.VLSEGrad[PCell],rC));
                WST[PCell]+=nu/dn*Amag*(V.Z-Dot(Grids->Variables.WLSEGrad[PCell],rC));
            }break;
        }
    }
}

void MallocMatrixes(struct Grid *Grids
                   ,QMatrix *LaspackMC
                   ,struct MatrixCoefficient *MC
                   ,struct MatrixCoefficient *UMC
                   ,struct MatrixCoefficient *VMC
                   ,struct MatrixCoefficient *WMC
                   ,double **UST
                   ,double **VST
                   ,double **WST
                   ,double **PST
                   ,Vector *ST
                   ,Vector *Fee)
{
    int NC=0,NCol=0,PCell;
    NC=Grids->NC;

    (*MC).Elem =(double**) malloc(sizeof(double*)*NC);
    (*UMC).Elem =(double**) malloc(sizeof(double*)*NC);
    (*VMC).Elem =(double**) malloc(sizeof(double*)*NC);
    (*WMC).Elem =(double**) malloc(sizeof(double*)*NC);

    for(PCell=0;PCell<NC;PCell++)
    {
        NCol=Grids->Cells[PCell].NoOfFace+1;
        (*MC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
        (*UMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
        (*VMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
        (*WMC).Elem[PCell]=(double*)malloc(sizeof(double)*NCol);
    }

    *UST=(double*) malloc(sizeof(double)*NC);
    *VST=(double*) malloc(sizeof(double)*NC);
    *WST=(double*) malloc(sizeof(double)*NC);
    *PST=(double*) malloc(sizeof(double)*NC);

    Q_Constr(LaspackMC,"LaspackMC",NC,False,Rowws,Normal,True);

    V_Constr(ST,"ST",NC,Normal,True);

    V_Constr(Fee,"Fee",NC,Normal,True);
}


void ResetLaspackMatrix(QMatrix *LaspackMC,struct Grid *Grids)
{
    int NC;
    NC=Grids->NC;
    Q_Destr(LaspackMC);
    Q_Constr(LaspackMC,"GMC",NC,False,Rowws,Normal,True);
}

void FreeMatrixes(struct Grid *Grids
                 ,QMatrix *LaspackMC
                 ,struct MatrixCoefficient *MC
                 ,struct MatrixCoefficient *UMC
                 ,struct MatrixCoefficient *VMC
                 ,struct MatrixCoefficient *WMC
                 ,double **UST
                 ,double **VST
                 ,double **WST
                 ,double **PST
                 ,Vector *ST
                 ,Vector *Fee)
{
    int NC=0,PCell;
    NC=Grids->NC;
    for(PCell=0;PCell<NC;PCell++)
    {
        free((*MC).Elem[PCell]);
        free((*UMC).Elem[PCell]);
        free((*VMC).Elem[PCell]);
        free((*WMC).Elem[PCell]);
    }

    free((*MC).Elem);
    free((*UMC).Elem);
    free((*VMC).Elem);
    free((*WMC).Elem);

    free(*UST);
    free(*VST);
    free(*WST);
    free(*PST);

    Q_Destr(LaspackMC);

    V_Destr(ST);

    V_Destr(Fee);
}

void VelsSurfIntegral(double *VelSI,struct Grid *Grids)
{
    struct MyVector Hf,Area;
    int i,PCell;
    int NoOfFace;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        VelSI[PCell]=0.0;
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        for (i=0;i<NoOfFace;i++)
        {
            Hf.X=FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,i);
            Hf.Y=FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,i);
            Hf.Z=FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,i);
            Area=Grids->Cells[PCell].Area[i];
            VelSI[PCell] += ((Grids->Cells[PCell].Neighbor[i]>=0)?(Dot(Hf,Area)):(Grids->Variables.NSF[PCell].Face[i]));
        }
    }
}

void MakePressureMatrix(struct MatrixCoefficient *SMC,double *PST,struct Grid *Grids,double *VelSI,double alfa)
{
    int PCell;
    int NoOfFace,s;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        double ad[8];
        SetSVect(ad,8,0.0);
        FindPressureCoeff(ad,Grids,PCell);

        SMC->Elem[PCell][0]=ad[0]/alfa;

        for (s=0;s<NoOfFace;s++)
            if(Grids->Cells[PCell].Neighbor[s]>=0)
                SMC->Elem[PCell][s+1]=ad[s+1];

        PST[PCell]=-ad[NoOfFace+1]+VelSI[PCell]/dt+((1.-alfa)/alfa)*(ad[0]*Grids->Variables.NP[PCell]);
    }
}

void UpdateVelocities_UStar(struct Grid *Grids)
{
    int PCell;
    double PDens;
    struct MyVector GradP,Vel,Added;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
        PDens=DensCal(Grids->Variables.PAlpha[PCell]);
        GradP=Gradient_Disc_AxNode(Grids->Variables.P,Grids->Variables.PressLSEGrad,PCell,Grids);
//        GradP=Gradient(Grids->Variables.P,PCell,Grids);

//        Added=SVP(dt/PDens,GradP);
//        PDens=SelectDens(Vel,Added,PDens);              /// Only In The Case Of Two Phase Flow...

        Grids->Variables.NU[PCell] += GradP.X*dt/PDens;
        Grids->Variables.NV[PCell] += GradP.Y*dt/PDens;
        Grids->Variables.NW[PCell] += GradP.Z*dt/PDens;
    }
}

void UpdateVelocities_UN(struct Grid *Grids)
{
    int PCell;
    double PDens;
    struct MyVector GradP,Added,Vel;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        Vel.X=Grids->Variables.U[PCell];Vel.Y=Grids->Variables.V[PCell];Vel.Z=Grids->Variables.W[PCell];
        PDens=DensCal(Grids->Variables.NAlpha[PCell]);
        GradP=Gradient_Disc_AxNode(Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,Grids);
//        GradP=Gradient(Grids->Variables.NP,PCell,Grids);

//        Added=SVP(dt/PDens,GradP);
//        PDens=SelectDens(Vel,Added,PDens);       ///        Only In The Case Of Two Phase Flow...

        Grids->Variables.NU[PCell] += -GradP.X*dt/PDens;
        Grids->Variables.NV[PCell] += -GradP.Y*dt/PDens;
        Grids->Variables.NW[PCell] += -GradP.Z*dt/PDens;
    }
}

double MassErrorCalc(struct Grid *Grids)
{
    double CellFlux=0.0
          ,TotalFlux=0.0;
    int PCell,j,NoOfFace;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
        NoOfFace=Grids->Cells[PCell].NoOfFace;
        CellFlux=0.0;
        for (j=0;j<NoOfFace;j++)
            CellFlux+=Grids->Variables.NSF[PCell].Face[j];
        TotalFlux+=fabs(CellFlux);
    }
    return TotalFlux;

}

void InitialConditionForFlowVariables(struct Grid *Grids)
{
    int PCell,FaceNum;
    struct MyVector V,Area,CellCenter;
    double A,VMag;
    for(PCell=0;PCell<Grids->NC;PCell++)
    {
//        V=FlowFunction(Grids->Cells[PCell].CellCenter,time);
        Grids->Variables.PU[PCell] = Grids->Variables.NU[PCell] = Grids->Variables.U[PCell] =0.0;
        Grids->Variables.PV[PCell] = Grids->Variables.NV[PCell] = Grids->Variables.V[PCell] =0.0;
        Grids->Variables.PW[PCell] = Grids->Variables.NW[PCell] = Grids->Variables.W[PCell] =0.0;

//        Grids->Variables.PPhi[PCell] =
//        Grids->Variables.NPhi[PCell] =
//        Grids->Variables.Phi[PCell]  = Grids->Cells[PCell].CellCenter.Z-1.615;

        CellCenter=Grids->Cells[PCell].CellCenter;
        Grids->Variables.PAlpha[PCell] =
        Grids->Variables.NAlpha[PCell] =
        Grids->Variables.Alpha[PCell] =1.0;//(((CellCenter.Z<-0.01))?1.0:0.0);

        Grids->Variables.NP[PCell] = Grids->Variables.P[PCell] = 0.0;//((Grids->Cells[PCell].CellCenter.Z<0.0)
//                                                                 ?(Grids->Cells[PCell].CellCenter.Z*Gravity.Z*DensW)
//                                                                 :(Grids->Cells[PCell].CellCenter.Z*Gravity.Z*DensA));
         ///saman
        VMag=Mag(V);
        Grids->Variables.PK[PCell] = Grids->Variables.NK[PCell] = Grids->Variables.K[PCell] =0.0;//1.5*pow(Ti*VMag,2.0);
//
//        Grids->Variables.NuTurb[PCell] =0.0;// Tm*ViscCal(Grids->Variables.PAlpha[PCell])/DensCal(Grids->Variables.PAlpha[PCell]);
//
        Grids->Variables.PEps[PCell] =
        Grids->Variables.NEps[PCell] =
        Grids->Variables.Eps[PCell] =0.0;//0.1643*pow(Grids->Variables.K[PCell],1.5)*2.0/(0.07);
//        Grids->Variables.Eps[PCell] =0.09*pow(1.5*pow(Ti*VMag,2.0),2.0)/(Grids->Variables.NuTurb[PCell]);
        Grids->Variables.POm[PCell] =
        Grids->Variables.NOm[PCell] =
        Grids->Variables.Om[PCell] =0.0;
        Grids->Variables.NuTurb[PCell] =0.0;// 0.09*pow(Grids->Variables.K[PCell],2)/Grids->Variables.Eps[PCell];
    }

///     UpdatePressLSEGradient(Grids);
//     for(PCell=0;PCell<Grids->NC;PCell++)
//        {
//          Grids->Variables.PressLSEGrad[PCell]=Gradient_Disc_AxNode(Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,Grids);
////            Grids->Variables.PressLSEGrad[PCell]=Gradient(Grids->Variables.NP,PCell,Grids);
//        }
//
////    UpdatePhiLSEGradient(Grids);
//    UpdateVelsLSEGradient(Grids);
//
//    for(PCell=0;PCell<Grids->NC;PCell++)
//    {
//        int NoOfFace;
//        NoOfFace=Grids->Cells[PCell].NoOfFace;
//        for(FaceNum=0;FaceNum<NoOfFace;FaceNum++)
//        {
//        Area=Grids->Cells[PCell].Area[FaceNum];
//        V.X = FaceValue_AxNode(Grids->Variables.NU,Grids->Variables.ULSEGrad,Grids,PCell,FaceNum);
//        V.Y = FaceValue_AxNode(Grids->Variables.NV,Grids->Variables.VLSEGrad,Grids,PCell,FaceNum);
//        V.Z = FaceValue_AxNode(Grids->Variables.NW,Grids->Variables.WLSEGrad,Grids,PCell,FaceNum);
//        Grids->Variables.SF[PCell].Face[FaceNum]=Grids->Variables.NSF[PCell].Face[FaceNum]=Dot(Area,V);
//        Grids->Variables.AF[PCell].Face[FaceNum]=Grids->Variables.NAF[PCell].Face[FaceNum]=FaceValue(Grids->Variables.NAlpha,Grids,PCell,FaceNum);
//        }
//    }
//    SetSurfFluxBoundary(Grids);
}

void SetDamping(struct Grid *Grids,double *UST,double *VST,double *WST,
                                   double *u  ,double *v  ,double *w)
{
    double DampTermu,DampTermv,DampTermw;
    int PCell;

    for (PCell=0;PCell<Grids->NC;PCell++)
    {
        double R1,R2,R3,DampTerm=0.0,alpha=-50.0,Vol;

        R1=Grids->Cells[PCell].CellCenter.X-Grids->DynaVariables.MassCenter.X;
        R2=0.0;
        R3=5.0;
        alpha=-10.0;
        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));

//        R1=Grids->Cells[PCell].CellCenter.X-Grids->DynaVariables.MassCenter.X;
//        R2=-4.0;
//        R3=-5.0;
//        alpha=-1.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Y-Grids->DynaVariables.MassCenter.Y;
//        R2=0.0;
//        R3=0.5;
//        alpha=-1.0;
//        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Y-Grids->DynaVariables.MassCenter.Y;
//        R2=0.0;
//        R3=-0.5;
//        alpha=-1.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Z-Grids->DynaVariables.MassCenter.Z;
//        R2=2.0;
//        R3=3.0;
//        alpha=-1.0;
//        DampTerm+=alpha*((R1<R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));
//
//        R1=Grids->Cells[PCell].CellCenter.Z-Grids->DynaVariables.MassCenter.Z;
//        R2=-1.0;
//        R3=-2.0;
//        alpha=-1.0;
//        DampTerm+=alpha*((R1>R2)?(0.0):(((R1-R2)/(R3-R2))*((R1-R2)/(R3-R2))));

        Vol=Grids->Cells[PCell].Volume;
        DampTermu=DampTermv=DampTermw=DampTerm;

        UST[PCell]+=DampTermu*u[PCell]*Vol;
        VST[PCell]+=DampTermv*v[PCell]*Vol;
        WST[PCell]+=DampTermw*w[PCell]*Vol;

    }
}

void SolveNSEquation(struct Grid *Grids
                    ,char ConvectionMethod
                    ,char TemporalTermDiscretizationMethodForMomentum
                    ,char TimeDiscretizationMethodForConvection
                    ,char TimeDiscretizationMethodForDiffusion)
{
    struct MatrixCoefficient USMC,VSMC,WSMC,MC;
    double *UST,*VST,*WST,*PST;
    double *VelSI;
    int i,j,dummyCNT,dInt,PCell,c;

    QMatrix LaspackMC;
    Vector ST;
    Vector Fee;
    MallocMatrixes(Grids,&LaspackMC,&MC,&USMC,&VSMC,&WSMC,&UST,&VST,&WST,&PST,&ST,&Fee);
    VelSI=(double *)malloc(sizeof(double)*Grids->NC);

    i=0;
    do
    {
        UpdateVelsLSEGradient(Grids);
        dInt=0;
        i++;
        MakeMomentumMatrix(&USMC,UST,Grids,'U',ConvectionMethod,
                            TemporalTermDiscretizationMethodForMomentum,
                            TimeDiscretizationMethodForConvection,
                            TimeDiscretizationMethodForDiffusion);

        MakeMomentumMatrix(&VSMC,VST,Grids,'V',ConvectionMethod,
                            TemporalTermDiscretizationMethodForMomentum,
                            TimeDiscretizationMethodForConvection,
                            TimeDiscretizationMethodForDiffusion);

        MakeMomentumMatrix(&WSMC,WST,Grids,'W',ConvectionMethod,
                            TemporalTermDiscretizationMethodForMomentum,
                            TimeDiscretizationMethodForConvection,
                            TimeDiscretizationMethodForDiffusion);

        SetSourceTerms(UST,VST,WST,Grids,TemporalTermDiscretizationMethodForMomentum);
        SetMomentumBC(&USMC,&VSMC,&WSMC,UST,VST,WST,Grids);
//        SetDamping(Grids,UST,VST,WST,Grids->Variables.U,Grids->Variables.V,Grids->Variables.W);


//        for(c=0;c<Grids->NC;++c)
//        {
//            printf("\n");
//            printf("U[ %d , %d ] =     %f\n",c,c,USMC.Elem[c][0]);
//
//            for(i=0;i<Grids->Cells[c].NoOfFace;i++)
//                if(Grids->Cells[c].Neighbor[i]>=0)
//                    printf("U[ %d , %d ] =     %f\n",c,Grids->Cells[c].Neighbor[i],USMC.Elem[c][i+1]);
//            printf("UST[ %d ] =     %f\n",c,UST[c]);
//        }
//        getchar();

        SetRTCAccuracy(1e-10);

        ResetLaspackMatrix(&LaspackMC,Grids);
        CopyToLaspackMatrix(&LaspackMC,&USMC,Grids);
        CopyToLaspackVector(&ST ,UST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NU,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
        dInt+=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NU,&Fee,Grids);

        ResetLaspackMatrix(&LaspackMC,Grids);
        CopyToLaspackMatrix(&LaspackMC,&VSMC,Grids);
        CopyToLaspackVector(&ST ,VST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NV,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
        dInt+=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NV,&Fee,Grids);

        ResetLaspackMatrix(&LaspackMC,Grids);
        CopyToLaspackMatrix(&LaspackMC,&WSMC,Grids);
        CopyToLaspackVector(&ST ,WST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NW,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,100,ILUPrecond,1.0);
        dInt+=GetLastNoIter();
        CopyFromLaspackVector(Grids->Variables.NW,&Fee,Grids);

        printf("##################  Momentum Solver Iterations %d\n",dInt);
    }while((dInt>0)&&(i<8));

    UpdateVelocities_UStar(Grids);
    UpdateVelsLSEGradient(Grids);
    VelsSurfIntegral(VelSI,Grids);

    i=0;
    do
    {
        i++;
        ResetLaspackMatrix(&LaspackMC,Grids);
        MakePressureMatrix(&MC,PST,Grids,VelSI,1.0-1e-10);
        SetRTCAccuracy(1e-10);
//        SetRTCAccuracy(1e-4/pow(10.0,(i<3)?i:3));
        CopyToLaspackMatrix(&LaspackMC,&MC,Grids);
        CopyToLaspackVector(&ST,PST,Grids);
        CopyToLaspackVector(&Fee,Grids->Variables.NP,Grids);
        Fee=*BiCGSTABIter(&LaspackMC,&Fee,&ST,1500,ILUPrecond,1.5);
//        Fee=*CGIter(&LaspackMC,&Fee,&ST,1500,ILUPrecond,1.5);
        dInt=GetLastNoIter();
        printf("##################  Pressure Solver Iterations %d with Error %e\n",dInt,GetLastAccuracy());
        double dummy=0.0;
        dummy = V_GetCmp(&Fee,1);
        CopyFromLaspackVectorByOffset(Grids->Variables.NP,&Fee,dummy,Grids);

//        UpdatePressLSEGradient(Grids);
        ///Two phase flow #4

        for(PCell=0;PCell<Grids->NC;PCell++)
        {
//            Grids->Variables.PressLSEGrad[PCell]=Gradient_Disc_AxNode(Grids->Variables.NP,Grids->Variables.PressLSEGrad,PCell,Grids);
            ///Two phase Flow standard (faster) #3

//            Grids->Variables.PressLSEGrad[PCell]=LSEGradient(Grids->Variables.NP,PCell,Grids);
            ///One phase flow (least square)  #2

            Grids->Variables.PressLSEGrad[PCell]=Gradient(Grids->Variables.NP,PCell,Grids);
            ///one phase flow standard (faster) #1
        }

    }while((dInt>10)&&(i<10));

    CompSurfFlux(Grids);
    SetSurfFluxBoundary(Grids);
//    SetOutletfFluxBoundary(Grids);
    UpdateVelocities_UN(Grids);

   double TotalFlux=0.0,CellFlux;
   for(PCell=0;PCell<Grids->NC;PCell++)
    {
        CellFlux=0.0;
        for(i=0;i<Grids->Cells[PCell].NoOfFace;i++)
                CellFlux+=Grids->Variables.NSF[PCell].Face[i];

        TotalFlux+=fabs(CellFlux);
    }
    printf("Total Flux = %e\n",TotalFlux);

    free(VelSI);
    FreeMatrixes(Grids,&LaspackMC,&MC,&USMC,&VSMC,&WSMC,&UST,&VST,&WST,&PST,&ST,&Fee);
}
