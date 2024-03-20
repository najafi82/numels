void WriteToCDF(struct Grid *Grids,char *FileName)
{
    int j;
    char path[80];
    FILE *dummyCDF;
    sprintf(path,"./result/%s.cdf",FileName);
    if( (dummyCDF=fopen(path,"wb"))==NULL )
    {
        puts("CDF Writing error.................Exiting ( Write2CDF() 000)");
        exit(-1);
    }
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////Writing Constants///////////////////////////////
/// ///////////////////////////////////////////////////////////////////
    fwrite(&geps,sizeof(double),1,dummyCDF);
    fwrite(&dt,sizeof(double),1,dummyCDF);
    fwrite(&time,sizeof(double),1,dummyCDF);
    fwrite(&FinalTime,sizeof(double),1,dummyCDF);
    fwrite(&DensW,sizeof(double),1,dummyCDF);
    fwrite(&DensA,sizeof(double),1,dummyCDF);
    fwrite(&ViscW,sizeof(double),1,dummyCDF);
    fwrite(&ViscA,sizeof(double),1,dummyCDF);
    fwrite(&Gravity,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&pi,sizeof(double),1,dummyCDF);
    fwrite(&VoFMethod,sizeof(char),1,dummyCDF);
    fwrite(&ConvectionMethod,sizeof(char),1,dummyCDF);
    fwrite(&TimeDiscretizationMethodForConvection,sizeof(char),1,dummyCDF);
    fwrite(&TimeDiscretizationMethodForDiffusion,sizeof(char),1,dummyCDF);
    fwrite(&TimeDiscretizationMethodForVoF,sizeof(char),1,dummyCDF);
    fwrite(&TemporalTermDiscretizationMethodForMomentum,sizeof(char),1,dummyCDF);
    fwrite(&TemporalTermDiscretizationMethodForVoF,sizeof(char),1,dummyCDF);
    fwrite(&NoOfTimeStep,sizeof(int),1,dummyCDF);
    fwrite(&SaveCounter,sizeof(int),1,dummyCDF);
/// /////////////////////////////////////////////////////////////////
/// ////////////////////Writing Grid/////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/////////////////////////Writing Nodes/////////////////////////////////
    for(j=0;j<Grids->NN;j++)
        fwrite(&(Grids->Nodes[j].Pos),sizeof(struct MyVector),1,dummyCDF);
/////////////////////////Writing Cells/////////////////////////////////
    for(j=0;j<Grids->NC;j++)
    {
        fwrite(&(Grids->Cells[j].CellCenter),sizeof(struct MyVector),1,dummyCDF);
        fwrite(Grids->Cells[j].FaceCenter,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Cells[j].FaceCoincide,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Cells[j].ProjectedCenter,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Cells[j].Area,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
    }
/// //////////////////////Writing Flow Variables/////////////////////////////////
    fwrite(Grids->Variables.NU,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NV,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NW,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.U,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.V,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.W,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.PU,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.PV,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.PW,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.Alpha,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NAlpha,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.PAlpha,sizeof(double),Grids->NC,dummyCDF);

    ///saman
    fwrite(Grids->Variables.PEps,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NEps,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.Eps,sizeof(double),Grids->NC,dummyCDF);

    fwrite(Grids->Variables.POm,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NOm,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.Om,sizeof(double),Grids->NC,dummyCDF);

    fwrite(Grids->Variables.PK,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NK,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.K,sizeof(double),Grids->NC,dummyCDF);

    fwrite(Grids->Variables.NuTurb,sizeof(double),Grids->NC,dummyCDF);

    fwrite(Grids->Variables.KLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.EpsLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);

//    fwrite(Grids->Variables.Phi,sizeof(double),Grids->NC,dummyCDF);
//    fwrite(Grids->Variables.NPhi,sizeof(double),Grids->NC,dummyCDF);
//    fwrite(Grids->Variables.PPhi,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.P,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.NP,sizeof(double),Grids->NC,dummyCDF);
//    fwrite(Grids->Variables.Dummy,sizeof(double),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.ULSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.VLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.WLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fwrite(Grids->Variables.PressLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
//    fwrite(Grids->Variables.PhiLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    for(j=0;j<Grids->NC;j++)
    {
        fwrite(Grids->Variables.SF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Variables.NSF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Variables.AF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fwrite(Grids->Variables.NAF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
    }
//////////////////////////////////////////////////////////////////////
    fwrite(&Grids->DynaVariables.LinAccel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewLinAccel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.AngAccel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewAngAccel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.LinVel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewLinVel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.AngVel,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewAngVel,sizeof(struct MyVector),1,dummyCDF);

    fwrite(&Grids->DynaVariables.MassCenter,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewMassCenter,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.BodyOr,sizeof(struct MyVector),1,dummyCDF);
    fwrite(&Grids->DynaVariables.NewBodyOr,sizeof(struct MyVector),1,dummyCDF);
//////////////////////////////////////////////////////////////////////
/////////////////////////Writing Ended.///////////////////////////////
//////////////////////////////////////////////////////////////////////
    fclose(dummyCDF);
    printf("data.cdf Writed Succesfully.\n");
}

void ReadFromCDF(struct Grid *Grids,char *FileName)
{
    int j;
    char path[80];
    FILE *dummyCDF;
    sprintf(path,"./result/%s.cdf",FileName);
    if( (dummyCDF=fopen(path,"rb"))==NULL )
    {
        puts("CDF Reading error.................Exiting ( ReadFromCDF() 000)");
        exit(-1);
    }
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////Reading Constants///////////////////////////////
/// ///////////////////////////////////////////////////////////////////
    fread(&geps,sizeof(double),1,dummyCDF);
    fread(&dt,sizeof(double),1,dummyCDF);
    fread(&time,sizeof(double),1,dummyCDF);
    fread(&FinalTime,sizeof(double),1,dummyCDF);
    fread(&DensW,sizeof(double),1,dummyCDF);
    fread(&DensA,sizeof(double),1,dummyCDF);
    fread(&ViscW,sizeof(double),1,dummyCDF);
    fread(&ViscA,sizeof(double),1,dummyCDF);
    fread(&Gravity,sizeof(struct MyVector),1,dummyCDF);
    fread(&pi,sizeof(double),1,dummyCDF);
    fread(&VoFMethod,sizeof(char),1,dummyCDF);
    fread(&ConvectionMethod,sizeof(char),1,dummyCDF);
    fread(&TimeDiscretizationMethodForConvection,sizeof(char),1,dummyCDF);
    fread(&TimeDiscretizationMethodForDiffusion,sizeof(char),1,dummyCDF);
    fread(&TimeDiscretizationMethodForVoF,sizeof(char),1,dummyCDF);
    fread(&TemporalTermDiscretizationMethodForMomentum,sizeof(char),1,dummyCDF);
    fread(&TemporalTermDiscretizationMethodForVoF,sizeof(char),1,dummyCDF);
    fread(&NoOfTimeStep,sizeof(int),1,dummyCDF);
    fread(&SaveCounter,sizeof(int),1,dummyCDF);
/// ///////////////////////////////////////////////////////////////////
/// /////////////////////Reading Grid//////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
///////////////////////Reading Nodes//////////////////////////////////
    for(j=0;j<Grids->NN;j++)
        fread(&(Grids->Nodes[j].Pos),sizeof(struct MyVector),1,dummyCDF);
///////////////////////Reading Cells//////////////////////////////////
    for(j=0;j<Grids->NC;j++)
    {
        fread(&(Grids->Cells[j].CellCenter),sizeof(struct MyVector),1,dummyCDF);
        fread(Grids->Cells[j].FaceCenter,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Cells[j].FaceCoincide,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Cells[j].ProjectedCenter,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Cells[j].Area,sizeof(struct MyVector),Grids->Cells[j].NoOfFace,dummyCDF);
    }
/// ///////////////////////////////////////////////////////////////////
/// /////////////////////Reading Flow Variables////////////////////////
/// ///////////////////////////////////////////////////////////////////
    fread(Grids->Variables.NU,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NV,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NW,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.U,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.V,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.W,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.PU,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.PV,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.PW,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.Alpha,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NAlpha,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.PAlpha,sizeof(double),Grids->NC,dummyCDF);

    ///saman
    fread(Grids->Variables.PEps,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NEps,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.Eps,sizeof(double),Grids->NC,dummyCDF);

    fread(Grids->Variables.POm,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NOm,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.Om,sizeof(double),Grids->NC,dummyCDF);

    fread(Grids->Variables.PK,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NK,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.K,sizeof(double),Grids->NC,dummyCDF);

    fread(Grids->Variables.NuTurb,sizeof(double),Grids->NC,dummyCDF);

    fread(Grids->Variables.KLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fread(Grids->Variables.EpsLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);

//    fread(Grids->Variables.Phi,sizeof(double),Grids->NC,dummyCDF);
//    fread(Grids->Variables.NPhi,sizeof(double),Grids->NC,dummyCDF);
//    fread(Grids->Variables.PPhi,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.P,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.NP,sizeof(double),Grids->NC,dummyCDF);
//    fread(Grids->Variables.Dummy,sizeof(double),Grids->NC,dummyCDF);
    fread(Grids->Variables.ULSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fread(Grids->Variables.VLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fread(Grids->Variables.WLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
    fread(Grids->Variables.PressLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
//    fread(Grids->Variables.PhiLSEGrad,sizeof(struct MyVector),Grids->NC,dummyCDF);
//    Grids->Variables.SF=(struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
//    Grids->Variables.NSF=(struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
//    Grids->Variables.AF=(struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
//    Grids->Variables.NAF=(struct SurfFlux *)malloc(sizeof(struct SurfFlux)*Grids->NC);
    for(j=0;j<Grids->NC;j++)
    {
        Grids->Variables.SF[j].Face=(double *)malloc(sizeof(double)*Grids->Cells[j].NoOfFace);
        Grids->Variables.NSF[j].Face=(double *)malloc(sizeof(double)*Grids->Cells[j].NoOfFace);
        Grids->Variables.AF[j].Face=(double *)malloc(sizeof(double)*Grids->Cells[j].NoOfFace);
        Grids->Variables.NAF[j].Face=(double *)malloc(sizeof(double)*Grids->Cells[j].NoOfFace);
        fread(Grids->Variables.SF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Variables.NSF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Variables.AF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
        fread(Grids->Variables.NAF[j].Face,sizeof(double),Grids->Cells[j].NoOfFace,dummyCDF);
    }
//////////////////////////////////////////////////////////////////////
    fread(&Grids->DynaVariables.LinAccel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewLinAccel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.AngAccel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewAngAccel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.LinVel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewLinVel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.AngVel,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewAngVel,sizeof(struct MyVector),1,dummyCDF);

    fread(&Grids->DynaVariables.MassCenter,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewMassCenter,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.BodyOr,sizeof(struct MyVector),1,dummyCDF);
    fread(&Grids->DynaVariables.NewBodyOr,sizeof(struct MyVector),1,dummyCDF);
//////////////////////////////////////////////////////////////////////
/////////////////////////Reading Ended.///////////////////////////////
//////////////////////////////////////////////////////////////////////
    fclose(dummyCDF);

    printf("data.cdf Readed succesfully.\n");

}
