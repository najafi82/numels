#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))
// #define isinf(a) (((a)>(1e308))||((a)<(-1e307)))
#define sign(a) ((a)>=0)?((isequal(a,0)?0:1)):(-1)
#define Dot(a,b) ((a.X*b.X)+(a.Y*b.Y)+(a.Z*b.Z))
#define sqr(a) ((a)*(a))
#define isequal(a,b) (fabs((a)-(b))<=geps)

/// DensCal(),ViscCal() and MolcNu() calculate density ,dynamic viscosity and kinematic viscosity of effective fluid respectively
#define DensCal(a) (((a)*DensW)+((1-(a))*DensA))
#define ViscCal(a) (((a)*ViscW)+((1-(a))*ViscA))
#define MolcNu(a)  (ViscCal((a))/DensCal((a)))

#define LSE_NodeShare
//#define LSE_Neighbour

const int BREAK_NODES[6/*No of Faces*/][5/*Max nodes per face + 1*/]=
{
    {4,0,1,5,4},//Face 0
    {4,1,3,7,5},//Face 1
    {4,3,2,6,7},//Face 2
    {4,2,0,4,6},//Face 3
    {4,1,0,2,3},//Face 4
    {4,4,5,7,6} //Face 5
};

const int WEDGE_NODES[5/*No of Faces*/][5/*Max nodes per face + 1*/]=
{
    {4,0,1,4, 3},//Face 0
    {4,1,2,5, 4},//Face 1
    {4,2,0,3, 5},//Face 2
    {3,0,2,1,-1},//Face 3
    {3,3,4,5,-1},//Face 4
};

const int TETRA_NODES[4/*No of Faces*/][4/*Max nodes per face + 1*/]=
{
    {3,1,0,2},//Face 0
    {3,0,1,3},//Face 1
    {3,1,2,3},//Face 2
    {3,2,0,3},//Face 3
};

const int PYRAMID_NODES[5/*No of Faces*/][5/*Max nodes per face + 1*/]=
{
    {4,0,2,3, 1},//Face 0
    {3,0,1,4,-1},//Face 1
    {3,1,3,4,-1},//Face 2
    {3,3,2,4,-1},//Face 3
    {3,2,0,4,-1},//Face 4
};

const int GAMBIT_2_VTK_BREAK  [8] = {0,1,3,2,4,5,7,6};
const int GAMBIT_2_VTK_WEDGE  [6] = {0,2,1,3,5,4};
const int GAMBIT_2_VTK_TETRA  [4] = {0,3,1,2};
const int GAMBIT_2_VTK_PYRAMID[5] = {0,1,3,2,4};

struct MyVector  ///Catresian 3D Vector
{
    double X,Y,Z;  ///Catresian 3D Vector component
};

struct RawCell  ///cell is an arbitrary shape element (brick,tetrahedrom,pyramid or wedge)
{
    int Type;  ///Type of Element
    int NodePerCell;  ///Number of nodes makes Cell
    int *NodeList;  ///List of Nodes According to GAMBIT NEUTRAL File Foramt.
};

struct Cell ///
{
    int NoOfFace;  //ok
    int *Neighbor;  //ok
    struct MyVector CellCenter;  //  ok
    struct MyVector *FaceCenter;  //  ok Center of faces
    struct MyVector *FaceCoincide;  //  Coincidence of cellcenters line and face
    struct MyVector *ProjectedCenter;  //  Projection of cellcenter to normal area line
    struct MyVector *Area;
    double Volume;
};

struct Node
{
    struct MyVector Pos;///
    int NumOfSharedCells;  /// Indicates that how many cells are shared in a node
    int *SharedCells;  /// cell number list
};


struct BoundaryCell  ///indicate type and value of boundary on each face
{
    int CellNum;  ///Element number
    int Type;        ///Type of Boundary
    int FaceNum;     ///Face Number
    double Value;    ///Value related to Boundary,(This value is important for symmetry type boundary condition)
};

struct RigidBoundaryCell  ///indicate number of element and its corresponding face for rigid body boundary
{
    int CellNum;  ///Element number
    int FaceNum;     ///Face Number
};

struct SurfFlux ///Include surface fluxes for each element
{
    double *Face;///surface fluxes for each element
};

struct DonAcc  //Includes number of donner and acceptor Element(Cell) for a face
{
    int d //  donner   cell number
       ,a;//  acceptor cell number
};

struct DonnerAcceptor
{
    struct DonAcc Face[6];
};

struct PosValue
{
    struct MyVector Pos;
    double Value;
    double Weight;
};

struct PosValueList
{
    int size;
    struct PosValue *Elem;
};

struct FlowVariables
{
    double *NU;
    double *NV;
    double *NW;
    double *U;
    double *V;
    double *W;
    double *PU;
    double *PV;
    double *PW;
    double *Alpha;
    double *NAlpha;
    double *PAlpha;

    ///saman
    double *NEps;
    double *Eps;
    double *PEps;

    double *NK;
    double *K;
    double *PK;
    double *NOm;
    double *Om;
    double *POm;

    double *NuTurb;

    struct MyVector *EpsLSEGrad;
    struct MyVector *KLSEGrad;
    struct MyVector *OmLSEGrad;

//    double *Phi;
//    double *NPhi;
//    double *PPhi;
    double *P;
    double *NP;
//    double *Dummy;
    struct MyVector *PressLSEGrad;
    struct MyVector *ULSEGrad;
    struct MyVector *VLSEGrad;
    struct MyVector *WLSEGrad;
//    struct MyVector *PhiLSEGrad;
    struct SurfFlux *SF;
    struct SurfFlux *NSF;
    struct SurfFlux *AF;
    struct SurfFlux *NAF;
};

struct DynamicsVariables
{
    double Mass;
    double Inertia[3][3];

    struct MyVector NewLinAccel;
    struct MyVector NewAngAccel;
    struct MyVector NewLinVel;
    struct MyVector NewAngVel;
    struct MyVector NewMassCenter;
    struct MyVector NewBodyOr;

    struct MyVector LinAccel;
    struct MyVector AngAccel;
    struct MyVector LinVel;
    struct MyVector AngVel;
    struct MyVector MassCenter;
    struct MyVector BodyOr;
};

struct Grid
{
    int NC;
    int NN;
    int NB;
    int NRB;

    int NRBG;
    int *RBNo;

    struct Node *Nodes;
    struct RawCell *RCells;
    struct Cell *Cells;
    struct BoundaryCell *BoundaryCells;
    struct RigidBoundaryCell *RigidBoundaryCells;
    struct FlowVariables Variables;
    struct DynamicsVariables DynaVariables;
};

struct MatrixCoefficient  ///Includes matrix coefficients and source terms for a domain with Hexagonal Element meshs
{
    double **Elem;  ///matrix coefficients
};

double geps=1e-10;

double pi=3.14159265358979323846;

char VoFMethod = 'H'
               // C : CICSAM
               // H : HRIC
               // S : STACS
    ,ConvectionMethod = 'C'
    ,KEConvectionMethod = 'U'
    ,KOConvectionMethod = 'U'
                      // U : UDS
                      // C : CDS
                      // T : UDS + CDS
                      // D : UDS by Deffered CDS terms
                      // G : Jasak

    ,TemporalTermDiscretizationMethodForMomentum = 'A'
    ,TemporalTermDiscretizationMethodForVoF      = 'A'
                              // A : Two Time Level
                              // B : Three Time Level

    ,TimeDiscretizationMethodForConvection = 'I'
    ,TimeDiscretizationMethodForDiffusion = 'I'
    ,TimeDiscretizationMethodForVoF      = 'C';
                              // I : Implicit
                              // C : Cranck-Nicholson
                              // E : Explicit
