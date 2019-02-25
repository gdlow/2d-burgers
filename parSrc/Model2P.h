#ifndef CLASS_MODEL2P
#define CLASS_MODEL2P

#include <mpi.h>

class Model {
public:
    Model(int argc, char* argv[]);
    ~Model();

    void PrintParameters();

    bool IsValid();

    // Getters
    bool   IsVerbose() const { return verbose; }
    bool   IsHelp()    const { return help; }
    double GetX0()     const { return x0; }
    double GetY0()     const { return y0; }
    double GetLx()     const { return Lx; }
    double GetLy()     const { return Ly; }
    double GetT()      const { return T; }
    int    GetNx()     const { return Nx; }
    int    GetNy()     const { return Ny; }
    int    GetNt()     const { return Nt; }
    double GetDx()     const { return dx; }
    double GetDy()     const { return dy; }
    double GetDt()     const { return dt; }
    double GetAx()     const { return ax; }
    double GetAy()     const { return ay; }
    double GetB()      const { return b; }
    double GetC()      const { return c; }
    MPI_Comm GetComm()       { return vu; }

    // Add any other getters here...

    // MPI parameters
    int GetP()         const { return p; }
    int GetRank()      const { return loc_rank; }
    int GetPx()        const { return Px; }
    int GetPy()        const { return Py; }
    int GetUp()        const { return up; }
    int GetDown()      const { return down; }
    int GetLeft()      const { return left; }
    int GetRight()     const { return right; }
    int GetCoordX();
    int GetCoordY();
    int GetLocNxr();
    int GetLocNyr();
    int GetDisplX();
    int GetDisplY();

private:
    void ParseParameters(int argc, char* argv[]);
    void ValidateParameters();

    // Private Setters
    void SetNumerics();
    void SetGridParameters();
    void SetCartesianGrid();
    void SetNeighbours();

    bool verbose;
    bool help;

    // Numerics: Everything here has to be predefined

    // ParseParameters:
    double Lx;
    double Ly;
    double T;
    // SetNumerics:
    double x0;
    double y0;
    int    Nx;
    int    Ny;
    int    Nt;
    double dx;
    double dy;
    double dt;

    // Physics
    double ax;
    double ay;
    double b;
    double c;

    // Add any additional parameters here...

    // MPI Parameters
    int p;
    int loc_rank;
    int Px;
    int Py;
    int* loc_coord;
    int* loc_Nxr;
    int* loc_Nyr;
    int* displs_x;
    int* displs_y;
    MPI_Comm vu;
    int up, down, left, right;
};

#endif //CLASS_MODEL2P