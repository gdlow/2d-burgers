#ifndef CLASS_MODEL2P
#define CLASS_MODEL2P

#include <mpi.h>

/**
 * @class Model
 * @brief Sets up the model instance specifying key parameters constructing the problem
 * */
class Model {
public:
    explicit Model(int argc, char* argv[]);
    ~Model();

    void PrintParameters();

    bool IsValid();

    /// Generic getters
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
    double GetBDx()    const { return bdx; }
    double GetBDy()    const { return bdy; }
    double GetBetaDx_2() const { return beta_dx_2; }
    double GetBetaDy_2() const { return beta_dy_2; }
    double GetBetaDx_Sum() const { return beta_dx_sum; }
    double GetBetaDy_Sum() const { return beta_dy_sum; }
    double GetAlpha_Sum() const { return alpha_sum; }

    // Add any other getters here...

    /// MPI getters
    int GetRank()      const { return loc_rank; }
    int GetPx()        const { return Px; }
    int GetPy()        const { return Py; }
    int GetUp()        const { return up; }
    int GetDown()      const { return down; }
    int GetLeft()      const { return left; }
    int GetRight()     const { return right; }
    int GetLocNxr()    const;
    int GetLocNyr()    const;
    int GetLocNyrNxr() const;
    int GetDisplX()    const;
    int GetDisplY()    const;
    int* GetDispls()         { return displs; }
    int* GetRecvCount()      { return recvcount; }
    int* GetRankNxrMap()     { return rankNxrMap; }
    int* GetRankNyrMap()     { return rankNyrMap; }
    int* GetRankDisplsXMap() { return rankDisplsXMap; }
    int* GetRankDisplsYMap() { return rankDisplsYMap; }
    MPI_Comm GetComm()       { return vu; }

private:
    void ParseParameters(int argc, char* argv[]);
    void ValidateParameters();

    /// Private setters
    void SetNumerics();
    void SetGridParameters();
    void SetCartesianGrid();
    void SetNeighbours();

    bool verbose;
    bool help;

    // Numerics: Everything here has to be predefined

    /// ParseParameters
    double Lx;
    double Ly;
    double T;
    /// SetNumerics
    double x0;
    double y0;
    int    Nx;
    int    Ny;
    int    Nt;
    double dx;
    double dy;
    double dt;

    /// Physics
    double ax;
    double ay;
    double b;
    double c;

    /// Constants for Burger problem
    double bdx;
    double bdy;

    double beta_dx_2;
    double beta_dy_2;
    double beta_dy_sum;
    double beta_dx_sum;
    double alpha_sum;


    // Add any additional parameters here...

    /// MPI Parameters
    int p;
    int loc_rank;
    int Px;
    int Py;
    int* loc_coord;
    int* loc_Nxr;
    int* loc_Nyr;
    int* displs_x;
    int* displs_y;
    int* displs;
    int* recvcount;
    int* rankNxrMap;
    int* rankNyrMap;
    int* rankDisplsXMap;
    int* rankDisplsYMap;
    MPI_Comm vu;
    int up, down, left, right;
};

#endif //CLASS_MODEL2P