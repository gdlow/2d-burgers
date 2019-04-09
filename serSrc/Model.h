#ifndef CLASS_MODEL
#define CLASS_MODEL

/**
 * @class Model
 * @brief Sets up the model instance specifying key parameters constructing the problem
 * */
class Model {
public:
    Model(int argc, char* argv[]);
    ~Model();

    void PrintParameters();

    bool IsValid();

    /// Getters
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

private:
    void ParseParameters(int argc, char* argv[]);
    void ValidateParameters();

    /// Private Setters
    void SetNumerics();

    bool verbose;
    bool help;

    /// Numerics: Everything here has to be predefined

    /// ParseParameters:
    double Lx;
    double Ly;
    double T;
    
    /// SetNumerics:
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
};

#endif //CLASS_MODEL