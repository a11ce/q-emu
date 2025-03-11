#include <iostream>
#include <vector>
#include <cmath>
#include <map>

using namespace std;

// Objects
class Complex {
    public:
    float x, y;
    Complex(float _x, float _y) {
        x = _x;
        y = _y;
    }
    std::string toString() { return to_string(x) + "+" + to_string(y) + "i"; }
};

using Matrix = std::vector<std::vector<Complex>>;

using StateVector = vector<Complex>;

// Basic Gate Operations
static const Matrix I = {{Complex(1, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(1, 0)}};
static const Matrix H = {{Complex(1, 0), Complex(1, 0)},
                         {Complex(1, 0), Complex(-1, 0)}};
  //this is not technically correct, needs the 1/sqrt(2) coeffec. but fine for testing purposes
static const Matrix X = {{Complex(0, 0), Complex(1, 0)},
                         {Complex(1, 0), Complex(0, 0)}};
static const Matrix Z = {{Complex(1, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(-1, 0)}};
static const Matrix S = {{Complex(1, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 1)}};

// Add still: T, Y, ...

//Projection Operations
static const Matrix Proj_0 = {{Complex(1, 0), Complex(0, 0)},
                              {Complex(0, 0), Complex(0, 0)}};
static const Matrix Proj_1 = {{Complex(0, 0), Complex(0, 0)},
                              {Complex(0, 0), Complex(1, 0)}};

//Helper Functions
void printStateVector(StateVector SV) {
  cout << "[\n";
  for (auto c : SV) {
    cout << "  " << c.toString() << "\n";
  }
  cout << "]\n";
}

void printMatrix(Matrix M) {
  cout << "[\n";
  for (auto r : M){
    for (auto c : r) {
      cout << "  " << c.toString();
    }
    cout << "  \n";
  }
  cout << "]\n";
}

Matrix initSquareMatrix(size_t size) {
  return Matrix(size, std::vector<Complex>(size, Complex(0,0)));
}

StateVector makeStateVector(size_t nQubits) {
  StateVector SV;
  Complex oneProb = Complex(1, 0);
  Complex zeroProb = Complex(0, 0);
  int size_sv = std::pow(2, nQubits);
  SV.push_back(oneProb);
  for (size_t idx = 0; idx < size_sv - 1; idx++) {
    SV.push_back(zeroProb);
  }
  return SV;
}

//Test Functions
int testMatrixEqual(Matrix M, Matrix M_expected, string test_name = "", int print_out = 0){
  //assuming square matrices (probably a better way to compare the matrixs but it is what it is)
  int passed = 1;
  if (M.size() != M_expected.size() || M.at(0).size() != M_expected.at(0).size()){
    passed = 0;
  } else {
    for(int r = 0; r < M.size(); r++) {
      for(int c = 0; c < M.size(); c++){
        if (M[r][c].x != M_expected[r][c].x){
          passed = 0;
        } else if (M[r][c].y != M_expected[r][c].y){
          passed = 0;
        }
      }
    }
  }
  if(print_out == 1) {
    printMatrix(M);
  }
  if(test_name != "") {
    cout << test_name << " => ";
  }
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}

int testStateVectorsEqual(StateVector SV, StateVector SV_expected, string test_name = "", int print_out = 0){
  int passed = 1;
  if (SV.size() != SV_expected.size()){
    passed = 0;
  } else {
    for(int r = 0; r < SV.size(); r++) {
      if (SV[r].x != SV_expected[r].x){
        passed = 0;
      } else if (SV[r].y != SV_expected[r].y){
        passed = 0;
      }
    }
  }
  if(print_out == 1 ){
    printStateVector(SV);
  }
  if(test_name != "") {
    cout << test_name << " => ";
  }
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}

//Operations
Matrix tensorProduct(Matrix M, Matrix N){
  //M tensor N = P

  //assuming M & N are square matrices; and M ⊗ N; 
  int size_M = M.size();
  int size_N = N.size();

  //init output matrix of size (size_M*size_N) x (size_M*size_N)
  // cout << "Size of matrix: " << size_M * size_N << " x " << size_M * size_N << "; \n";
  Matrix P = initSquareMatrix(size_M * size_N);

  //iter matrix M
  for(int col_M_i = 0; col_M_i < size_M; col_M_i++) {
    for(int row_M_i = 0; row_M_i < size_M; row_M_i++){
      Complex M_val = M[row_M_i][col_M_i];

      //iter matrix N
      for(int col_N_i = 0; col_N_i < size_N; col_N_i++) {
        for(int row_N_i = 0; row_N_i < size_N; row_N_i++) {
          Complex N_val = N[row_N_i][col_N_i];

          //find index in P
          // P Row = (Size N x Row Index M) + Row Index N
          // P Col = (Size N x Col Index M) + Col Index M
          int P_row = (size_N * row_M_i) + row_N_i;
          int P_col = (size_N * col_M_i) + col_N_i;
          P[P_row][P_col] = Complex(M_val.x * N_val.x, M_val.y * N_val.y);
        }
      }
    }
  }

  return P;
}

Matrix tensorSeries(vector<Matrix> M_Vec) {
    Matrix U = M_Vec.at(0);
    for(int i = 1; i < M_Vec.size(); i++) {
        U = tensorProduct(U, M_Vec.at(i));
    }
    return U;
}

StateVector matrixVectorMultiply(StateVector SV, Matrix M) { 
  //init output state vector to all zeros
  StateVector O_SV;
  int n_Rows_M = M.size();

  int n_Cols_M = M[0].size();

  for (int row_i = 0; row_i < n_Rows_M; row_i++) {
    Complex prod = Complex(0,0);
    //iter rows of each col in M
    for (int col_i = 0; col_i < n_Cols_M; col_i++) {
      //sparse matrix, add check for 0 in M[row_i][col_i], skip next for if 0 (fine here)
      //iter rows of SV vetor
      prod.x += M[row_i][col_i].x * SV[col_i].x;
      prod.y += M[row_i][col_i].y * SV[col_i].y;
    }
    O_SV.push_back(prod);
  }
  return O_SV;
}

Matrix matrixAddition(Matrix M, Matrix N) {
  //ensure that matrices are the same size & square (they should always be square in our case but still doesnt hurt to check)
  //TODO
  Matrix P = initSquareMatrix(M.size());
  for(int r = 0; r < M.size(); r++){
    for(int c = 0; c < M.size(); c++){
      P[r][c].x = M[r][c].x + N[r][c].x;
      P[r][c].y = M[r][c].y + N[r][c].y;
    }
  }
  return P;
}

/*
  the process here allows us to support any configuration of a single control/target controlled matrix.
  it has to be done in this goofy ass/mildly hardcoded way.

  it's somewhat of an odd process to create arbitrary controlled operations that scales exponentially with the # of control qubits you have, so for right 
  now (until I can figure out if it's feasible to offer inherent support for something like the Toffoli operation) any multi-control gate will have to be
  decomposed into a set of operations that are supported.
*/
Matrix collapseControlledMatrixVector(vector<Matrix> M_Vec, Matrix U, int control, int target) {
    vector<Matrix> TS_Ctrl(M_Vec);
    vector<Matrix> TS_Targ(M_Vec);
    //this essentially encodes two scenarios:
      //1.) when the control qubit is 0, the target function should be the identity (i.e. leave it alone)
      //2.) when the control qubit is 1, the target function should be the given matrix U
    //these two "scenarios" have to then be tensored and added together to correctly encode any 2 qubit controlled operation
    TS_Ctrl.at(control) = Proj_0;
    TS_Targ.at(control) = Proj_1;
    TS_Targ.at(target) = U;
    return matrixAddition(tensorSeries(TS_Ctrl), tensorSeries(TS_Targ));
}

// Classes
class Gate {
    public:
        struct Interface
        {
            virtual string type() const = 0;
            virtual Matrix toMatrix() const = 0;
            virtual int wire() const = 0;

            virtual ~Interface() = default;

            virtual int targetWire() const {
                return -1;
            }

            virtual void print() const {
                printMatrix(toMatrix());
            }
        };
        std::shared_ptr<const Interface> _p;
    public:
        Gate(Interface* p) : _p(p){
        }
        string type() const {
            return _p->type();
        }
        Matrix toMatrix() const {
            return _p->toMatrix();
        }
        int wire() const {
            return _p->wire();
        }
        int targetWire() const {
            return _p->targetWire();
        }
        void print() const {
            _p->print();
        }
};

class SingleQubitGate : public Gate::Interface {
    public: 
        Matrix matrix;
        int wireIdx;

        SingleQubitGate(Matrix _matrix, int _wireIdx) 
            : matrix(_matrix), wireIdx(_wireIdx) {};

        string type() const override {
            return "single";
        };

        Matrix toMatrix() const override {
            return matrix;
        };

        int wire() const override {
            return wireIdx;
        }
};

class ControlledGate : public Gate::Interface {
    public:
        Matrix matrix;
        int wireIdx;
        int targetWireIdx;

        ControlledGate(Matrix _matrix, int _wireIdx, int _targetWireIdx) 
            : matrix(_matrix), wireIdx(_wireIdx), targetWireIdx(_targetWireIdx) {};
        
        string type() const override {
            return "controlled";
        }

        Matrix toMatrix() const override {
            return matrix;
        }

        int wire() const override {
            return wireIdx;
        }

        int targetWire() const override {
            return targetWireIdx;
        }
};

class TimeSlice {
  public:
    vector<Gate> gates;
    vector<Matrix> matrices;
    int nQubits;

    TimeSlice(vector<Gate> _gates, int _nQubits) : gates(_gates), nQubits(_nQubits){};

    Matrix processGates(vector<Gate> _gates, int num_of_wires) {
      vector<Matrix> M_Vec(num_of_wires, I);
      //currently testing an implementation allowing for more than one controlled operation per slice
      vector<Gate> control_queue;
      for(Gate g : _gates) {
        // if gate is control (not best way to do this but it works for now)
        if(g.type() == "controlled"){
          control_queue.emplace_back(g);
        } else {
          M_Vec.at(g.wire()) = g.toMatrix();
        }
      }
      if(!control_queue.empty()) {
        // again, leaving this as a "queue"/vector for the above reasons even though it currently will only pull the first
        Gate g = control_queue.at(0);
        return collapseControlledMatrixVector(M_Vec, g.toMatrix(), g.wire(), g.targetWire());
      } else {
        return tensorSeries(M_Vec);
      }
    }

    /* currently tinkering with a better/more effic. method for providing the unitary here that will wrap the processGates function,
      hence the arbitrary function wrap here. */

    //returns the combined unitary for this timeslice 
    Matrix toMatrix() {
        return processGates(gates, nQubits);
    }
};

class Circuit {
  public:
    int nQubits;
    vector<TimeSlice> slices;

    Circuit(int _nQubits) {
        int nQubits = _nQubits;
    }

    Circuit(vector<TimeSlice> _slices, int _nQubits) {
      nQubits = _nQubits;
      slices = _slices;
    }

    // Transformation //
    // Allows to 
    StateVector applyTransformationToSlice(StateVector SV, int sliceIdx) {
        StateVector SV_t(SV);    
        for(int i = 0; i < (sliceIdx + 1); i++) {
            SV_t = matrixVectorMultiply(SV_t, slices.at(i).toMatrix());
        }
        return SV_t;
    }

    //Computes the final state vector of the circuit
    StateVector applyTransformations(StateVector _SV) {
        StateVector SV(_SV);
        for(TimeSlice TS : slices) {
            StateVector t(SV);
            SV = matrixVectorMultiply(t, TS.toMatrix());
        }
        return SV;
    }

    // Helpers // 
    void print() {
        int i = 1;
        for (TimeSlice TS : slices) {
            cout << "Timeslice " << i << ":\n";
            i+=1;
            for (Gate G : TS.gates) {
                G.print();
            };
        }
    }
};

int main(void) {
    // Slices & Tensoring //
    
    TimeSlice TS_0({}, 3);
    Matrix TestIII = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
    testMatrixEqual(TS_0.toMatrix(), TestIII, "I ⊗ I ⊗ I");
    
    TimeSlice TS_1({new SingleQubitGate(X, 2)}, 3);
    Matrix TestIIX = {{Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};

    
    testMatrixEqual(TS_1.toMatrix(), TestIIX, "I ⊗ I ⊗ X");

    TimeSlice TS_2({new SingleQubitGate(X, 0), new ControlledGate(X, 2, 1)}, 3);
    Matrix TestXICX = {{Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)}};
    testMatrixEqual(TS_2.toMatrix(), TestXICX, "X ⊗ CX(2, 1)");

    // Baby Circuits //
    Circuit Circ_1({TimeSlice({new SingleQubitGate(X, 0)}, 2),
                   TimeSlice({new SingleQubitGate(X, 1)}, 2)}, 2);
    StateVector Ket11 = {Complex(0,0), Complex(0,0), Complex(0,0),Complex(1,0)};
    testStateVectorsEqual(Circ_1.applyTransformations(makeStateVector(2)), Ket11, "(I ⊗ X((X ⊗ I)|00>))");

    TimeSlice Circ_2_TS_0({new SingleQubitGate(X, 1)}, 2);
    TimeSlice Circ_2_TS_1({new ControlledGate(X, 0, 1)}, 2);
    TimeSlice Circ_2_TS_2({new ControlledGate(X, 1, 0)}, 2);
    Circuit Circ_2({Circ_2_TS_0, Circ_2_TS_1, Circ_2_TS_2}, 2);
    StateVector Ket01 = {Complex(0,0), Complex(1,0), Complex(0,0),Complex(0,0)};
    testStateVectorsEqual(Circ_2.applyTransformationToSlice(makeStateVector(2),0), Ket01, "(I ⊗ X)|00>))");
    testStateVectorsEqual(Circ_2.applyTransformationToSlice(makeStateVector(2),1), Ket01, "(CX(0,1))|01>))");
    testStateVectorsEqual(Circ_2.applyTransformationToSlice(makeStateVector(2),2), Ket11, "(CX(1,0)|01>))");

}
