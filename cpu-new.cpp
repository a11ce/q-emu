#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Complex {
    public:
    float x, y;
    Complex(float _x, float _y) {
        x = _x;
        y = _y;
    }
    std::string toString() { return to_string(x) + "+" + to_string(y) + "i"; }
};

// Objects
using Matrix = std::vector<std::vector<Complex>>;

using StateVector = vector<Complex>;

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

// Matrix tensorSeries(vector<)

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

// Basic Gate Operations
static const Matrix I = {{Complex(1, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(1, 0)}};
static const Matrix H = {{Complex(1, 0), Complex(1, 0)},
                         {Complex(1, 0), Complex(-1, 0)}};
static const Matrix X = {{Complex(0, 0), Complex(1, 0)},
                         {Complex(1, 0), Complex(0, 0)}};

//Projection Operations
static const Matrix Proj_0 = {{Complex(1, 0), Complex(0, 0)},
                              {Complex(0, 0), Complex(0, 0)}};
static const Matrix Proj_1 = {{Complex(0, 0), Complex(0, 0)},
                              {Complex(0, 0), Complex(1, 0)}};

// Classes
class Gate {
public:
  Matrix matrix;
  int is_control = 0;
  Gate(Matrix _matrix): matrix(_matrix) {};
  Gate(Matrix _matrix, int _is_control): matrix(_matrix), is_control(_is_control) {};

  virtual void toMatrix() {
    cout << "Base Class Gate Matrix";
  }

  void print() {
    printMatrix(matrix);
  }
};

// class Gate {
//     virtual void toMatrix() = 0;
//     virtual void print();
// }

class ControlledGate : public Gate {
public:
    int controlWireIdx;
    int targetWireIdx;

    void toMatrix() override {
        cout << "Controlled Gate Matrix";
    }

    ControlledGate(Matrix _matrix, int _controlWireIdx, int _targetWireIdx) : Gate(_matrix, 1) {
        controlWireIdx = _controlWireIdx;
        targetWireIdx = _targetWireIdx;
    }
};

class OneQubitGate : public Gate {
public:
  int wireIdx;
//   int symbol;
  OneQubitGate(Matrix _matrix, int _wireIdx) : Gate(_matrix, 0) {
    wireIdx = _wireIdx;
  };
};

class TimeSlice {
  public:
    vector<Gate> gates;
    int nQubits;
    Matrix matrix;

    //if including explicit identity gates
    TimeSlice(vector<Gate> _gates) {
        gates = _gates;
        nQubits = _gates.size();
    };

    // TimeSlice(vector<Gate> _gates, int nQubits) {
    //     gates = vector<Gate>(nQubits, I);
    //     for(Gate g : _gates){
    //         gates.at(g.)
    //     }
    // }

    Matrix at(int index) {
      return gates.at(index).matrix;
    }

    int size() {
      return gates.size();
    }

    //returns the combined unitary for this timeslice 
    Matrix toMatrix() {
        int n = gates.size();
        Matrix U = gates.at(0).matrix;
        if(n < 2){
            return U;
        } else {
            for(int i = 1; i < gates.size(); i++) {
                // auto gate = gates;
                // if(gates.at(i).is_control == 1) {
                //     U = controlledOperation(gates.at(i).matrix, gates.at(i))
                // }
                U = tensorProduct(U, gates.at(i).matrix);
            }
        //     for(auto gate : gates) { 
        //         if(gate.is_control == 1) {
        //             ControlledGate c_gate = (ControlledGate)gate;
        //             // u = controlledOperation(gate.matrix, gates.at(i))
        //             TimeSlice TS_Ctrl(vector<Gate>(nQubits, Gate(I)));
        //             TimeSlice TS_Targ(vector<Gate>(nQubits, Gate(I)));

        //             TS_Ctrl.gates.at(c_gate.controlWireIdx) = Gate(Proj_0);
        //             TS_Targ.gates.at(c_gate.controlWireIdx) = Gate(Proj_1);
        //             TS_Targ.gates.at(c_gate.targetWireIdx) = Gate(c_gate.matrix);

        //             U = matrixAddition(TS_Ctrl.toMatrix(), TS_Targ.toMatrix());
        //             // matrix = U;
        //             return U;
        //         } else {
        //             U = tensorProduct(U, gate.matrix);
        //         }
        //     }
        }
        matrix = U;
        return U;
    }
};

// vector<Matrix> collapseTimeslice(TimeSlice TS) {
//     vector<Matrix> ts_m;
//     for(Gate G: TS.gates) {
//         if(G.is_control == 1) {
            
//         } else {
//             ts_m.push(G.matrix);
//         }
//     }
// }

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

    Gate gateAt(int sliceIdx, int wireIdx) {
        return slices.at(sliceIdx).at(wireIdx);
    }

    // Transformation //

    StateVector applyTransformationToSlice(StateVector SV, int sliceIdx) {
        // StateVector SV_t(SV);
        // for(int i = 0; i < (sliceIdx + 1); i++) {
        //     StateVector t(SV_t);
        //     SV_t = matrixVectorMultiply(t, slices.at(i).toMatrix());
        // }
        // return SV_t;

        StateVector SV_t(SV);    
        for(int i = 0; i < (sliceIdx + 1); i++) {
            // StateVector t(SV_t);
            SV_t = matrixVectorMultiply(SV_t, slices.at(i).toMatrix());
        }
        return SV_t;
    }

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

// Matrix controlledOperation(Matrix u, int control_index, int target_index, int num_of_wires){
//   //this needs to be update but since we usually only perform these multi-qubit controlled operations in their own timeslice, it might be okay?
//   // Matrix I = {{Complex(1, 0), Complex(0, 0)},
//   //             {Complex(0, 0), Complex(1, 0)}};

//   // Timeslice TS_Ctrl = Timeslice(num_of_wires, I);
//   // Timeslice TS_Targ = Timeslice(num_of_wires, I);

//   TimeSlice TS_Ctrl(vector<Gate>(num_of_wires, Gate(I)));
//   TimeSlice TS_Targ(vector<Gate>(num_of_wires, Gate(I)));

//   // TS_Ctrl.gates = vector<Matrix>(num_of_wires, I);
//   // TS_Targ.gates = vector<Matrix>(num_of_wires, I);

//   TS_Ctrl.gates.at(control_index) = projectionMatrixOf(0);
//   TS_Targ.gates.at(control_index) = projectionMatrixOf(1);
//   TS_Targ.gates.at(target_index) = u;

//   // Matrix TS_Contr = matrixAddition(TS_Ctrl.unitary(), TS_Targ.unitary());
//   Matrix TS_Contr = matrixAddition(tensorSeries(TS_Ctrl), tensorSeries(TS_Targ));

//   return TS_Contr;
// }

int main(void) {
    /**/
    // Slices & Tensoring //
    TimeSlice TS_0({OneQubitGate(I, 0), OneQubitGate(I, 1), OneQubitGate(I, 2)});
    Matrix TestIII = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
    testMatrixEqual(TS_0.toMatrix(), TestIII, "I ⊗ I ⊗ I");

    TimeSlice TS_1({OneQubitGate(I, 0), OneQubitGate(I, 1), OneQubitGate(X, 2)});
    Matrix TestIIX = {{Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};

    
    testMatrixEqual(TS_1.toMatrix(), TestIIX, "I ⊗ I ⊗ X");

    // Circuits //

    Circuit Circ_1({TimeSlice({OneQubitGate(H, 0), OneQubitGate(I, 1)}),
                    TimeSlice({OneQubitGate(I, 0), OneQubitGate(H, 1)})}, 2);
    // Circ_1.print();

    auto SV_0 = Circ_1.applyTransformationToSlice(makeStateVector(2), 1);
    auto SV_1 = Circ_1.applyTransformations(makeStateVector(2));
    testStateVectorsEqual(SV_0, SV_1, "((I ⊗ H)((H ⊗ I)|000>))");
    // auto SV_1 = Circ_1.applyTransformationAtSlice(SV_0, 1);

    // cout << typeid(OneQubitGate(I, 0)).name() << "\n";
    // cout << abs(-10);
    //
    // ControlledGate CG(X, 1,2 );
    // TimeSlice TS_C({CG});
    // for(auto g : TS_C.gates) {
    //     g.toMatrix();
    // };
}
