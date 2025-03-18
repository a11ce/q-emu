#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <vector>

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

const float INV_SQRT_2 = 0.70710678118654752440084436210484903928483593768847;

// Basic Gate Operations
const Matrix I = {{Complex(1, 0), Complex(0, 0)},
                  {Complex(0, 0), Complex(1, 0)}};

const Matrix H = {{Complex(INV_SQRT_2, 0), Complex(INV_SQRT_2, 0)},
                  {Complex(INV_SQRT_2, 0), Complex(-1 * INV_SQRT_2, 0)}};

const Matrix X = {{Complex(0, 0), Complex(1, 0)},
                  {Complex(1, 0), Complex(0, 0)}};

const Matrix Z = {{Complex(1, 0), Complex(0, 0)},
                  {Complex(0, 0), Complex(-1, 0)}};

const Matrix S = {{Complex(1, 0), Complex(0, 0)},
                  {Complex(0, 0), Complex(0, 1)}};

// Add still: T, Y, ...

// Projection Operations
const Matrix Proj_0 = {{Complex(1, 0), Complex(0, 0)},
                       {Complex(0, 0), Complex(0, 0)}};
const Matrix Proj_1 = {{Complex(0, 0), Complex(0, 0)},
                       {Complex(0, 0), Complex(1, 0)}};

// Helper Functions
void printStateVector(StateVector SV) {
  cout << "[\n";
  for (auto c : SV) {
    cout << "  " << c.toString() << "\n";
  }
  cout << "]\n";
}

void printMatrix(Matrix M) {
  cout << "[\n";
  for (auto r : M) {
    for (auto c : r) {
      cout << "  " << c.toString();
    }
    cout << "  \n";
  }
  cout << "]\n";
}

Matrix initSquareMatrix(size_t size) {
  return Matrix(size, std::vector<Complex>(size, Complex(0, 0)));
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

StateVector makeTargetStateVector(string state) {
  StateVector SV;
  for (int i = 0; i < state.size(); i++) {
    SV.push_back(Complex(state[i] - '0', 0));
  }
  return SV;
}

// Test Functions
int testMatrixEqual(Matrix M, Matrix M_expected, string test_name = "",
                    int print_out = 0) {
  // assuming square matrices (probably a better way to compare the matrixs but
  // it is what it is)
  int passed = 1;
  if (M.size() != M_expected.size() ||
      M.at(0).size() != M_expected.at(0).size()) {
    passed = 0;
  } else {
    for (int r = 0; r < M.size(); r++) {
      for (int c = 0; c < M.size(); c++) {
        if (M[r][c].x != M_expected[r][c].x) {
          passed = 0;
        } else if (M[r][c].y != M_expected[r][c].y) {
          passed = 0;
        }
      }
    }
  }
  if (print_out == 1) {
    printMatrix(M);
  }
  if (test_name != "") {
    cout << test_name << " => ";
  }
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}

int testStateVectorsEqual(StateVector SV, StateVector SV_expected,
                          string test_name = "", int print_out = 0) {
  int passed = 1;
  if (SV.size() != SV_expected.size()) {
    passed = 0;
  } else {
    for (int r = 0; r < SV.size(); r++) {
      if (SV[r].x != SV_expected[r].x) {
        passed = 0;
      } else if (SV[r].y != SV_expected[r].y) {
        passed = 0;
      }
    }
  }
  if (print_out == 1) {
    printStateVector(SV);
  }
  if (test_name != "") {
    cout << test_name << " => ";
  }
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}

// Operations
/* Note:
  These functions (tensorProduct, tensorSeries, matrixVectorMultiply, and
  possibly matrixAddition but possibly not necessary) are the functions we
  should basically push to be processed effic. via CUDA. The idea here was to
  abstract these from the actual circuit running code and thus should be pretty
  straight forward to just parallize in CUDA and plug right back in here
  hopefully. That would remove the bulk of the heavy lifting from the CPU and
  give us plenty of different ways to optimize via CUDA to satisfy the project
  requirements.
*/
Matrix tensorProduct(Matrix M, Matrix N) {
  // M tensor N = P

  // assuming M & N are square matrices; and M ⊗ N;
  int size_M = M.size();
  int size_N = N.size();

  // init output matrix of size (size_M*size_N) x (size_M*size_N)
  // cout << "Size of matrix: " << size_M * size_N << " x " << size_M * size_N
  // << "; \n";
  Matrix P = initSquareMatrix(size_M * size_N);

  // iter matrix M
  for (int col_M_i = 0; col_M_i < size_M; col_M_i++) {
    for (int row_M_i = 0; row_M_i < size_M; row_M_i++) {
      Complex M_val = M[row_M_i][col_M_i];

      // iter matrix N
      for (int col_N_i = 0; col_N_i < size_N; col_N_i++) {
        for (int row_N_i = 0; row_N_i < size_N; row_N_i++) {
          Complex N_val = N[row_N_i][col_N_i];

          // find index in P
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
  for (int i = 1; i < M_Vec.size(); i++) {
    U = tensorProduct(U, M_Vec.at(i));
  }
  return U;
}

Matrix matrixMultiply(Matrix M, Matrix N) {
  // implied : if M.size() = M.at(0).size() & (same for N) & M.size() ==
  // N.size()
  Matrix P = initSquareMatrix(M.size());
  int numRows = M.size();
  int numCols = M.size();

  for (int i = 0; i < numRows; i++) {
    for (int j = 0; j < numCols; j++) {
      for (int k = 0; k < numCols; k++) {
        P[i][j].x += M[i][k].x * N[k][j].x;
        P[i][j].y += M[i][k].y * N[k][j].y;
      }
    }
  }

  return P;
}

StateVector matrixVectorMultiply(StateVector SV, Matrix M) {
  // init output state vector to all zeros
  StateVector O_SV;
  int n_Rows_M = M.size();

  int n_Cols_M = M[0].size();

  for (int row_i = 0; row_i < n_Rows_M; row_i++) {
    Complex prod = Complex(0, 0);
    // iter rows of each col in M
    for (int col_i = 0; col_i < n_Cols_M; col_i++) {
      // sparse matrix, add check for 0 in M[row_i][col_i], skip next for if 0
      // (fine here) iter rows of SV vetor
      prod.x += M[row_i][col_i].x * SV[col_i].x;
      prod.y += M[row_i][col_i].y * SV[col_i].y;
    }
    O_SV.push_back(prod);
  }
  return O_SV;
}

Matrix matrixAddition(Matrix M, Matrix N) {
  // ensure that matrices are the same size & square (they should always be
  // square in our case but still doesnt hurt to check)
  // TODO
  Matrix P = initSquareMatrix(M.size());
  for (int r = 0; r < M.size(); r++) {
    for (int c = 0; c < M.size(); c++) {
      P[r][c].x = M[r][c].x + N[r][c].x;
      P[r][c].y = M[r][c].y + N[r][c].y;
    }
  }
  return P;
}

/* Note:
  the process here allows us to support any configuration of a single
  control/target controlled matrix. it has to be done in this goofy ass/mildly
  hardcoded way.

  it's somewhat of an odd process to create arbitrary controlled operations that
  scales exponentially with the # of control qubits you want to encode for a
  given slice, so for right now (until I can figure out if it's feasible to
  offer inherent support for something like the Toffoli operation) any
  multi-control gate will have to be decomposed into a set of operations that
  are supported.


*/
Matrix collapseControlledMatrixVector(vector<Matrix> M_Vec, Matrix U,
                                      int control, int target) {
  vector<Matrix> TS_Ctrl(M_Vec);
  vector<Matrix> TS_Targ(M_Vec);
  // this essentially encodes two scenarios:
  // 1.) when the control qubit is 0, the target function should be the identity
  // (i.e. leave it alone) 2.) when the control qubit is 1, the target function
  // should be the given matrix U
  // these two "scenarios" have to then be tensored and added together to
  // correctly encode any 2 qubit controlled operation
  TS_Ctrl.at(control) = Proj_0;
  TS_Targ.at(control) = Proj_1;
  TS_Targ.at(target) = U;
  return matrixAddition(tensorSeries(TS_Ctrl), tensorSeries(TS_Targ));
}

vector<vector<int>> cartesian_product_01(int n) {
  int num_combinations = pow(2, n);
  vector<vector<int>> carts(num_combinations, vector<int>(n));
  for (int i = 0; i < num_combinations; ++i) {
    for (int j = n - 1; j >= 0; --j) {
      carts.at(i).at(j) = ((i >> j) & 1);
      cout << ((i >> j) & 1) << " ";
    }
    cout << endl;
  }
  return carts;
}

Matrix collapseMCMT(vector<Matrix> M_Vec, Matrix U, vector<size_t> control,
                    vector<size_t> target) {
  // here to handle if the ctrl and targs are just 1
  if (control.size() == 1 && target.size() == 1) {
    return collapseControlledMatrixVector(M_Vec, U, control.at(0),
    target.at(0));
  }

  vector<vector<int>> carts = cartesian_product_01(control.size());

  //# of tensor vectors is 2^(# of controls)
  vector<vector<Matrix>> tensors(carts.size(), M_Vec);

  for (int t = 0; t < tensors.size(); t++) {
    // preset all targets to the unitary
    for (int tar : target) {
      tensors.at(t).at(tar) = U;
    }
    // handle setting controls
    for (int i = 0; i < control.size(); i++) {
      int ctrl_num = control.at(i);
      // tensors.at(t).at(ctrl_num) = ;
      if (carts.at(t).at(i) == 0) {
        tensors.at(t).at(ctrl_num) = Proj_0;
        for (int tar : target) {
          tensors.at(t).at(tar) = I;
        }
      } else {
        tensors.at(t).at(control.at(i)) = Proj_1;
      }
    }
  }
  // collapse tensors to single tensor
  Matrix tot(pow(2, M_Vec.size()),
             vector<Complex>(pow(2, M_Vec.size()), {0, 0}));
  for (auto sub_Vec : tensors) {
    auto t = tensorSeries(sub_Vec);
    tot = matrixAddition(tot, t);
  }
  return tot;
}

class Gate {
public:
  virtual Matrix toMatrix() const = 0;
  virtual string toGateString() const = 0;
  virtual string toString() const = 0;
};

class OneQubitGate : public Gate {
public:
  size_t wireIdx;

  virtual string toString() const override {
    ostringstream oss;
    oss << "(" << toGateString() << " " << wireIdx << ")";
    return oss.str();
  }

protected:
  OneQubitGate(){};
  OneQubitGate(size_t w) : wireIdx(w){};
};

class ControlledGate : public Gate {
public:
  vector<size_t> controlWireIdx;
  vector<size_t> targetWireIdx;
  //TODO fix the to string to handle multiple
  virtual string toString() const override {
    ostringstream oss;
    oss << "(" << toGateString() << " " << controlWireIdx.at(0) << " "
        << targetWireIdx.at(0) << ")";
    return oss.str();
  }

protected:
  ControlledGate(){};
  ControlledGate(size_t c, size_t t) {
    controlWireIdx = vector<size_t>(c);
    targetWireIdx = vector<size_t>(t);
  };
  ControlledGate(vector<size_t> c, vector<size_t> t) : controlWireIdx(c), targetWireIdx(t){};
};

class H_Gate : public OneQubitGate {
public:
  H_Gate(){};
  H_Gate(size_t w) : OneQubitGate(w){};

  virtual string toGateString() const override { return "H"; };
  virtual Matrix toMatrix() const override {
    return {{Complex(INV_SQRT_2, 0), Complex(INV_SQRT_2, 0)},
            {Complex(INV_SQRT_2, 0), Complex(-1 * INV_SQRT_2, 0)}};
  };
};

class X_Gate : public OneQubitGate {
public:
  X_Gate(){};
  X_Gate(size_t w) : OneQubitGate(w){};

  virtual string toGateString() const override { return "X"; };
  virtual Matrix toMatrix() const override {
    return {{{0, 0}, {1, 0}}, //
            {{1, 0}, {0, 0}}};
  };
};

class Z_Gate : public OneQubitGate {
public:
  Z_Gate(){};
  Z_Gate(size_t w) : OneQubitGate(w){};

  virtual string toGateString() const override { return "Z"; };
  virtual Matrix toMatrix() const override {
    return {{{1, 0}, {0, 0}}, {{0, 0}, {-1, 0}}};
  };
};

// Single Qubit Gates to add: Y, S (phase), T (transversal)

class CX_Gate : public ControlledGate {
public:
  CX_Gate(){};
  CX_Gate(size_t c, size_t t) : ControlledGate(c, t){};
  virtual string toGateString() const override { return "CX"; };
  virtual Matrix toMatrix() const override {
    return {{{0, 0}, {1, 0}}, //
            {{1, 0}, {0, 0}}};
  };
};

class CZ_Gate : public ControlledGate {
public:
  CZ_Gate(){};
  CZ_Gate(size_t c, size_t t) : ControlledGate(c, t){};
  virtual string toGateString() const override { return "CZ"; };
  virtual Matrix toMatrix() const override {
    return {{{1, 0}, {0, 0}}, {{0, 0}, {-1, 0}}};
  };
};

class TimeSlice {
  //  virtual TimeSlice() = 0;
public:
  virtual string toString() const = 0;
};

class PeekTimeSlice : public TimeSlice {
public:
  PeekTimeSlice(){};
  virtual string toString() const override { return "(peek)"; }
};

class OpTimeSlice : public TimeSlice {
public:
  virtual Matrix toTransformation() const = 0;
};

class CompiledTimeSlice : public OpTimeSlice {
public:
  Matrix theMatrix;
  //  CompiledTimeSlice(Matrix m) : theMatrix(m){};
  CompiledTimeSlice(vector<OpTimeSlice *> TS) {
    // Calculating matrix at timeslice from R to L
    // CompiledTimeslice C_TS;
    // vector<Matrix> TS_Vec;
    components = TS.size();

    Matrix M = TS.at(TS.size() - 1)->toTransformation();
    for (int i = (TS.size() - 2); i >= 0; i--) {
      M = matrixMultiply(M, TS.at(i)->toTransformation());
    }
    theMatrix = M;
  }

private:
  size_t components;

  virtual string toString() const override {
    return "(compiled[" + to_string(components) + "])";
  }
  virtual Matrix toTransformation() const override { return theMatrix; }
};

class GateTimeSlice : public OpTimeSlice {
public:
  vector<Gate *> gates;
  size_t nQubits;

  GateTimeSlice(vector<Gate *> _gates, size_t _nQubits)
      : gates(_gates), nQubits(_nQubits){};

  GateTimeSlice(size_t nQ) : nQubits(nQ){};

  void addGate(Gate *g) { gates.push_back(g); }

  /* currently tinkering with a better/more effic. method for providing the
unitary here that will wrap the processGates function, hence the arbitrary
function wrap here. */

  // returns the combined unitary for this timeslice

  virtual Matrix toTransformation() const override {
    return toTransformationAux(/*gates, nQubits*/);
  }

  virtual string toString() const override {
    ostringstream oss;
    for (auto gate : gates) {
      oss << gate->toString() << " ";
    }
    return oss.str();
  }

private:
  Matrix
  toTransformationAux(/*vector<Gate *> _gates, int num_of_wires*/) const {
    vector<Matrix> M_Vec(nQubits, I);
    // currently testing an implementation allowing for more than one controlled
    // operation per slice
    vector<ControlledGate *> control_queue;
    for (auto g : this->gates) {
      // if gate is control (not best way to do this but it works for now)
      if (auto *oneQG = dynamic_cast<OneQubitGate *>(g)) {
        M_Vec.at(oneQG->wireIdx) = oneQG->toMatrix();
      } else if (auto *twoQG = dynamic_cast<ControlledGate *>(g)) {
        control_queue.push_back(twoQG);
      } else {
        cout << "unknown gate type in toTransformation, somehow?\n";
        exit(1);
      }
    }
    if (!control_queue.empty()) {
      // again, leaving this as a "queue"/vector for the above reasons even
      // though it currently will only pull the first
      ControlledGate *g = control_queue.at(0);
      return collapseMCMT(
            M_Vec, g->toMatrix(), g->controlWireIdx, g->targetWireIdx);
    } else {
      return tensorSeries(M_Vec);
    }
  }
};

class Circuit {
public:
  size_t nQubits;
  vector<TimeSlice *> program;

  Circuit(size_t n) : nQubits(n){};
  Circuit(vector<TimeSlice *> s, size_t n) : nQubits(n), program(s){};

  void addTimeSlice(TimeSlice *TS) { program.push_back(TS); }

  StateVector runToPosition(StateVector SV, int sliceIdx) {
    StateVector SV_t(SV);
    for (int i = 0; i < sliceIdx; i++) {
      if (const auto *op = dynamic_cast<OpTimeSlice *>(program.at(i))) {
        SV_t = matrixVectorMultiply(SV_t, op->toTransformation());
      }
    }
    return SV_t;
  }

  // Computes the final state vector of the circuit
  StateVector run(StateVector SV) {
    return runToPosition(SV, program.size());
  }

  void compile() {
    vector<TimeSlice *> newProgram;
    vector<OpTimeSlice *> acc;
    for (auto slice : program) {
      if (auto *op = dynamic_cast<OpTimeSlice *>(slice)) {
        acc.push_back(op);
      } else {
        newProgram.push_back(new CompiledTimeSlice(acc));
        acc = {};
        newProgram.push_back(slice);
      }
    }
    newProgram.push_back(new CompiledTimeSlice(acc));
    program = newProgram;
  }

  // Helpers //
  void print() {
    cout << "Circuit with " << nQubits << " cubits and " << program.size()
         << " timesteps:\n";
    for (const auto &ts : program) {
      cout << ts->toString() << "\n";
    }
    return;
  }
};

/// parsing

optional<Gate *> tryParseOneQubitGate(char c) {
  switch (c) {
  case 'H':
    return new H_Gate();
  case 'X':
    return new X_Gate();
  default:
    return {};
  }
}

optional<Gate *> tryParseControlledGate(char c) {
  switch (c) {
  // TODO this is wrong because uncontrolled z exists
  case 'Z':
    return new CX_Gate();
  default:
    return {};
  }
}

template <typename T> void printVec(vector<T> v) {
  cout << "<";
  for (auto e : v) {
    cout << e << " ";
  }
  cout << ">";
}

bool isNonGateCircuitChar(char c) { return c == '.' || c == '-'; }

size_t findControlMark(size_t gateIdx, char gateChar, vector<char> slice) {
  for (size_t idx = 0; idx < slice.size(); idx++) {
    if (gateChar == 'x' && idx != gateIdx) {
      return idx;
    } else if (slice[idx] == '.') {
      return idx;
    }
  }
  cout << "no control mark in slice: ";
  printVec(slice);
  cout << "\n";
  exit(1);
}

char trySingleCharOfDiagramSlice(vector<char> slice) {
  char same = slice[0];
  for (auto c : slice) {
    if (c != same) {
      return '\0';
    }
  }
  return same;
}

bool isNonSemanticDiagramSlice(vector<char> slice) {
  char same = trySingleCharOfDiagramSlice(slice);
  return (same == '|' || same == '0' || same == '>' || same == '-');
}

// a simple diagram slice parses to a single time slice
optional<TimeSlice *> tryParseSimpleDiagramSlice(vector<char> slice) {

  if (trySingleCharOfDiagramSlice(slice) == '!') {
    return new PeekTimeSlice();
  }

  GateTimeSlice *TS = new GateTimeSlice(slice.size());

  for (size_t idx = 0; idx < slice.size(); idx++) {
    if (auto oQG = tryParseOneQubitGate(slice[idx])) {
      auto gate = static_cast<OneQubitGate *>(oQG.value());
      gate->wireIdx = idx;
      TS->gates.push_back(gate);
    } else if (auto tQG = tryParseControlledGate(slice[idx])) {
      auto gate = static_cast<ControlledGate *>(tQG.value());
      gate->targetWireIdx = {idx};
      gate->controlWireIdx = {findControlMark(idx, slice[idx], slice)};
      TS->gates.push_back(gate);
    }
  }

  if (TS->gates.empty()) {
    return {};
  }
  return TS;
}

optional<vector<TimeSlice *>> tryParseSwapGate(vector<char> slice) {

  vector<size_t> swapQubits;
  bool seenOtherGate;
  for (size_t idx = 0; idx < slice.size(); idx++) {
    if (slice[idx] == 'x') {
      swapQubits.push_back(idx);
    } else if (slice[idx] != '-') {
      seenOtherGate = true;
    }
  }
  if (swapQubits.empty()) {
    return {};
  }
  if (seenOtherGate) {
    cout << "swap gate with something else in same slice: ";
    printVec(slice);
    cout << "\n";
    exit(1);
  }
  if (swapQubits.size() != 2) {
    cout << "too many swap gate markers: ";
    printVec(slice);
    cout << "\n";
    exit(1);
  }

  size_t A = swapQubits[0];
  size_t B = swapQubits[1];
  auto G = new CX_Gate(A, B);
  cout << "Swap Qubits A: " << A << "\n";
  cout << "Swap Qubits B: " << B << "\n";
  cout << "Swap Gate Control Wires Size: " << G->controlWireIdx.size() <<"\n";
  cout << "Swap Gate Target Wires Size: " << G->targetWireIdx.size() << "\n";
  return {{
      new GateTimeSlice({new CX_Gate(A, B)}, slice.size()),
      new GateTimeSlice({new CX_Gate(B, A)}, slice.size()),
      new GateTimeSlice({new CX_Gate(A, B)}, slice.size()),
  }};
}

// a complex diagram slice compiles to multiple time slices
optional<vector<TimeSlice *>> tryParseComplexDiagramSlice(vector<char> slice) {
  // special casing is probably fine? use this when we have more complex gates
  /*
  if (auto swap = tryParseSwapGate(slice)) {
    return swap.value();
  }
  return {};
  */
  return tryParseSwapGate(slice);
}

Circuit parseCircuitDiagram(string D) {
  vector<string> wireStrs;
  auto str = stringstream(D);
  for (string line; getline(str, line, '\n');) {
    wireStrs.push_back(line);
  }

  vector<vector<char>> slices;
  for (size_t hIdx = 0; hIdx < wireStrs[0].length(); hIdx++) {
    slices.push_back(vector<char>());
    for (auto wireStr : wireStrs) {
      slices[hIdx].push_back(wireStr[hIdx]);
    }
  }

  Circuit circuit(wireStrs.size());

  for (auto slice : slices) {
    if (isNonSemanticDiagramSlice(slice)) {
      ;
    } else if (auto cts = tryParseComplexDiagramSlice(slice)) {
      for (auto ts : cts.value()) {
        circuit.program.push_back(ts);
      }
    } else if (auto ts = tryParseSimpleDiagramSlice(slice)) {
      circuit.program.push_back(ts.value());
    } else {
      cout << "bad diagram slice: ";
      printVec(slice);
      cout << "\n";
      exit(1);
    }
  }
  // cout << "end of parse";
  return circuit;
}

string gateToString(Gate *gate) {
  std::ostringstream oss;

  if (const auto *oneQG = dynamic_cast<OneQubitGate *>(gate)) {
    oss << "(" << oneQG->toGateString() << " " << to_string(oneQG->wireIdx)
        << ")\n";
  } else if (const auto *twoQG = dynamic_cast<ControlledGate *>(gate)) {
    //TODO - fix this to account for MCMT gates
    oss << "(" << twoQG->toGateString() << " " << twoQG->controlWireIdx.at(0) << " "
        << twoQG->targetWireIdx.at(0) << ")\n";
        cout << "never gets to this point here \n";
  } else {
    oss << "(unknown)\n";
  }
  return oss.str();
}

int main(void) {
  // auto circuit = parseCircuitDiagram("|0>--.\n"
  //                                    "|0>--X\n");
  auto circuit = parseCircuitDiagram("|0>-H-.!---x\n"
                                     "|0>-H-Z!-H-x\n");

  circuit.print();
  cout << "adter print \n";
  printStateVector(circuit.run(makeStateVector(2)));

  circuit.compile();

  circuit.print();

  printStateVector(circuit.run(makeStateVector(2)));
  return 0;
}