#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


// class Circuit {
//   size_t nCubits;
//   vector<TimeSlice> program;
// }

// class TimeSlice {
//   vector<Gate> gates;
//   Matrix unitary() {
//     return tensorSeries(this.gates);
//   }; 
// }
/*
class Gate {

}

class OneQubitGate {
  size_t wireIdx;
  Matrix toMatrix() virtual;

}

class TwoQubitGate {
  size_t controlWireIdx;
  size_t targetWireIdx;
  Matrix toMatrix() virtual;

}

class H OneQubitGate {
  Matrix toMatrix() override;
}

class CNOT TwoQubitGate {
  Matrix toMatrix() override;
}
*/
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

using Timeslice = std::vector<Matrix>;

// Helper Functions
//modified this as the size of the state vector is 2^# of qubits
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

//assuming square matrices (probably a better way to compare the matrixs but it is what it is)
int testMatrixEqual(Matrix M, Matrix M_expected, int print_out = 0){
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
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}

int testStateVectorsEqual(StateVector SV, StateVector SV_expected){
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
  cout << "Passed?: " << (passed == 1 ? "TRUE" : "FALSE") << "\n";
  return passed;
}
//idk if this is totally necessary, but it made me sad to type once so i just threw it in a function
//to avoid having to type it again down the line
Matrix initSquareMatrix(size_t size) {
  return Matrix(size, std::vector<Complex>(size, Complex(0,0)));
}

Matrix projectionMatrixOf(int Ket) {
  //this is absolute dogshit rn but for testing,
  // what we want here is a special matrix for the outer product of |0> (i.e. |0><0|) or |1> (i.e. |1><1|)
  if (Ket == 0) {
    return {{Complex(1, 0), Complex(0, 0)},
            {Complex(0, 0), Complex(0, 0)}};
  } else if (Ket == 1) {
    return {{Complex(0, 0), Complex(0, 0)},
            {Complex(0, 0), Complex(1, 0)}};
  // This also is not right but its fine for now
  } else {
    return {{Complex(1, 0), Complex(0, 0)},
            {Complex(0, 0), Complex(1, 0)}};
  }
}

/////
// Transformation/Operation Functions
StateVector applyTransformation(StateVector SV, Matrix M) { 
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

Matrix tensorSeries(Timeslice TS) {
  //the TimeSlice represents a given timesplice in which elements are the operations are ordered in the vector from top wire to bottom wire
  int n = TS.size();
  Matrix U;
  //again this is dogshit handle this better but works for now
  if(n < 2) {
    U = TS.at(0);
    return U;
  } else {
    U = TS.at(0);
    // printMatrix(U);
    //set the unitary to the first element of the array
    for(int i = 1; i < TS.size(); i++) {
      //check that matrices are compatible size
      //TODO

      //multiply these two matrices
      U = tensorProduct(U, TS.at(i));
      // printMatrix(U);
    }
    return U;
  }
}

Matrix controlledOperation(Matrix u, int control_index, int target_index, int num_of_wires){
  //this needs to be update but since we usually only perform these multi-qubit controlled operations in their own timeslice, it might be okay?
  Matrix I = {{Complex(1, 0), Complex(0, 0)},
              {Complex(0, 0), Complex(1, 0)}};

  Timeslice TS_Ctrl = Timeslice(num_of_wires, I);
  Timeslice TS_Targ = Timeslice(num_of_wires, I);

  // Timeslice TS_Ctrl;
  // Timeslice TS_Targ;

  // TS_Ctrl.gates = vector<Matrix>(num_of_wires, I);
  // TS_Targ.gates = vector<Matrix>(num_of_wires, I);

  TS_Ctrl.at(control_index) = projectionMatrixOf(0);
  TS_Targ.at(control_index) = projectionMatrixOf(1);
  TS_Targ.at(target_index) = u;

  Matrix TS_Contr = matrixAddition(tensorSeries(TS_Ctrl), tensorSeries(TS_Targ));

  return TS_Contr;
}

int main(void) {
  //Gate Operations
  Matrix I = {{Complex(1, 0), Complex(0, 0)},
              {Complex(0, 0), Complex(1, 0)}};
  
  Matrix H = {{Complex(1, 0), Complex(1, 0)},
              {Complex(1, 0), Complex(-1, 0)}};

  Matrix X = {{Complex(0, 0), Complex(1, 0)},
              {Complex(1, 0), Complex(0, 0)}};
  
  
  //Single-Qubit Operation Test
  
  cout << "1-Qubit Transformation: \n";
  auto SV = makeStateVector(1);
  cout << "Input State Vector |psi0>: \n";
  printStateVector(SV);

  cout << "Unitary Transformation, U: \n";

  cout << "Output State Vector, |psi1>: \n";
  StateVector O_SV = applyTransformation(SV, H);

  // printStateVector(O_SV);
  StateVector O_SV_expected = {Complex(1,0), Complex(1,0)};
  testStateVectorsEqual(O_SV, O_SV_expected);

  // I tensor I Test
  cout << "\nTesting Tensor Product of I ⊗ I: \n";
  cout << "Identity Matrix: \n";

  auto I_Tensor_I = tensorProduct(I, I);

  cout << "I ⊗ I: \n";
  // printMatrix(I_Tensor_I);

  Matrix expected_II = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
                        
  testMatrixEqual(I_Tensor_I, expected_II);

  // H tensor I test
  cout << "\nTesting Tensor Product of H ⊗ I: \n";
  cout << "Hadamard Matrix: \n";

  auto H_Tensor_I = tensorProduct(H, I);
  cout << "H ⊗ I: \n";

  Matrix expected_HI = {{Complex(1, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(1, 0)},
                        {Complex(1, 0), Complex(0, 0), Complex(-1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(-1, 0)}};
  testMatrixEqual(H_Tensor_I, expected_HI);

  // 2 Qubit Operation Test
  cout << "\n2-Qubit Transformation: \n";
  auto SV_2 = makeStateVector(2);
  cout << "Input State Vector |psi> = |0 0>: \n";
  // printStateVector(SV_2);

  cout << "Unitary Transformation, U = H ⊗ I \n";
  
  auto O_SV_2 = applyTransformation(SV_2, tensorProduct(H, I));
  StateVector expected_O_SV_2 = {Complex(1, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)};

  cout << "(H ⊗ I)|00>: \n";
  // printStateVector(O_SV_2);

  testStateVectorsEqual(O_SV_2, expected_O_SV_2);

  // state1 = {{Complex(1, 0)}, {Complex(0, 0)}}; // [1
  //                                            //  0]
  // state2 = {{Complex(1, 0), Complex(0, 0)}}; // [1 0]

  // printMatrix(controlMatrixForKet(1));
  
  // Testing tensoring matrices of different sizes
  
  cout << "\nIn-Series Tensoring: \n";
  cout << "I ⊗ I ⊗ I: \n";

  
  // Matrix I ⊗ I ⊗ I
  Matrix ItensorI = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
  Matrix TestIII = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
  Matrix ItensorItensorI = tensorProduct(ItensorI, I);
  testMatrixEqual(ItensorItensorI, TestIII);
  
  // Timeslice I ⊗ I ⊗ X
  cout << "I ⊗ I ⊗ X: \n";
  Timeslice TS1 = {I, I, X};

  // Matrix XII = tensorSeries(TS1);

  Matrix TestIIX = {{Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                     {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};
  testMatrixEqual(tensorSeries(TS1), TestIIX);

  
  cout << "\nControlled Operation Testing\n";
  cout << "CX(0,1) w/ 2 qubits:\n";

  //what this says is:
  //  if wire 0 is 0, leave wire 1 alone
  //  if wire 0 is 1, apply the X operator to wire 1

  // CX (cont. wire = 0; target = 1) =>
  // (|1><1| ⊗ I) ⊗ X =>
  // (|0><0| ⊗ I) + (|1><1| ⊗ X)

  auto CX_0_1 = controlledOperation(X, 0, 1, 2);
  Matrix CX_0_1_expec = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};
  testMatrixEqual(CX_0_1, CX_0_1_expec);

  cout << "CX(1,2) w/ 3 qubits:\n";
  auto CX_1_2 = controlledOperation(X, 1, 2, 3);
  Matrix CX_1_2_expec = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};
  testMatrixEqual(CX_1_2, CX_1_2_expec);

  cout << "CX(2,1) w/ 3 qubits:\n";
  auto CX_2_1 = controlledOperation(X, 2, 1, 3);
  Matrix CX_2_1_expec = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                         {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)}};
  testMatrixEqual(CX_2_1, CX_2_1_expec);


  // printMatrix(matrixAddition(tensorProduct(controlMatrixForKet(0), I), tensorProduct(controlMatrixForKet(1), X)));
  // Matrix CX_0 = tensorSeries({controlMatrixForKet(0), I});
  // Matrix CX_1 = tensorSeries({controlMatrixForKet(1), X});
  // printMatrix(matrixAddition(CX_0, CX_1));

  // this but if the control and target wires are flipped
  //    if wire 1 is 0, leave wire 0 alone
  //    if wire 1 is 1, apply the X operator to wire 1

  // Matrix R_CX_0 = tensorSeries({I, controlMatrixForKet(0)});
  // Matrix R_CX_1 = tensorSeries({X, controlMatrixForKet(1)});
  // printMatrix(matrixAddition(R_CX_0, R_CX_1));

  // printMatrix(controlledOperation(X, 1, 2, 3));
  // printMatrix(controlledOperation(X, 0, 1, 3));
  // printMatrix(controlledOperation(X, 1, 0, 3));
}
