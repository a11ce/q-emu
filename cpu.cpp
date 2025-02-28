#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/*
class Circuit {
  size_t nCubits;
  vector<TimeSlice> program;
}

class TimeSlice {
  vector<Gate> gates;
  Matrix toTransformation(); // tensors them together
}

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

using Matrix = std::vector<std::vector<Complex>>;

using StateVector = vector<Complex>;

//idk if this is totally necessary, but it made me sad to type once so i just threw it in a function
//to avoid having to type it again down the line
Matrix initSquareMatrix(size_t size) {
  return Matrix(size, std::vector<Complex>(size, Complex(0,0)));
}

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

Matrix tensorProduct(Matrix M, Matrix N){
  //M tensor N = P

  //assuming M & N are square matrices
  int size_M = M.size();
  int size_N = N.size();

  //this is just handling tensoring two matrices of the same size, should add some constraints and adjust for arbitray valid sized matrices
  //init output matrix of size (size_M*size_N) x (size_M*size_N)
  Matrix P = initSquareMatrix(size_M * size_M);
  
  //iter matrix M
  for(int col_M_i = 0; col_M_i < size_M; col_M_i++) {
    for(int row_M_i = 0; row_M_i < size_N; row_M_i++){
      Complex M_val = M[row_M_i][col_M_i];

      //iter matrix N
      for(int col_N_i = 0; col_N_i < size_M; col_N_i++) {
        for(int row_N_i = 0; row_N_i < size_N; row_N_i++) {
          Complex N_val = N[row_N_i][col_N_i];

          //find index in P
          int P_row = (size_M * row_M_i) + row_N_i;
          int P_col = (size_M * col_M_i) + col_N_i;
          P[P_row][P_col] = Complex(M_val.x * N_val.x, M_val.y * N_val.y);
        }
      }
    }
  }

  return P;
}

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
int testMatrixEqual(Matrix M, Matrix M_expected){
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

int main(void) {
  //Single-Qubit Operation Test
  cout << "1-Qubit Transformation: \n";
  auto SV = makeStateVector(1);
  cout << "Input State Vector |psi0>: \n";
  printStateVector(SV);

  Matrix H = {{Complex(1, 0), Complex(1, 0) },
               {Complex(1, 0), Complex(-1, 0) }};

  cout << "Unitary Transformation, U: \n";
  printMatrix(H);
  
  cout << "Output State Vector, |psi1>: \n";
  StateVector O_SV = applyTransformation(SV, H);

  printStateVector(O_SV);
  StateVector O_SV_expected = {Complex(1,0), Complex(1,0)};
  testStateVectorsEqual(O_SV, O_SV_expected);

  // I tensor I Test
  cout << "\nTesting Tensor Product of I ⊗ I: \n";
  Matrix I = {{Complex(1, 0), Complex(0, 0)},
              {Complex(0, 0), Complex(1, 0)}};
  cout << "Identity Matrix: \n";
  printMatrix(I);

  auto I_Tensor_I = tensorProduct(I, I);

  cout << "I ⊗ I: \n";
  printMatrix(I_Tensor_I);

  Matrix expected_II = {{Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
                        
  testMatrixEqual(I_Tensor_I, expected_II);

  // H tensor I test
  cout << "\nTesting Tensor Product of H ⊗ I: \n";
  cout << "Hadamard Matrix: \n";
  printMatrix(H);
  auto H_Tensor_I = tensorProduct(H, I);
  cout << "H ⊗ I: \n";
  printMatrix(H_Tensor_I);

  Matrix expected_HI = {{Complex(1, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(1, 0)},
                        {Complex(1, 0), Complex(0, 0), Complex(-1, 0), Complex(0, 0)},
                        {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(-1, 0)}};
  testMatrixEqual(H_Tensor_I, expected_HI);

  // 2 Qubit Operation Test
  cout << "\n2-Qubit Transformation: \n";
  auto SV_2 = makeStateVector(2);
  cout << "Input State Vector |psi> = |0 0>: \n";
  printStateVector(SV_2);

  cout << "Unitary Transformation, U = H ⊗ I \n";
  
  auto O_SV_2 = applyTransformation(SV_2, tensorProduct(H, I));
  StateVector expected_O_SV_2 = {Complex(1, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)};

  cout << "(H ⊗ I)|00>: \n";
  printStateVector(O_SV_2);

  testStateVectorsEqual(O_SV_2, expected_O_SV_2);
}
