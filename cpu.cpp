#include <iostream>
#include <vector>

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

StateVector applyTransformation(StateVector SV, Matrix M) { return "MEOW"; }

StateVector makeStateVector(size_t nCubits) {
  StateVector SV;
  Complex oneProb = Complex(1, 0);
  Complex zeroProb = Complex(0, 0);
  SV.push_back(oneProb);
  for (size_t idx = 0; idx < nCubits - 1; idx++) {
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

int main(void) {
  auto SV = makeStateVector(2);

  printStateVector(SV);
}
