#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

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

using Matrix = std::vector<std::vector<Complex>>;

using StateVector = vector<Complex>;

class Gate {
public:
  virtual Matrix toMatrix() const = 0;
  virtual char toGateChar() const = 0;
};

class TimeSlice {
public:
  vector<Gate *> gates;
  Matrix toTransformation(); // tensors them together
};

class Circuit {
public:
  size_t nCubits;
  vector<TimeSlice> program;
};

class OneQubitGate : public Gate {
public:
  size_t wireIdx;
};

class TwoQubitGate : public Gate {
public:
  size_t controlWireIdx;
  size_t targetWireIdx;
};

class H_Gate : public OneQubitGate {
public:
  H_Gate(){};
  virtual char toGateChar() const override { return 'H'; };
  virtual Matrix toMatrix() const override { return {{{0, 0}}}; };
};

class CNOT_Gate : public TwoQubitGate {
public:
  CNOT_Gate(){};
  virtual char toGateChar() const override { return 'Z'; };
  virtual Matrix toMatrix() const override { return {{{0, 0}}}; };
};

class SWAP_Gate : public TwoQubitGate {
public:
  SWAP_Gate(){};
  virtual char toGateChar() const override { return 'x'; };
  virtual Matrix toMatrix() const override { return {{{0, 0}}}; };
};

// StateVector applyTransformation(StateVector SV, Matrix M) { return "MEOW"; }

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

optional<Gate *> tryParseOneQubitGate(char c) {
  switch (c) {
  case 'H':
    return new H_Gate();
  default:
    return {};
  }
}

optional<Gate *> tryParseTwoQubitGate(char c) {
  switch (c) {
  case 'Z':
    return new CNOT_Gate();
  case 'x':
    return new SWAP_Gate();
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

optional<TimeSlice> parseTimeSlice(vector<char> slice) {
  char same = slice[0];
  for (auto c : slice) {
    if (c != same) {
      same = -1;
      break;
    }
  }

  if (same == '|' || same == '0' || same == '>' || same == '-') {
    return {};
  }

  TimeSlice TS;

  for (size_t idx = 0; idx < slice.size(); idx++) {
    if (auto oQG = tryParseOneQubitGate(slice[idx])) {
      auto gate = static_cast<OneQubitGate *>(oQG.value());
      gate->wireIdx = idx;
      TS.gates.push_back(gate);
    } else if (auto tQG = tryParseTwoQubitGate(slice[idx])) {
      auto gate = static_cast<TwoQubitGate *>(tQG.value());
      gate->targetWireIdx = idx;
      gate->controlWireIdx = findControlMark(idx, slice[idx], slice);
      // mild hack so we dont get (x 0 1) and (x 1 0)
      if (!(slice[idx] == 'x' && gate->targetWireIdx < gate->controlWireIdx)) {
        TS.gates.push_back(gate);
      }
    } else if (!isNonGateCircuitChar(slice[idx])) {
      cout << "bad char in diagram: " << slice[idx] << "\n";
      exit(1);
    }
  }
  return TS;
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

  Circuit circuit;
  circuit.nCubits = wireStrs.size();

  for (auto slice : slices) {
    auto ts = parseTimeSlice(slice);
    if (ts) {
      circuit.program.push_back(ts.value());
    }
  }

  return circuit;
}

string gateToString(Gate *gate) {
  std::ostringstream oss;

  if (const auto *oneQG = dynamic_cast<OneQubitGate *>(gate)) {
    oss << "(" << oneQG->toGateChar() << " " << to_string(oneQG->wireIdx)
        << ")\n";
  } else if (const auto *twoQG = dynamic_cast<TwoQubitGate *>(gate)) {
    oss << "(" << twoQG->toGateChar() << " " << twoQG->controlWireIdx << " "
        << twoQG->targetWireIdx << ")\n";
  } else {
    oss << "(unknown)\n";
  }
  return oss.str();
}

string timeSliceToString(const TimeSlice &ts) {
  std::ostringstream oss;
  for (auto gate : ts.gates) {
    oss << gateToString(gate);
  }
  return oss.str();
}

void printCircuit(const Circuit &c) {
  cout << "Circuit with " << c.nCubits << " cubits and " << c.program.size()
       << " timesteps:\n";
  for (const auto &ts : c.program) {
    cout << timeSliceToString(ts) << "\n";
  }
  return;
}

int main(void) {

  auto circuit = parseCircuitDiagram("|0>-H-.---x\n"
                                     "|0>-H-Z-H-x");
  //  cout << circuit.program.size();
  printCircuit(circuit);
  return 0;
}
