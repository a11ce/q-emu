#include "cpu.cpp"

int main(void) {
  /*
  - Xparse swap into 3 cnots
  - Xparse bangs into peek timeslice
  - compile slices
    -  split compilation at peek times


  */

  // parse test
  auto circuit = parseCircuitDiagram("|0>-H-.!---x!\n"
                                     "|0>-H-Z!-H-x!\n");

  circuit.print();
  printStateVector(circuit.run(makeStateVector(2)));
  cerr << "compiling in main\n";
  circuit.compile();
  cerr << "compiled in main\n";

  circuit.print();
  return 0;

  printStateVector(circuit.run(makeStateVector(2)));

  auto c2 = parseCircuitDiagram("|0>-H\n"
                                "|0>--");
  c2.print();
  printStateVector(c2.run(makeStateVector(2)));

  auto c3 = parseCircuitDiagram("|0>-H-x\n"
                                "|0>---x");
  c3.print();
  printStateVector(c3.run(makeStateVector(2)));

  // Slices & Tensoring //

  GateTimeSlice TS_0({}, 3);
  Matrix TestIII = {
      {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)}};
  testMatrixEqual(TS_0.toTransformation(), TestIII, "I ⊗ I ⊗ I");

  GateTimeSlice TS_1({new X_Gate(2)}, 3);
  Matrix TestIIX = {
      {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)}};

  testMatrixEqual(TS_1.toTransformation(), TestIIX, "I ⊗ I ⊗ X");

  GateTimeSlice TS_2({new X_Gate(0), new CX_Gate(2, 1)}, 3);
  Matrix TestXICX = {
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(1, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(1, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(0, 0), Complex(1, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)},
      {Complex(0, 0), Complex(1, 0), Complex(0, 0), Complex(0, 0),
       Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)}};
  testMatrixEqual(TS_2.toTransformation(), TestXICX, "X ⊗ CX(2, 1)");

  // Baby Circuits //
  Circuit Circ_1({new GateTimeSlice({new X_Gate(0)}, 2),  //
                  new GateTimeSlice({new X_Gate(1)}, 2)}, //
                 2);
  StateVector Ket11 = {Complex(0, 0), Complex(0, 0), Complex(0, 0),
                       Complex(1, 0)};
  testStateVectorsEqual(Circ_1.run(makeStateVector(2)), Ket11,
                        "(I ⊗ X((X ⊗ I)|00>))");

  GateTimeSlice *Circ_2_TS_0 = new GateTimeSlice({new X_Gate(1)}, 2);
  GateTimeSlice *Circ_2_TS_1 = new GateTimeSlice({new CX_Gate(0, 1)}, 2);
  GateTimeSlice *Circ_2_TS_2 = new GateTimeSlice({new CX_Gate(1, 0)}, 2);

  Circuit Circ_2({Circ_2_TS_0, Circ_2_TS_1, Circ_2_TS_2}, 2);

  StateVector Ket01 = {Complex(0, 0), Complex(1, 0), Complex(0, 0),
                       Complex(0, 0)};
  testStateVectorsEqual(Circ_2.runToPosition(makeStateVector(2), 1), Ket01,
                        "(I ⊗ X)|00>))");
  testStateVectorsEqual(Circ_2.runToPosition(makeStateVector(2), 2), Ket01,
                        "(CX(0,1))|01>))");
  testStateVectorsEqual(Circ_2.runToPosition(makeStateVector(2), 3), Ket11,
                        "(CX(1,0)|01>))");

  // test 2
  //  auto circuit_2 = parseCircuitDiagram("|0>-Z-.!---X\n"
  //                                     "|0>-X-Z!-X-Z");

  GateTimeSlice *Circ_3_TS_0 =
      new GateTimeSlice({new Z_Gate(0), new X_Gate(1)}, 2);
  GateTimeSlice *Circ_3_TS_1 = new GateTimeSlice({new CZ_Gate(0, 1)}, 2);

  GateTimeSlice *Circ_3_TS_2 = new GateTimeSlice({new X_Gate(1)}, 2);
  GateTimeSlice *Circ_3_TS_3 =
      new GateTimeSlice({new X_Gate(0), new Z_Gate(1)}, 2);

  // printMatrix(Circ_3_TS_0->toTransformation());
  // printMatrix(Circ_3_TS_1->toTransformation());

  Matrix C_TS_0 = matrixMultiply(Circ_3_TS_1->toTransformation(),
                                 Circ_3_TS_0->toTransformation());
  Matrix C_TS_1 = matrixMultiply(Circ_3_TS_3->toTransformation(),
                                 Circ_3_TS_2->toTransformation());
  Circuit Circ_3_Total({Circ_3_TS_0, Circ_3_TS_1, Circ_3_TS_2, Circ_3_TS_3}, 2);

  // printMatrix(C_TS_0);
  // printMatrix(C_TS_1);
  // printMatrix(matrixMultiply(C_TS_1, C_TS_0));
  testStateVectorsEqual(Circ_3_Total.runToPosition(makeStateVector(2), 1),
                        matrixVectorMultiply(makeStateVector(2), C_TS_0),
                        "(CZ(0,1)(Z ⊗ X(|00>)))");

  vector<Matrix> M_Vec = {I, I, I};
  vector<size_t> control = {0};
  vector<size_t> target = {1};

  auto tensors = collapseMCMT(M_Vec, Z, control, target);
  printMatrix(tensors);
  // for (auto s : tensors) {
  //   for (auto m : s) {
  //     printMatrix(m);
  //   }
  //   cout << "\n --------- \n";
  // }

  /// Grovers
  auto SV = makeTargetStateVector("100");
  testStateVectorsEqual(SV, {Complex(1, 0), Complex(0, 0), Complex(0, 0)},
                        "Test Target State Vector", 1);

  printStateVector(circuit.run(SV));
}
