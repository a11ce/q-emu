#include "cpu.cpp"
#include <cuda/std/complex>

using ResultVector = vector<StateVector>;

__device__ void cuMatrixVectorMultiply(size_t sLen, Complex *inV, Complex *outV,
                                       Complex *M) {
  for (size_t idy = 0; idy < sLen; idy++) {
    outV[idy] = {0, 0};
    // Complex r =  {0,0};
    for (size_t idx = 0; idx < sLen; idx++) {
      outV[idy].x += M[(idy * sLen) + idx].x * inV[idx].x;
      outV[idy].y += M[(idy * sLen) + idx].y * inV[idx].y;
    }
  }
}

__global__ void QProgramKernel(size_t nOps, size_t sLen, Complex *matrices,
                               Complex *inVector, Complex *outVectors) {
  // first is special case because of inVector
  cuMatrixVectorMultiply(sLen, inVector, outVectors, matrices);

  for (size_t idx = 1; idx < nOps; idx++) {
    Complex *matrixBp = matrices + (idx * sLen * sLen);
    Complex *prevVBp = outVectors + ((idx - 1) * sLen);
    Complex *thisVBp = outVectors + (idx * sLen);

    cuMatrixVectorMultiply(sLen, prevVBp, thisVBp, matrixBp);
  }
  return;
}

void lowerMatrix(Matrix m, Complex *matrices, size_t baseOffset) {
  Complex *bp = matrices + baseOffset;
  cerr << "bp is " << bp << "\n";
  for (size_t idy = 0; idy < m.size(); idy++) {
    for (size_t idx = 0; idx < m.size(); idx++) {
      bp[(idy * m.size()) + idx] = m[idy][idx];
    }
  }
}

Complex *lowerStateVector(StateVector SV) {
  Complex *sv = (Complex *)malloc(SV.size() * sizeof(Complex));
  for (size_t idx = 0; idx < SV.size(); idx++) {
    sv[idx] = SV[idx];
  }
  return sv;
}

size_t *lowerPeekPoints(vector<size_t> p) {
  size_t *rp = (size_t *)malloc(p.size() * sizeof(size_t));
  for (size_t idx = 0; idx < p.size(); idx++) {
    rp[idx] = p[idx];
  }
  return rp;
}

void printOutVectors(Complex *OV, size_t sLen, size_t nP) {
  for (size_t idx = 0; idx < nP; idx++) {
    Complex *bp = OV + (sLen * idx);
    cout << "[\n";
    for (size_t idv = 0; idv < sLen; idv++) {
      cout << "  " << bp[idv].toString() << "\n";
    }
    cout << "]\n";
  }
  return;
}

ResultVector runCircuitOnGPU(Circuit C, StateVector SV) {
  C.compile();

  cerr << "compiled\n";
  size_t matrixSideLength = pow(2, C.nQubits);

  size_t sizeOfComplex = sizeof(Complex); // 2 * sizeof(float);
  size_t sizeOfMatrix = sizeOfComplex * matrixSideLength * matrixSideLength;
  size_t sizeOfMatrices = sizeOfMatrix * C.operationsCount();

  Complex *matrices = (Complex *)malloc(sizeOfMatrix * C.operationsCount());

  cerr << "base matrices is " << matrices << "\n";

  size_t pc = 0;

  cerr << "lowering matrices\n";

  for (auto slice : C.program) {
    if (auto *op = dynamic_cast<OpTimeSlice *>(slice)) {
      lowerMatrix(op->toTransformation(), matrices,
                  pc * matrixSideLength * matrixSideLength);
      pc++;
    } else {
      cout << "no-op after compilation\n";
      exit(1);
    }
  }

  cerr << "lowered matrices\n";

  Complex *h_inVector =
      lowerStateVector(SV); // malloc(sizeof(Complex) * matrixSideLength);
                            //    lowerStateVector(h_inVector, SV);
  Complex *d_inVector;
  cudaMalloc((void **)&d_inVector, sizeof(Complex) * matrixSideLength);
  cudaMemcpy(d_inVector, h_inVector, sizeof(Complex) * matrixSideLength,
             cudaMemcpyHostToDevice);

  size_t nOps = C.program.size();

  // matrices now contains the floats
  Complex *d_matrices;
  cudaMalloc((void **)&d_matrices, sizeOfMatrices);
  cudaMemcpy(d_matrices, matrices, sizeOfMatrices, cudaMemcpyHostToDevice);

  size_t totalOutVectorSize = nOps * matrixSideLength * sizeof(Complex);
  Complex *d_outVectors;
  cudaMalloc((void **)&d_outVectors, totalOutVectorSize);
  cerr << "launching\n";

  dim3 dimGrid(1, 1);
  dim3 dimBlock(1, 1);
  QProgramKernel<<<dimGrid, dimBlock>>>(nOps, matrixSideLength, d_matrices,
                                        d_inVector, d_outVectors);
  cerr << "launched\n";

  Complex *h_outVectors = (Complex *)malloc(totalOutVectorSize);

  cudaMemcpy(h_outVectors, d_outVectors, totalOutVectorSize,
             cudaMemcpyDeviceToHost);

  cerr << "memcpyd\n";

  printOutVectors(h_outVectors, matrixSideLength, nOps);
  return {};
}

int main() {
  auto circuit = parseCircuitDiagram("|0>-H!-.!---x\n"
                                     "|0>-H!-Z!-H-x\n");

  printStateVector(circuit.run(makeStateVector(2)));

  auto res = runCircuitOnGPU(circuit, makeStateVector(2));
}
