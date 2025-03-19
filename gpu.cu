#include "cpu.cpp"
#include <cuda/std/complex>

using ResultVector = vector<StateVector>;

__device__ void cuMatrixVectorMultiply(size_t sLen, Complex *inV, Complex *outV,
                                       Complex *M) {
  size_t idy = threadIdx.x;
  if (idy >= sLen) {
    return;
  }
  outV[idy] = {0, 0};
  for (size_t idx = 0; idx < sLen; idx++) {
    outV[idy].x += M[(idy * sLen) + idx].x * inV[idx].x;
    outV[idy].y += M[(idy * sLen) + idx].y * inV[idx].y;
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

__device__ void cuTensorProduct(size_t sLen, Complex *M, Complex *N, Complex *P) {
  int mRow = blockIdx.y * 32 + threadIdx.y;
  int mCol = blockIdx.x * 32 + threadIdx.x;

  // dimensions of P = len(M) * len(N)
  //since len(N) always = 2 & len(M) = sLen --> P = sLen * 2
  int pLen = sLen * 2;

  Complex mVal = M[(mRow * sLen) +mCol];
  ///TODO check that it doesnt overflow the thisVal part
  // if(){

  // }
  for(int col_i = 0; col_i < 2; col_i++) {
    for(int row_i = 0; row_i < 2; row_i++) {
      Complex nVal = N[(row_i * 2) + col_i];
      int P_row = (2 * mRow) + row_i;
      int P_col = (2 * mCol) + col_i;
      P[(P_row * pLen) + P_col] = Complex(mVal.x * nVal.x, mVal.y * nVal.y);
    }
  }
}

__global__ void QTensorKernel(size_t nMatrices, size_t matrixSideLength,
                              Complex *matrices, Complex *outMatrix, Complex *scratchMatrix) {
  // outMatrix already contains matrices[0]
  for (size_t idx = 1; idx < nMatrices; idx++) {
    Complex *matrixBP = matrices + (idx * matrixSideLength * matrixSideLength);
    //Set the matrix side length to = sizeOutMatrix * 2
    // M tensor N = P (where P is in our scratch/temp matrix)
    cuTensorProduct(matrixSideLength, outMatrix, matrixBP, scratchMatrix);
    __syncthreads();
    //copy the scratch/P matrix into the outMatrix and move to the next tensor
    memcpy(outMatrix, scratchMatrix, pow(matrixSideLength*2, 2));
  }
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

void copyRMatrix(Complex *src, Complex *dst, size_t sideLen) {
  memcpy(dst, src, sideLen * sideLen * sizeof(Complex));
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

int ceilDiv(int a, int b)
{
    return ceil((float) a / float(b));
}


Complex *gpu_toTransformationMatrix(
    pair<optional<ControlledGate *>, vector<Matrix>> pti) {
  if (pti.first) {
    // TODO on gpu
    exit(1);
    /*auto g = pti.first.value();
    return collapseMCMT(pti.second, g->toMatrix(), g->controlWireIdx,
                        g->targetWireIdx);*/
  }
  vector<Matrix> series = pti.second;
  size_t matrixSideLength = series[0].size();
  size_t sizeOfMatrix = sizeof(Complex) * matrixSideLength * matrixSideLength;
  size_t nMatrices = series.size();
  size_t sizeOfMatrices = sizeOfMatrix * nMatrices;

  Complex *matrices = (Complex *)malloc(sizeOfMatrix * nMatrices);
  cerr << "lowering for tensor\n";

  size_t pc = 0;
  for (auto mat : series) {
    lowerMatrix(mat, matrices, pc * matrixSideLength * matrixSideLength);
  }

  cerr << "end lowering for tensor\n";

  Complex *d_matrices;
  cudaMalloc((void **)&d_matrices, sizeOfMatrices);
  cudaMemcpy(d_matrices, matrices, sizeOfMatrices, cudaMemcpyHostToDevice);

  size_t outMatrixSideLength = pow(2, series.size());
  size_t totalOutMatrixSize =
      sizeof(Complex) * outMatrixSideLength * outMatrixSideLength;

  Complex *outMatrix = (Complex *)malloc(totalOutMatrixSize);

  Complex *d_outMatrix;
  Complex *d_scratchMatrix;
  cudaMalloc((void **)&d_outMatrix, totalOutMatrixSize);
  cudaMalloc((void **)&d_scratchMatrix, totalOutMatrixSize);
  cudaMemcpy(d_outMatrix, matrices, totalOutMatrixSize, cudaMemcpyHostToDevice);

  // TODO dims
  dim3 dimBlock(32,32);
  dim3 dimGrid(ceilDiv(outMatrixSideLength, 32), ceilDiv(outMatrixSideLength, 32));
  QTensorKernel<<<dimGrid, dimBlock>>>(nMatrices, matrixSideLength, d_matrices,
                                       d_outMatrix, d_scratchMatrix);
  cudaMemcpy(outMatrix, d_outMatrix, totalOutMatrixSize,
             cudaMemcpyDeviceToHost);
  return outMatrix;
}

ResultVector runCircuitOnGPU(Circuit C, StateVector SV) {
  /* cerr << "compiling\n";

   C.compile();

   cerr << "compiled\n";*/
  size_t matrixSideLength = pow(2, C.nQubits);

  size_t sizeOfComplex = sizeof(Complex); // 2 * sizeof(float);
  size_t sizeOfMatrix = sizeOfComplex * matrixSideLength * matrixSideLength;
  size_t sizeOfMatrices = sizeOfMatrix * C.operationsCount();

  Complex *matrices = (Complex *)malloc(sizeOfMatrix * C.operationsCount());

  cerr << "base matrices is " << matrices << "\n";

  size_t pc = 0;

  cerr << "lowering matrices\n";

  for (auto slice : C.program) {
    if (auto *op = dynamic_cast<GateTimeSlice *>(slice)) {
      cerr << "start lowering\n";

      Complex *transformationMatrix =
          gpu_toTransformationMatrix(op->toPreTensorInfo());
      /*
  lowerMatrix(transformationMatrix, op->toTransformation(), matrices,
              pc * matrixSideLength * matrixSideLength);*/
      copyRMatrix(transformationMatrix,
                  matrices + (pc * matrixSideLength * matrixSideLength),
                  matrixSideLength);

      cerr << "end lowering\n";
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

  // dim3 dimGrid(1, 1);
  dim3 dimGrid(1);
  // TODO wrong add multi blocks
  dim3 dimBlock(1024);

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

int main(int argc, char *argv[]) {
  int v = atoi(argv[1]); /*auto circuit = parseCircuitDiagram("|0>-H-.---x\n"
                                    "|0>-H-Z-H-x\n");*/
  auto c = parseCircuitDiagram("|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n"
                               "|0>HHH\n");
  // c.compile();
  if (v == 1) {
    cerr << "cpu\n";
    // c.compile();
    c.run(makeStateVector(13));
  } else if (v == 0) {
    cerr << "gpu\n";
    runCircuitOnGPU(c, makeStateVector(13));
  } else {
    return 0;
  }
  // sixteen.run(makeStateVector(16));
  // printStateVector(eight.run(makeStateVector(8)));

  //  auto res = runCircuitOnGPU(eight, makeStateVector(8));
}
