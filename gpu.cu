#include "cpu.cpp"
#include <cuda/std/complex>

using ResultVector = vector<StateVector>;

__device__ void cuMatrixVectorMultiply(size_t sLen, Complex *inV, Complex *outV,
                                       Complex *M) {
  size_t idy = (blockIdx.x * 1024) + threadIdx.x;
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

__device__ void cuTensorProduct(size_t sLen, Complex *M, Complex *N,
                                Complex *P) {
  int mRow = blockIdx.y * 32 + threadIdx.y;
  int mCol = blockIdx.x * 32 + threadIdx.x;
  if (mRow > sLen || mCol > sLen) {
    return;
  }

  // dimensions of P = len(M) * len(N)
  // since len(N) always = 2 & len(M) = sLen --> P = sLen * 2
  int pLen = sLen * 2;

  Complex mVal = M[(mRow * sLen) + mCol];
  /// TODO check that it doesnt overflow the thisVal part

  for (int col_i = 0; col_i < 2; col_i++) {
    for (int row_i = 0; row_i < 2; row_i++) {
      Complex nVal = N[(row_i * 2) + col_i];
      int P_row = (2 * mRow) + row_i;
      int P_col = (2 * mCol) + col_i;
      P[(P_row * pLen) + P_col] = Complex(mVal.x * nVal.x, mVal.y * nVal.y);
    }
  }
}

__global__ void QTensorKernel(size_t nMatrices, size_t matrixSideLength,
                              Complex *matrices, Complex *outMatrix,
                              Complex *scratchMatrix) {
  // outMatrix already contains matrices[0]
  for (size_t idx = 1; idx < nMatrices; idx++) {
    Complex *matrixBP = matrices + (idx * matrixSideLength * matrixSideLength);
    // Set the matrix side length to = sizeOutMatrix * 2
    //  M tensor N = P (where P is in our scratch/temp matrix)
    cuTensorProduct(matrixSideLength, outMatrix, matrixBP, scratchMatrix);
    __syncthreads();
    // copy the scratch/P matrix into the outMatrix and move to the next tensor
    //    memcpy(outMatrix, scratchMatrix, pow(matrixSideLength * 2, 2));
  }
}

void lowerMatrix(Matrix m, Complex *matrices, size_t baseOffset) {
  Complex *bp = matrices + baseOffset;
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

int ceilDiv(int a, int b) { return ceil((float)a / float(b)); }

Complex *gpu_toTransformationMatrix(
    pair<optional<ControlledGate *>, vector<Matrix>> pti) {
  vector<Matrix> series = pti.second;
  size_t matrixSideLength = series[0].size();
  size_t sizeOfMatrix = sizeof(Complex) * matrixSideLength * matrixSideLength;
  size_t nMatrices = series.size();
  size_t sizeOfMatrices = sizeOfMatrix * nMatrices;

  Complex *matrices = (Complex *)malloc(sizeOfMatrix * nMatrices);

  size_t pc = 0;
  for (auto mat : series) {
    lowerMatrix(mat, matrices, pc * matrixSideLength * matrixSideLength);
  }

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

  dim3 dimBlock(32, 32);
  dim3 dimGrid(ceilDiv(outMatrixSideLength, 32),
               ceilDiv(outMatrixSideLength, 32));
  QTensorKernel<<<dimGrid, dimBlock>>>(nMatrices, matrixSideLength, d_matrices,
                                       d_outMatrix, d_scratchMatrix);
  cudaMemcpy(outMatrix, d_outMatrix, totalOutMatrixSize,
             cudaMemcpyDeviceToHost);
  return outMatrix;
}

ResultVector runCircuitOnGPU(Circuit C, StateVector SV, int onlyTensor,
                             size_t iters) {
  size_t matrixSideLength = pow(2, C.nQubits);

  size_t sizeOfComplex = sizeof(Complex); // 2 * sizeof(float);
  size_t sizeOfMatrix = sizeOfComplex * matrixSideLength * matrixSideLength;
  size_t sizeOfMatrices = sizeOfMatrix * C.operationsCount();

  Complex *matrices = (Complex *)malloc(sizeOfMatrix * C.operationsCount());

  size_t pc = 0;

  cerr << "tensoring\n";

  for (auto slice : C.program) {
    if (auto *op = dynamic_cast<GateTimeSlice *>(slice)) {

      Complex *transformationMatrix =
          gpu_toTransformationMatrix(op->toPreTensorInfo());

      if (!onlyTensor) {
        copyRMatrix(transformationMatrix,
                    matrices + (pc * matrixSideLength * matrixSideLength),
                    matrixSideLength);
      }
      pc++;
    }
  }
  if (onlyTensor) {
    return {};
  }
  cerr << "running\n";

  size_t nOps = C.program.size();

  // matrices now contains the floats
  Complex *d_matrices;
  cudaMalloc((void **)&d_matrices, sizeOfMatrices);
  cudaMemcpy(d_matrices, matrices, sizeOfMatrices, cudaMemcpyHostToDevice);

  size_t totalOutVectorSize = nOps * matrixSideLength * sizeof(Complex);
  Complex *d_outVectors;
  cudaMalloc((void **)&d_outVectors, totalOutVectorSize);

  for (size_t idx = 0; idx < iters; idx++) {

    Complex *h_inVector = lowerStateVector(SV);
    Complex *d_inVector;
    cudaMalloc((void **)&d_inVector, sizeof(Complex) * matrixSideLength);
    cudaMemcpy(d_inVector, h_inVector, sizeof(Complex) * matrixSideLength,
               cudaMemcpyHostToDevice);

    dim3 dimGrid(ceilDiv(matrixSideLength, 1024));
    dim3 dimBlock(1024);

    QProgramKernel<<<dimGrid, dimBlock>>>(nOps, matrixSideLength, d_matrices,
                                          d_inVector, d_outVectors);

    Complex *h_outVectors = (Complex *)malloc(totalOutVectorSize);

    cudaMemcpy(h_outVectors, d_outVectors, totalOutVectorSize,
               cudaMemcpyDeviceToHost);
  }
  return {};
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "usage: ./q-emu [0|1] [0|1]. First 1 is use GPU, second 1 is only "
            "tensor\n";
    return 0;
  }

  int useGPU = atoi(argv[1]);
  int onlyTensor = atoi(argv[2]);
  Circuit c = groversCircuit_NR(10, "1101010101");

  if (useGPU == 0) {
    cerr << "cpu\n";
    cerr << "tensoring\n";
    c.preTensor(); // this changes the stored program
    if (onlyTensor) {
      return 0;
    }
    cerr << "running\n";
    for (size_t idx = 0; idx < 100; idx++) {
      c.run(makeStateVector(13));
    }
  } else if (useGPU == 1) {
    cerr << "gpu\n";
    runCircuitOnGPU(c, makeStateVector(10), onlyTensor, 100);
  } else {
    cerr << "bad arg. run with no args for help\n";
  }
}
