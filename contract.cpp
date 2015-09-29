/*
 *  contract.cpp a MEX-file: Contraction of two arrays for Matlab 5.3
 *
 *  Ville Bergholm 2002-2008-2015
 */

#include <new>
#include <stdlib.h>
#include "mex.h"

void newhandler()
{
  mexErrMsgTxt("Couldn't allocate memory using 'new'!");
  abort();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  std::set_new_handler(newhandler);

  long int i, j, k; // loop variables

  // Check for proper number of input and output arguments.
  if (nrhs != 4) mexErrMsgTxt("Input: array A, array B, index a, index b");

  if (nlhs > 1) mexErrMsgTxt("One output argument only!");

  // Check data type of input arguments.
  if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])))
    mexErrMsgTxt("Input arrays must be of type double.");

  if (!(mxIsDouble(prhs[2])) || !(mxIsDouble(prhs[3])))
    mexErrMsgTxt("Contraction indices must be of type double.");


  bool compa = mxIsComplex(prhs[0]);
  bool compb = mxIsComplex(prhs[1]);
  bool comp = compa || compb;

  int ndim_a = mxGetNumberOfDimensions(prhs[0]);
  int ndim_b = mxGetNumberOfDimensions(prhs[1]);


  const int *dim_a, *dim_b;
  dim_a = mxGetDimensions(prhs[0]);
  dim_b = mxGetDimensions(prhs[1]);


  // contr. indices are the truncated real parts of the first elements in inputs 3 and 4
  int ca = (int)(mxGetPr(prhs[2])[0]);
  int cb = (int)(mxGetPr(prhs[3])[0]);
  if (ca < 1 || ca > ndim_a)
    mexErrMsgTxt("Contracted index of array A out of bounds!");
  if (cb < 1 || cb > ndim_b)
    mexErrMsgTxt("Contracted index of array B out of bounds!");


  // indices into c style
  ca--; cb--;
  // dim[0] is the first dimension
  int lc = dim_a[ca];
  if (lc !=  dim_b[cb])
    mexErrMsgTxt("Contracted indices must be of equal length!");


  const double *par, *pai, *pbr, *pbi;  // data array pointers
  // get pointers to data
  par = mxGetPr(prhs[0]);
  pbr = mxGetPr(prhs[1]);
  if (comp == true) {
    // Q: what if only one of a,b is complex? one im. data pointer is NULL then.
    // A: only the im. data pointer of the complex one is ever dereferenced. no problemo.
    pai = mxGetPi(prhs[0]);
    pbi = mxGetPi(prhs[1]);
  }

  // output array dimensions
  int ndim_o = ndim_a + ndim_b - 2;

  // could also use malloc()
  int *dim_o = new int[ndim_o];

  for (i=0; i < ca; i++)        dim_o[i]   = dim_a[i];
  for (i=ca+1; i < ndim_a; i++) dim_o[i-1] = dim_a[i];
  for (i=0; i < cb; i++)        dim_o[i+ndim_a-1] = dim_b[i];
  for (i=cb+1; i < ndim_b; i++) dim_o[i+ndim_a-2] = dim_b[i];

  // index and counter tables
  int *idim_a = new int[ndim_a+1];
  int *idim_b = new int[ndim_b+1];
  int *cdim_a = new int[ndim_a+1];
  int *cdim_b = new int[ndim_b+1];
  idim_a[0] = idim_b[0] = 1;
  cdim_a[0] = cdim_b[0] = 0;
  for (i=1; i <= ndim_a; i++) {
    idim_a[i] = idim_a[i-1] * dim_a[i-1];
    cdim_a[i] = 0;
  }
  // here were the bugs
  for (i=1; i <= ndim_b; i++) {
    idim_b[i] = idim_b[i-1] * dim_b[i-1];
    cdim_b[i] = 0;
  }

  // we try to avoid slow multiplications and favor fast additions,
  // that's why we use these counters and index tables

  // create the output array
  mxArray *out = mxCreateNumericArray(ndim_o, dim_o, mxDOUBLE_CLASS, comp ? mxCOMPLEX : mxREAL);
  if (!out)
    mexErrMsgTxt("Could not create output array.");

  double *por, *poi;
  por = mxGetPr(out);
  if (comp == true) poi = mxGetPi(out);


  // and finally the actual contraction

  // resulting array size
  long int size=1;
  for (i = 0; i<ndim_o; i++) size *= dim_o[i];

  // mexPrintf("a");

  // the main contraction loop. i is the absolute index to the resulting array
  long int ia = 0, ib = 0;
  for (i = 0; i < size; i++) {
    long int tia = ia;
    long int tib = ib;

    for (k=0; k < lc; k++) {
      if (compa == true) { //complex data in a
	if (compb == true) { // both complex
	  por[i] += par[tia] * pbr[tib] - pai[tia] * pbi[tib];
	  poi[i] += par[tia] * pbi[tib] + pai[tia] * pbr[tib];
	} else { // a complex, b real
	  por[i] += par[tia] * pbr[tib];
	  poi[i] += pai[tia] * pbr[tib];
	}
      } else { // real data in a
	if (compb == true) { // a real, b complex
	  por[i] += par[tia] * pbr[tib];
	  poi[i] += par[tia] * pbi[tib];
	} else {  // plain real data
	por[i] += par[tia] * pbr[tib];
	}
      }
      tia += idim_a[ca];
      tib += idim_b[cb];
    }
    //mexPrintf("%d",i);
    //this keeps track of the current position in a and b
    for (j=0; ; j++) {
      if (j == ca) j++; //skip contracted index
      if (j >= ndim_a) { // go on to b
	for (j=0; ; j++) {
	  if (j == cb) j++; //skip contracted index
	  if (j >= ndim_b) break; // finally ends here...
	  if (++cdim_b[j] >= dim_b[j]) {
	    // move to next index
	    ib -= idim_b[j+1] - idim_b[j];
	    cdim_b[j] = 0;
	  } else {
	    ib += idim_b[j];
	    break; //infinite for
	  }
	}
	break;  // ...and here
      }
      if (++cdim_a[j] >= dim_a[j]) {
	// move to next index
	ia -= idim_a[j+1] - idim_a[j];
	cdim_a[j] = 0;
      } else {
	ia += idim_a[j];
	break; //infinite for
      }
    }
  }

  // mexPrintf("l\n");
  plhs[0] = out;
  delete [] dim_o;
  delete [] idim_a;
  delete [] idim_b;
  delete [] cdim_a;
  delete [] cdim_b;
  return;
}

