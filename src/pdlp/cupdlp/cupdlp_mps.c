#include "cupdlp_mps.h"

#include <math.h>

#include "cupdlp_cs.h"

/* Implement hash mapping from

 https://stackoverflow.com/questions/4384359/quick-way-to-implement-dictionary-in-c

 */

struct cupdlp_hash_internal {
  struct cupdlp_hash_internal *next;
  char key[128];
  unsigned int val;
};

typedef struct cupdlp_hash_internal cupdlp_hash;

typedef struct {
  int nMaxElem;
  int nElem;

  cupdlp_hash **hash;

} cupdlp_dict;

static unsigned int hash(char *str, int nHashElem) {
  unsigned int iHash = 0;

  for (; *str != '\0'; ++str) {
    iHash += *str + 31 * iHash;
    if (iHash > 16777216) {
      iHash = iHash % nHashElem;
    }
  }

  return iHash % nHashElem;
}

static cupdlp_hash *get(cupdlp_dict *dict, char *key) {
  unsigned int iHash = hash(key, dict->nMaxElem);
  cupdlp_hash *elem = dict->hash[iHash];

  for (; elem != NULL; elem = elem->next) {
    if (strcmp(key, elem->key) == 0) {
      return elem;
    }
  }

  return NULL;
}

static int freehash(cupdlp_hash *hash, int nfreed) {
  if (hash->next) {
    nfreed = freehash(hash->next, nfreed);
  }

  cupdlp_free(hash);
  return nfreed + 1;
}

static int rowIdxsplit(cupdlp_int m, cupdlp_int n, cupdlp_int *Ap,
                       cupdlp_int *Ai, double *Ax, cupdlp_int *rowIndex,
                       cupdlp_int **pBp, cupdlp_int **pBi, double **pBx,
                       cupdlp_int **pCp, cupdlp_int **pCi, double **pCx,
                       double *b) {
  /*
     Split an csc matrix into two according to the value of rowIndex:
     Rows corresponding to 0 in rowIndex will be put in matrix B
     Rows corresponding to non-zero in rowIndex will be put in matrix C
  */
  cupdlp_int retcode = RETCODE_OK;

  int nBrow = 0;
  int nCrow = 0;

  for (int i = 0; i < m; ++i) {
    if (rowIndex[i]) {
      nCrow += 1;
    } else {
      nBrow += 1;
    }
  }

  /* We are mapping the rows to a smaller set from 1 to # of rows*/
  int *BrowMap = NULL;
  int *CrowMap = NULL;
  double *bRow = NULL;

  CUPDLP_INIT(BrowMap, m);
  CUPDLP_INIT(CrowMap, m);
  CUPDLP_INIT(bRow, m);

  int iBrow = 0;
  int iCrow = 0;
  for (int i = 0; i < m; ++i) {
    if (rowIndex[i]) {
      CrowMap[i] = iCrow;
      bRow[nBrow + iCrow] = b[i];
      iCrow += 1;
    } else {
      BrowMap[i] = iBrow;
      bRow[iBrow] = b[i];
      iBrow += 1;
    }
  }

  int nBnz = 0;
  int nCnz = 0;

  /* First iterate through the matrix to get the number of nonzeros*/
  for (int i = 0; i < n; ++i) {
    for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
      int iRow = Ai[j];
      if (rowIndex[iRow]) {
        nCnz += 1;
      } else {
        nBnz += 1;
      }
    }
  }

  assert(nBnz + nCnz == Ap[n]);

  /* Allocate memory for B and C */
  cupdlp_int *Bp = NULL;
  cupdlp_int *Bi = NULL;
  double *Bx = NULL;

  cupdlp_int *Cp = NULL;
  cupdlp_int *Ci = NULL;
  double *Cx = NULL;

  /* We allocate one more unit of memory in case there is no B or C */
  CUPDLP_INIT(Bp, n + 1);
  CUPDLP_INIT(Bi, nBnz + 1);
  CUPDLP_INIT(Bx, nBnz + 1);

  CUPDLP_INIT(Cp, n + 1);
  CUPDLP_INIT(Ci, nCnz + 1);
  CUPDLP_INIT(Cx, nCnz + 1);

  int iBnz = 0;
  int iCnz = 0;

  /* Iterate again to fill in the data */
  for (int i = 0; i < n; ++i) {
    for (int j = Ap[i]; j < Ap[i + 1]; ++j) {
      int iRow = Ai[j];

      if (rowIndex[iRow]) {
        Ci[iCnz] = CrowMap[iRow];
        Cx[iCnz] = Ax[j];
        iCnz += 1;
      } else {
        Bi[iBnz] = BrowMap[iRow];
        Bx[iBnz] = Ax[j];
        iBnz += 1;
      }
    }

    Bp[i + 1] = iBnz;
    Cp[i + 1] = iCnz;
  }

  *pBp = Bp;
  *pBi = Bi;
  *pBx = Bx;
  *pCp = Cp;
  *pCi = Ci;
  *pCx = Cx;

  cupdlp_copy(b, bRow, double, m);

exit_cleanup:

  if (retcode != RETCODE_OK) {
    if (Bp) {
      cupdlp_free(Bp);
    }

    if (Bi) {
      cupdlp_free(Bi);
    }

    if (Bx) {
      cupdlp_free(Bx);
    }

    if (Cp) {
      cupdlp_free(Cp);
    }

    if (Ci) {
      cupdlp_free(Ci);
    }

    if (Cx) {
      cupdlp_free(Cx);
    }
  }

  if (bRow) {
    cupdlp_free(bRow);
  }

  if (BrowMap) {
    cupdlp_free(BrowMap);
  }

  if (CrowMap) {
    cupdlp_free(CrowMap);
  }

  return retcode;
}

extern cupdlp_int cupdlpDictCreate(cupdlp_dict **pDict, int nMaxElem) {
  cupdlp_int retcode = RETCODE_OK;

  if (!pDict) {
    return retcode;
  }

  cupdlp_dict *dict = NULL;
  CUPDLP_INIT(dict, 1);

  /* Balance load of access */
  dict->nMaxElem = (int)(nMaxElem / 0.700);
  dict->nElem = 0;

  CUPDLP_INIT(dict->hash, dict->nMaxElem);

  *pDict = dict;

exit_cleanup:

  return retcode;
}

extern cupdlp_int cupdlpDictAddElem(cupdlp_dict *dict, char *key, int val) {
  cupdlp_int retcode = RETCODE_OK;

  if (dict->nElem >= dict->nMaxElem) {
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  cupdlp_hash *elem = get(dict, key);

  if (!elem) {
    CUPDLP_INIT(elem, 1);

    unsigned int hashval = hash(key, dict->nMaxElem);

    elem->next = dict->hash[hashval];
    elem->val = val;
    cupdlp_copy(elem->key, key, char, strlen(key));
    dict->hash[hashval] = elem;
  } else {
    /* Two keys are the same. Now allowed in the LP context */
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  dict->nElem += 1;

exit_cleanup:

  return retcode;
}

extern unsigned int cupdlpDictMapElem(cupdlp_dict *dict, char *key) {
  cupdlp_hash *hash = get(dict, key);

  if (!hash) {
    return -1;
  } else {
    return hash->val;
  }
}

extern void cupdlpDictClear(cupdlp_dict *dict) {
  if (!dict) {
    return;
  }

  int iHashElem = 0;
  for (int i = 0; i < dict->nMaxElem; ++i) {
    if (dict->hash[i]) {
      int nFreedElem = freehash(dict->hash[i], 0);
      iHashElem += nFreedElem;
    }
  }

  assert(dict->nElem == iHashElem);
  cupdlp_free(dict->hash);
  cupdlp_zero(dict, cupdlp_dict, 1);

  return;
}

extern void cupdlpDictDestroy(cupdlp_dict **pDict) {
  if (!pDict) {
    return;
  }

  cupdlpDictClear(*pDict);
  cupdlp_free(*pDict);

  return;
}

/* LP-related
 I used
 https://www.ibm.com/docs/en/icos/12.8.0.0?topic=standard-records-in-mps-format
 for the mps standard format
 */
#define INDICATOR_NAME ("NAME")
#define INDICATOR_ROWS ("ROWS")
#define INDICATOR_COLS ("COLUMNS")
#define INDICATOR_RHS ("RHS")
#define INDICATOR_RANGE ("RANGES")
#define INDICATOR_BOUNDS ("BOUNDS")
#define INDICATOR_END ("ENDATA")

#define CONSTR_SENSE_OBJ ('N')
#define CONSTR_SENSE_EQ ('E')
#define CONSTR_SENSE_LEQ ('L')
#define CONSTR_SENSE_GEQ ('G')

#define BOUND_SENSE_UP ("UP")
#define BOUND_SENSE_LOW ("LO")

#define LINE_BUFFER (512)
#define str_begin_with(pre, str) (strncmp((pre), (str), strlen((pre))) == 0)

/* Implement an LP mps file reader
   Ignore all comments and names, only serving purpose of extracting LP data
 */
cupdlp_int cupdlpMpsRead(char *fname, char *name, int *pnRow, int *pnEqRow,
                         int *pnInEqRow, int *pnCol, int *pnElem,
                         int **pfullMatBeg, int **pfullMatIdx,
                         double **pfullMatElem, int **peqMatBeg,
                         int **peqMatIdx, double **peqMatElem,
                         int **pIneqMatBeg, int **pIneqMatIdx,
                         double **pIneqMatElem, double **prowRHS,
                         double **pcolObj, int *pnColUb, int **pcolUbIdx,
                         double **pcolUbElem) {
  cupdlp_int retcode = RETCODE_OK;

  FILE *mps = NULL;

  int nLine = 0;
  char probName[LINE_BUFFER] = "?";
  int nRow = 0;
  int nEqRow = 0;
  int nInEqRow = 0;
  int nCol = 0;
  int nElem = 0;
  int nBound = 0;

  /* LP data */
  int *eqMatBeg = NULL;
  int *eqMatIdx = NULL;
  double *eqMatElem = NULL;

  int *inEqMatBeg = NULL;
  int *inEqMatIdx = NULL;
  double *inEqMatElem = NULL;

  int *colUbIdx = NULL;
  double *colUbElem = NULL;

  double *rowRHS = NULL;
  double *colObj = NULL;

  cupdlp_dcs *colMat = NULL;
  cupdlp_dcs *cscMat = NULL;

  /* We use -1 to denote >=, 0 to denote ==, and 1 to denote <= */
  int *rowSenses = NULL;

  /* Variable and constraint hash */
  cupdlp_dict *rowHash = NULL;
  cupdlp_dict *colHash = NULL;

  printf("Reading specialized standard form mps %s \n", fname);
  mps = fopen(fname, "r");

  if (!mps) {
    printf("Failed to open file \"%s\". \n", fname);
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  /* Get number of constraints and variables */
  char thisLine[LINE_BUFFER] = "*";
  fgets(thisLine, LINE_BUFFER, mps);
  nLine += 1;

  if (!str_begin_with(INDICATOR_NAME, thisLine)) {
    printf("Line [%d] contains no NAME argument. \n", nLine);
    printf("Line content: %s \n", thisLine);
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  /* Get problem name */
  strncpy(probName, thisLine + 5, LINE_BUFFER - 5);
  /* Remove \n */
  probName[strcspn(probName, "\n")] = '\0';

  /* Moving on to ROW */
  fgets(thisLine, LINE_BUFFER, mps);
  nLine += 1;

  /* First count number of rows and columns */
  int nLineBefore = nLine;
  if (!str_begin_with(INDICATOR_ROWS, thisLine)) {
    printf("Line [%d] contains no %s argument. \n", nLine, INDICATOR_ROWS);
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  for (nRow = 0; !feof(mps); ++nRow) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;
    if (str_begin_with(INDICATOR_COLS, thisLine)) {
      break;
    }
    /* Till here nRow contains the objective row */
  }

  /* Go on to columns */
  int nget = 0;
  int nNz = 0;
  char rname[128] = "*";
  char cname[128] = "*";
  char cname2[128] = "*";
  char objname[128] = "*";
  double dElem = 0.0;
  for (nCol = 0; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_RHS, thisLine)) {
      break;
    }

    nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);

    if (nget != 3) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (strcmp(cname, cname2) != 0) {
      nCol += 1;
      strcpy(cname2, cname);
    }

    nNz += 1;
  }

  /* Move on to the upperbounds */
  for (; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_BOUNDS, thisLine)) {
      break;
    }
  }

  char bound[4] = "*";
  for (nBound = 0; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_END, thisLine)) {
      break;
    }

    nget = sscanf(thisLine, "%s %s %s %lg", bound, rname, cname, &dElem);

    if (nget != 4) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (strcmp(bound, BOUND_SENSE_UP) != 0 && dElem != 0.0) {
      printf("Non 'UP' sense detected. \n");
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    nBound += 1;
  }

  /* Till now the number of rows (including c) and columns are both known */
  /* Return to the start of file */
  fseek(mps, 0, SEEK_SET);
  for (nLine = 0; nLine < nLineBefore; ++nLine) {
    fgets(thisLine, LINE_BUFFER, mps);
  }

  /* Subtract the objective row off */
  nRow -= 1;

  /* Build up Hash mapping for rows and columns */
  CUPDLP_CALL(cupdlpDictCreate(&rowHash, nRow));
  CUPDLP_CALL(cupdlpDictCreate(&colHash, nCol));

  /* Prepare matrix data */
  colMat = cupdlp_dcs_spalloc(nRow, nCol, nNz, 1, 1);

  if (!colMat) {
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  /* Prepare vector data */
  CUPDLP_INIT(rowRHS, nRow);
  CUPDLP_INIT(colObj, nCol);
  CUPDLP_INIT(rowSenses, nRow);
  CUPDLP_INIT(colUbIdx, nBound + 1);
  CUPDLP_INIT(colUbElem, nBound + 1);

  /* Build up hash and go through senses */
  int iRhs = 0;
  char sense = '\0';

  for (iRhs = 0; !feof(mps); ++iRhs) {
    fgets(thisLine, LINE_BUFFER, mps);
    nget = sscanf(thisLine, " %c %s", &sense, rname);
    nLine += 1;

    if (str_begin_with(INDICATOR_COLS, thisLine)) {
      break;
    }

    if (nget != 2) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (sense == CONSTR_SENSE_OBJ) {
      /* There is a row of objective */
      strcpy(objname, rname);
      iRhs -= 1;
      continue;
    } else {
      CUPDLP_CALL(cupdlpDictAddElem(rowHash, rname, iRhs));
      if (sense == CONSTR_SENSE_GEQ) {
        rowSenses[iRhs] = -1;
        nInEqRow += 1;
      } else if (sense == CONSTR_SENSE_LEQ) {
        rowSenses[iRhs] = 1;
        nInEqRow += 1;
      } else {
        nEqRow += 1;
      }
    }
  }

  assert(iRhs == nRow && nRow == nEqRow + nInEqRow);

  /* Collect variable data */
  int iCol = 0;
  cname2[0] = '\0';

  for (iCol = 0; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_RHS, thisLine)) {
      break;
    }

    nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);

    if (nget != 3) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (strcmp(cname, cname2) != 0) {
      CUPDLP_CALL(cupdlpDictAddElem(colHash, cname, iCol));
      iCol += 1;
      strcpy(cname2, cname);
    }

    /* Objective vector */
    if (strcmp(rname, objname) == 0) {
      int iCol = cupdlpDictMapElem(colHash, cname);
      colObj[iCol] = dElem;
    } else {
      int iCol = cupdlpDictMapElem(colHash, cname);
      int iRow = cupdlpDictMapElem(rowHash, rname);

      assert(iCol >= 0 && iRow >= 0);

      /* Revert the sense for >= constraint */
      if (rowSenses[iRow] == -1) {
        dElem = -dElem;
      }

      if (cupdlp_dcs_entry(colMat, iRow, iCol, dElem) != 1) {
        retcode = RETCODE_FAILED;
        goto exit_cleanup;
      }
    }
  }

  /* Collect RHS */
  for (; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_BOUNDS, thisLine)) {
      break;
    }

    if (str_begin_with(INDICATOR_END, thisLine)) {
      break;
    }

    nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);

    if (nget != 3) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    /* If found obj shift */
    if (strcmp(rname, objname) == 0) {
      printf("Shifting model objective by %5.3e \n", -dElem);
      continue;
    }

    int iRow = cupdlpDictMapElem(rowHash, rname);

    if (rowSenses[iRow] == -1) {
      rowRHS[iRow] = -dElem;
    } else {
      rowRHS[iRow] = dElem;
    }
  }

  /* Collect bounds */
  int iBound = 0;
  for (; !feof(mps);) {
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;

    if (str_begin_with(INDICATOR_END, thisLine)) {
      break;
    }

    nget = sscanf(thisLine, "%s %s %s %lg", bound, rname, cname, &dElem);

    if (nget != 4) {
      printf("Error at line %d. \n", nLine);
      retcode = RETCODE_FAILED;
      goto exit_cleanup;
    }

    if (dElem < INFINITY && dElem > -INFINITY) {
      int iCol = cupdlpDictMapElem(colHash, cname);
      colUbIdx[iBound] = iCol;
      colUbElem[iBound] = dElem;
      iBound += 1;
    } else {
      printf("Warning: ignored upperbound %5.2e. \n", dElem);
    }
  }

  assert(iBound == nBound);

  /* Finished */
  cscMat = cupdlp_dcs_compress(colMat);

  if (!cscMat) {
    retcode = RETCODE_FAILED;
    goto exit_cleanup;
  }

  /* Get the results done */
  if (name) {
    strcpy(name, probName);
  }

  /* Split rows of inequality and equality */
  CUPDLP_CALL(rowIdxsplit(nRow, nCol, cscMat->p, cscMat->i, cscMat->x,
                          rowSenses, &eqMatBeg, &eqMatIdx, &eqMatElem,
                          &inEqMatBeg, &inEqMatIdx, &inEqMatElem, rowRHS));

  nElem = cscMat->p[nCol];

#if 0
    cupdlp_dcs A, B;
    A.p = eqMatBeg;
    A.i = eqMatIdx;
    A.x = eqMatElem;
    A.nz = -1;
    A.m = nEqRow;
    A.n = nCol;
    cupdlp_dcs_print(&A, 0);
    
    B.p = inEqMatBeg;
    B.i = inEqMatIdx;
    B.x = inEqMatElem;
    B.nz = -1;
    B.m = nInEqRow;
    B.n = nCol;
    cupdlp_dcs_print(&B, 0);
#endif

  *pnRow = nRow;
  *pnEqRow = nEqRow;
  *pnInEqRow = nInEqRow;
  *pnColUb = nBound;
  *pnCol = nCol;
  *pnElem = nElem;
  *prowRHS = rowRHS;
  *pcolObj = colObj;
  *pcolUbIdx = colUbIdx;
  *pcolUbElem = colUbElem;
  *peqMatBeg = eqMatBeg;
  *peqMatIdx = eqMatIdx;
  *peqMatElem = eqMatElem;
  *pIneqMatBeg = inEqMatBeg;
  *pIneqMatIdx = inEqMatIdx;
  *pIneqMatElem = inEqMatElem;
  *pfullMatBeg = cscMat->p;
  *pfullMatIdx = cscMat->i;
  *pfullMatElem = cscMat->x;

exit_cleanup:

  if (retcode != RETCODE_OK) {
    if (eqMatBeg) {
      cupdlp_free(eqMatBeg);
    }

    if (eqMatIdx) {
      cupdlp_free(eqMatIdx);
    }

    if (eqMatElem) {
      cupdlp_free(eqMatElem);
    }

    if (inEqMatBeg) {
      cupdlp_free(inEqMatBeg);
    }

    if (inEqMatIdx) {
      cupdlp_free(inEqMatIdx);
    }

    if (inEqMatElem) {
      cupdlp_free(inEqMatElem);
    }

    if (colUbIdx) {
      cupdlp_free(colUbIdx);
    }

    if (colUbElem) {
      cupdlp_free(colUbElem);
    }

    if (rowRHS) {
      cupdlp_free(rowRHS);
    }

    if (colObj) {
      cupdlp_free(colObj);
    }
  }

  cupdlpDictDestroy(&rowHash);
  cupdlpDictDestroy(&colHash);

  cupdlp_free(rowSenses);

  cupdlp_dcs_spfree(colMat);
  // cupdlp_dcs_spfree(cscMat);

  if (mps) {
    fclose(mps);
  }

  return retcode;
}
